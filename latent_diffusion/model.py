"""Cosine-schedule v-prediction diffusion with a time-modulated residual MLP."""
import math

import torch
from torch import nn
import torch.nn.functional as F


class Block(nn.Module):
    def __init__(self, width):
        super().__init__()
        self.norm = nn.LayerNorm(width, elementwise_affine=False)
        self.modulate = nn.Linear(width, 2 * width)
        self.net = nn.Sequential(nn.Linear(width, 2 * width), nn.SiLU(), nn.Linear(2 * width, width))

    def forward(self, x, time):
        scale, shift = self.modulate(F.silu(time)).chunk(2, -1)
        return (x + self.net(self.norm(x) * (1 + scale) + shift)) / math.sqrt(2)


class Denoiser(nn.Module):
    def __init__(self, dim, width=512, depth=6):
        super().__init__()
        if dim < 1 or width < 4 or width % 2 or depth < 1:
            raise ValueError("Require positive dim/depth and an even width >= 4")
        self.register_buffer("frequencies", torch.exp(-math.log(10000) * torch.arange(width // 2) / (width // 2 - 1)))
        self.time = nn.Sequential(nn.Linear(width, width), nn.SiLU(), nn.Linear(width, width))
        self.input = nn.Linear(dim, width)
        self.blocks = nn.ModuleList(Block(width) for _ in range(depth))
        self.output = nn.Sequential(nn.LayerNorm(width), nn.SiLU(), nn.Linear(width, dim))
        nn.init.zeros_(self.output[-1].weight)
        nn.init.zeros_(self.output[-1].bias)

    def forward(self, z, t):
        angles = t[:, None].float() * self.frequencies[None]
        time = self.time(torch.cat((angles.sin(), angles.cos()), -1))
        h = self.input(z)
        for block in self.blocks:
            h = block(h, time)
        return self.output(h)


class Diffusion:
    def __init__(self, steps=1000, device="cpu"):
        if steps < 2:
            raise ValueError("At least two diffusion steps are required")
        t = torch.linspace(0, 1, steps + 1, dtype=torch.float64)
        cumulative = torch.cos((t + 0.008) / 1.008 * math.pi / 2).square()
        cumulative = cumulative / cumulative[0]
        betas = (1 - cumulative[1:] / cumulative[:-1]).clamp(max=0.999)
        self.alpha = torch.cumprod(1 - betas, 0).float().to(device)
        self.steps = steps

    def coefficients(self, t):
        alpha = self.alpha[t, None]
        return alpha.sqrt(), (1 - alpha).sqrt()

    def corrupt(self, x, t, noise):
        a, s = self.coefficients(t)
        return a * x + s * noise, a * noise - s * x

    def recover(self, noisy, velocity, t):
        a, s = self.coefficients(t)
        return a * noisy - s * velocity, s * noisy + a * velocity

    def loss(self, model, x, generator=None):
        t = torch.randint(self.steps, (len(x),), device=x.device, generator=generator)
        noise = torch.randn(x.shape, device=x.device, generator=generator)
        noisy, target = self.corrupt(x, t, noise)
        return F.mse_loss(model(noisy, t), target)

    @torch.inference_mode()
    def sample(self, model, n, dim, sampling_steps=100, generator=None):
        """Deterministic DDIM updates; randomness comes from the initial noise."""
        if not 2 <= sampling_steps <= self.steps:
            raise ValueError("sampling_steps must be between 2 and diffusion steps")
        grid = torch.linspace(self.steps - 1, 0, sampling_steps).round().long().tolist()
        x = torch.randn((n, dim), device=self.alpha.device, generator=generator)
        for i, time in enumerate(grid):
            t = torch.full((n,), time, device=x.device, dtype=torch.long)
            clean, noise = self.recover(x, model(x, t), t)
            if i + 1 == len(grid):
                x = clean
            else:
                alpha = self.alpha[grid[i + 1]]
                x = alpha.sqrt() * clean + (1 - alpha).sqrt() * noise
        if not torch.isfinite(x).all():
            raise FloatingPointError("Non-finite diffusion sample")
        return x
