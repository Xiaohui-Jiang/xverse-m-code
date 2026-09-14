"""Strict, read-only integration with an external xVERSE checkout."""
import hashlib
import importlib
import sys
from pathlib import Path

import torch


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def imports(repo):
    repo = Path(repo).resolve()
    sys.path.insert(0, str(repo))
    model_module = importlib.import_module("main.utils_model")
    if Path(model_module.__file__).resolve() != repo / "main/utils_model.py":
        raise RuntimeError("A different main.utils_model is already imported")
    return model_module, importlib.import_module("main.data")


def load_frozen(repo, checkpoint_path, device="cpu"):
    """Rebuild the original model, with no replacement sample embeddings."""
    module, _ = imports(repo)
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    saved = checkpoint.get("training_args", {})
    required = ("normalize_z_bio", "recon_distribution", "nb_dispersion_mode",
                "expr_encoder_depth", "mask_encoder_depth", "decoder_depth")
    missing = [key for key in required if key not in saved]
    if missing:
        raise ValueError(f"Checkpoint lacks required architecture metadata: {missing}")
    state = {k.removeprefix("module."): v for k, v in checkpoint["model_state_dict"].items()}
    sample = state.get("sample_emb.weight")
    text = state.get("celltype_text_embeddings")
    panel = state.get("panel_contrast_head.3.weight")
    kwargs = {key: saved[key] for key in required}
    for key, default in (("dropout", 0.1), ("dispersion_rank", 16),
                         ("dispersion_min", 1e-4), ("dispersion_max", 1e4)):
        kwargs[key] = saved.get(key, default)
    model = module.XVerseModel(
        num_samples=None if sample is None else sample.shape[0],
        batch_emb_dim=64 if sample is None else sample.shape[1],
        hidden_dim=state["bio_encoder.latent_head.weight"].shape[0],
        total_gene=state["decoder.out.weight"].shape[0],
        expr_hidden_dim=state["bio_encoder.expr_encoder.block.0.weight"].shape[0],
        mask_hidden_dim=state["bio_encoder.mask_encoder.block.0.weight"].shape[0],
        dec_hidden_dim=state["decoder.fc1.weight"].shape[0],
        total_cell_type=text.shape[0] if text is not None else state["cell_type_head.4.weight"].shape[0],
        celltype_text_embeddings=text,
        celltype_text_temperature=saved.get("celltype_text_temp", 0.2),
        panel_contrast_projection_dim=0 if panel is None else panel.shape[0],
        **kwargs,
    )
    model.load_state_dict(state, strict=True)
    model.requires_grad_(False).eval().to(device)
    return model


@torch.inference_mode()
def encode(model, counts, observed_mask):
    raw, _ = model.bio_encoder(counts, observed_mask=observed_mask, apply_mask_aug=False)
    return model.embedding_z(raw)


@torch.inference_mode()
def decode(model, z):
    """Use the frozen biological decoder and the actual cell-dependent likelihood."""
    mu = model.decode(z)
    if model.recon_distribution == "poisson":
        counts = torch.poisson(mu)
        mean = mu
    else:
        theta = model.predict_dispersion(z).expand_as(mu)
        # torch NB has mean total_count * exp(logits).
        distribution = torch.distributions.NegativeBinomial(
            total_count=theta, logits=mu.clamp_min(1e-12).log() - theta.log())
        counts = distribution.sample()
        mean = mu
        if model.recon_distribution == "zinb":
            dropout = torch.sigmoid(model.zero_inflation_head(z).float().clamp(-20, 20))
            counts = counts * (torch.rand_like(dropout) >= dropout)
            mean = mu * (1 - dropout)
    if not torch.isfinite(mean).all() or not torch.isfinite(counts).all():
        raise FloatingPointError("Non-finite decoded values")
    return mean, counts
