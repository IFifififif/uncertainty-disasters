"""Runtime configuration helpers for MODEL entrypoints."""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass(frozen=True)
class ModelRuntimeConfig:
    simplified: bool
    do_estimation: bool
    irf_t: int
    irf_n_sims: int


def _env_flag(name: str, default: str = "0") -> bool:
    return os.getenv(name, default) in {"1", "true", "True"}


def load_model_runtime_config_from_env() -> ModelRuntimeConfig:
    """Load MODEL runtime settings from environment variables."""
    return ModelRuntimeConfig(
        simplified=_env_flag("MODEL_SIMPLIFIED", "0"),
        do_estimation=_env_flag("MODEL_DO_ESTIMATION", "0"),
        irf_t=int(os.getenv("MODEL_IRF_T", "40")),
        irf_n_sims=int(os.getenv("MODEL_IRF_N_SIMS", "100")),
    )

