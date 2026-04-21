"""src/model/__main__.py"""
from .solve import MicroMacroModel
from .runtime import load_model_runtime_config_from_env

if __name__ == '__main__':
    cfg = load_model_runtime_config_from_env()
    model = MicroMacroModel(simplified=cfg.simplified)
    model.run_all(
        do_estimation=cfg.do_estimation,
        irf_T=cfg.irf_t,
        irf_n_sims=cfg.irf_n_sims,
    )
