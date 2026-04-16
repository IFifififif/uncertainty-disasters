"""src/model/__main__.py"""
import os
from .solve import MicroMacroModel

if __name__ == '__main__':
    simplified = os.getenv("MODEL_SIMPLIFIED", "1") not in {"0", "false", "False"}
    do_estimation = os.getenv("MODEL_DO_ESTIMATION", "0") in {"1", "true", "True"}
    irf_t = int(os.getenv("MODEL_IRF_T", "40"))
    irf_sims = int(os.getenv("MODEL_IRF_N_SIMS", "100"))
    model = MicroMacroModel(simplified=simplified)
    model.run_all(do_estimation=do_estimation, irf_T=irf_t, irf_n_sims=irf_sims)
