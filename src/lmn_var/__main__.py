"""src/lmn_var/__main__.py"""
import os
from .estimation import LMNVAR

if __name__ == '__main__':
    lmn = LMNVAR()
    n_draws = int(os.getenv("LMN_N_DRAWS", "1500000"))
    lmn.run_all(n_draws=n_draws)
