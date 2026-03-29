"""src/lmn_var/__main__.py"""
from .estimation import LMNVAR

if __name__ == '__main__':
    lmn = LMNVAR()
    lmn.run_all()
