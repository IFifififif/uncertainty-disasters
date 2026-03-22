"""src/iv_var/__main__.py"""
from .estimation import IVVAR

if __name__ == '__main__':
    ivvar = IVVAR()
    ivvar.run_all()
