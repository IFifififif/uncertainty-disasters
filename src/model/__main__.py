"""src/model/__main__.py"""
from .solve import MicroMacroModel

if __name__ == '__main__':
    model = MicroMacroModel()
    model.run_all()
