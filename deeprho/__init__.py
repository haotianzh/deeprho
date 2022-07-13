from .version import __version__
import importlib
import sys


class LazyLoader:
    """
        LazyLoader for loading modules when they are called at the first time.
    """
    def __init__(self, lib_name):
        self.lib_name = lib_name
        self._mod = None
    def __getattr__(self, name):
        if self._mod is None:
            self._mod = importlib.import_module(self.lib_name)
        return getattr(self._mod, name)


def lazy_import(name):
    """
        Implementation in importlib docs.
        Issue is importlib.util.find_spec() will execute __init__.py once.
    """
    spec = importlib.util.find_spec(name)
    loader = importlib.util.LazyLoader(spec.loader)
    spec.loader = loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module