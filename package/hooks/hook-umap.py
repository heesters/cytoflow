from PyInstaller.utils.hooks import copy_metadata

# umap and pynndescent both look up their own version with importlib.metadata
# at import time.
datas = copy_metadata('umap-learn') + copy_metadata('pynndescent')

hiddenimports = ['pynndescent']
