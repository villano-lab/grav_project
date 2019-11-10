import subprocess
import tempfile
import os

def _exec_notebook(path):
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=10000",
                "--ExecutePreprocessor.kernel_name=python3",
                "--output", fout.name, path]
        subprocess.check_call(args)


def test():
    print('Testing Jupyter notebooks...')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-BlackHole.ipynb')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-Halo.ipynb')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-Bulge.ipynb')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-Disk.ipynb')
