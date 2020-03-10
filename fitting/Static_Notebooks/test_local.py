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
    _exec_notebook('Fitting-BlackHole.ipynb')
    _exec_notebook('Fitting-Halo.ipynb')
    _exec_notebook('Fitting-Bulge.ipynb')
    _exec_notebook('Fitting-Disk.ipynb')
    _exec_notebook('FittingB-D.ipynb')
    #_exec_notebook('FittingB-H.ipynb') #This one will not enter testing because it takes too long.
    #Iff the bulge process is sped up then it can be added. 
    #It does, however, still need to finish a complete run so that runtime is up-to-date.
    #
    #_exec_notebook('FittingBH-B.ipynb) #Still waiting to run this; bulge file needs to be available.
    _exec_notebook('FittingBH-D.ipynb')
    _exec_notebook('NGC5533_Noordermeer_plot.ipynb')
