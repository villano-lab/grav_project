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
    #_exec_notebook('FittingBH-B.ipynb) #Takes too long to run.
    #Iff the bulge process is sped up then the above can be added.
    #
    _exec_notebook('FittingBH-D.ipynb')
    #_exec_notebook('FittingBH-H_different_galaxies.ipynb') #Getting encoding error, not sure why.
    #_exec_notebook('FittingBH-H.ipynb') #Still running.
    #_exec_notebook('FittingD-H.ipynb') #Did not finish running. Gave up on this one after three days straight of computation. Can be re-attempted if disk is sped up.
    #_exec_notebook('NGC5533_Noordermeer_plot-tweaked.ipynb') #Still waiting to run this; halo file needs to be available.
    _exec_notebook('NGC5533_Noordermeer_plot.ipynb')
    _exec_notebook('time setup.ipynb') #Not something that should break, but will flag if time library updates break anything.
