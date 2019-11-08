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
#    for file in [f for f in os.listdir('../fitting/Static_Notebooks/') if f.endswith('.ipynb')]: #I found this online, can't get it to work
#        _exec_notebook('../fitting/Static_Notebooks/')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-BlackHole.ipynb')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-Halo.ipynb')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-Bulge.ipynb')
    _exec_notebook('../fitting/Static_Notebooks/Fitting-Disk.ipynb')
    _exec_notebook('BlackHole_Velocity.ipynb')
    _exec_notebook('Bulge_RotationCurve_n2_7.ipynb')
    _exec_notebook('Disk_Velocity_kpc.ipynb')
    _exec_notebook('Halo_Velocity.ipynb')
    _exec_notebook('NGC5533_Total.ipynb')
    _exec_notebook('../LibraryAutoTest.ipynb')
