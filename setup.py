import sys
from cx_Freeze import setup, Executable

"""
import os.path
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6') 
"""

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {
    "packages": ["os"],
    "includes": ["scipy.spatial.ckdtree", "scipy._distributor_init"],
                     }

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "HEA Calc",
        version = "1.0",
        description = "Finding new HEA alloys",
        options = {"build_exe": build_exe_options},
        executables = [Executable("hea_modeliranje_AI.py", base=base)])
