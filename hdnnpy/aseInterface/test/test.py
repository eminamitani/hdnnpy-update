import sys
import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/../' )
from aseInterface import *
from ase.build import bulk, make_supercell
from ase.optimize import BFGS

hdnnpy=hdnnpy()
hdnnpy.set_label('./Crystal')
a1 = bulk('Si', 'diamond', a=3.6)
super=make_supercell(a1,[[2,1,1],[1,2,1],[1,1,2]])
super.set_calculator(hdnnpy)
e=super.get_potential_energy()
print(e)
f=super.get_forces()
print(f)
dyn=BFGS(super)
dyn.run(fmax=0.05)
