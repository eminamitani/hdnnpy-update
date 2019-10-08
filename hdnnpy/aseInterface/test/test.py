import sys
import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/../' )
from aseInterface import *
from ase.build import bulk

hdnnpy=hdnnpy()
hdnnpy.set_label('./Crystal')
a1 = bulk('Si', 'diamond', a=3.6)
a1.set_calculator(hdnnpy)
e=a1.get_potential_energy()
print(e)
f=a1.get_forces()
print(f)


