import sys
import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/../' )
from aseInterface import *
from ase.build import bulk

hdnnpy=hdnnpy()
a1 = bulk('Si', 'diamond', a=3.6)
hdnnpy.write_input(atoms=a1)


