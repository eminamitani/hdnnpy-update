from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
import os

class hdnnpy(FileIOCalculator):

    '''
    hdnnpy predict to single structure --> getting total energy & force
    '''

    '''
    hdnnpy calculator object
    
    
    
    '''

    implemented_properties=['energy','forces']
    command='hdnnpy predict'

    default_parameters = dict(verbose = True, file = 'disp.xyz',dump_format = '.npz',)

    def __init__(self, restart=None, ignore_bad_restart_file=False, label='hdnnpy', atoms=None):

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters

        #generate structure:'disp.xyz'
        xyz = 'disp.xyz'
        atoms.info['tag'] = self.prefix + atoms.get_chemical_formula()
        ase.io.write(xyz, atoms, format='xyz')



        #generate config gile:'prediction_config.py'

        s='c.PredictionApplication.verbose = True \n'
        s+='c.PredictionConfig.data_file = \'disp.xyz\' \n'
        s += 'c.PredictionConfig.dump_format  = \'.npz\' \n'
        s += 'c.PredictionConfig.load_dir   = \'. \' \n'
        s += 'c.PredictionConfig.order   = 1 \n'
        s += 'c.PredictionConfig.tags = [\'*\']  \n'


        with open('prediction_config.py', 'w') as config:
            config.write(s)

