from ase import Atoms
from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
from ase.io import write, read
import numpy as np
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

    default_parameters = dict(verbose = True, file = 'disp.xyz',dump_format = '.npz')


    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None, atoms=None, **kwargs):

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)



    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters

        #generate structure:'disp.xyz'
        xyz = 'disp.xyz'
        atoms.info['tag'] = self.prefix + atoms.get_chemical_formula()
        write(xyz, atoms, format='xyz')



        #generate config gile:'prediction_config.py'
        #default situation
        #there is output directory and trained nnp data (master_nnp.npz, postprocess/pca.npz, training_result.yaml is stored)
        #this script generate 'disp.xyz' & prediction_config.py in the present directory

        s='c.PredictionApplication.verbose = False \n'
        s+='c.PredictionConfig.data_file = \'disp.xyz\' \n'
        s += 'c.PredictionConfig.dump_format  = \'.npz\' \n'
        s += 'c.PredictionConfig.load_dir   = \'./output\' \n'
        s += 'c.PredictionConfig.order   = 1 \n'
        s += 'c.PredictionConfig.tags = [\'*\']  \n'


        with open('prediction_config.py', 'w') as config:
            config.write(s)

    def read_results(self):
        force = []
        energy=0.0
        prdata = np.load('./output/prediction_result.npz')
        keys = prdata.files

        for ik in keys:
            if (ik.find('force') > 0):
                force_data = prdata[ik]
                force.append(force_data)
            if(ik.find('energy')>0):
                #Need to cast to float for plot trajectory
                energy=float(prdata[ik][0][0])

        self.results['forces'] = np.array(
            force).reshape((-1, 3))
        self.results['energy']=energy


    def get_forces(self, atoms):

        return FileIOCalculator.get_forces(self, atoms)







