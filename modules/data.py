# -*- coding: utf-8 -*-

# define variables
from config import file_
from config import mpi

# import python modules
from os import path
from re import match
from collections import defaultdict
from itertools import product
import dill
import numpy as np
from mpi4py import MPI
from chainer.datasets import TupleDataset
from chainer.datasets import split_dataset_random
from chainer.datasets import get_cross_validation_datasets_random
from quippy import Atoms
from quippy import AtomsReader
from quippy import AtomsList
from quippy import AtomsWriter
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import VaspToCm

from .util import pprint, mkdir
from .util import DictAsAttributes


def get_simple_function(name, nsample=1000):
    def make_complex(nsample):
        mesh = int(round(nsample**(1./3)))
        lin = np.linspace(0.1, 1.0, mesh, dtype=np.float32)
        x, y, z = np.meshgrid(lin, lin, lin)
        x, y, z = x.reshape(-1), y.reshape(-1), z.reshape(-1)
        input = np.c_[x, y, z]
        label = (x**2 + np.sin(y) + 3.*np.exp(z) - np.log(x*y)/2 - y/z).reshape(-1, 1)
        dinput = np.identity(3, dtype=np.float32)[None, :, :].repeat(mesh**3, axis=0)
        dlabel = np.c_[2**x - 1/(2*x), np.cos(y) - 1/(2*y) - 1/z, 3.*np.exp(z) + y/z**2].reshape(-1, 3, 1)
        return TupleDataset(input, dinput, label, dlabel)

    def make_LJ(nsample):
        input = np.linspace(0.1, 1.0, nsample, dtype=np.float32).reshape(-1, 1)
        label = 0.001/input**4 - 0.009/input**3
        dinput = np.ones((nsample, 1, 1), dtype=np.float32)
        dlabel = (0.027/input**4 - 0.004/input**5).reshape(-1, 1, 1)
        return TupleDataset(input, dinput, label, dlabel)

    def make_sin(nsample):
        input = np.linspace(-2*3.14, 2*3.14, nsample, dtype=np.float32).reshape(-1, 1)
        label = np.sin(input)
        dinput = np.ones((nsample, 1, 1), dtype=np.float32)
        dlabel = np.cos(input).reshape(-1, 1, 1)
        return TupleDataset(input, dinput, label, dlabel)

    if name == 'complex':
        dataset = make_complex(nsample)
    elif name == 'LJ':
        dataset = make_LJ(nsample)
    elif name == 'sin':
        dataset = make_sin(nsample)
    else:
        raise ValueError("function '{}' is not implemented.")
    dataset.config = name
    return dataset


class AtomicStructureDataset(TupleDataset):
    def __init__(self, hp):
        self._hp = hp

    def __len__(self):
        return len(self._datasets[0])

    @property
    def phonopy(self):
        return self._phonopy

    @property
    def composition(self):
        return self._composition

    @property
    def config(self):
        return self._config

    @property
    def input(self):
        return self._datasets[0]

    @property
    def dinput(self):
        return self._datasets[1]

    def load_xyz(self, xyz_file):
        self._data_dir = path.dirname(xyz_file)
        self._config = path.basename(self._data_dir)
        with open(path.join(self._data_dir, 'composition.dill')) as f:
            self._composition = dill.load(f)
        self._nsample = len(AtomsReader(xyz_file))
        self._natom = len(self._composition.element)
        self._atoms_objs = AtomsList(xyz_file, start=mpi.rank, step=mpi.size) if mpi.rank < self._nsample else []

        self._configure_mpi()
        Es, Fs = self._make_label(save=True)
        Gs, dGs = self._make_input(save=True)
        self._datasets = (Gs, dGs, Es, Fs)

    def load_poscar(self, poscar, dimension=[[2, 0, 0], [0, 2, 0], [0, 0, 2]], distance=0.03, save=True, scale=1.0):
        self._data_dir = path.dirname(poscar)
        self._config = path.basename(self._data_dir)
        unitcell, = AtomsList(poscar, format='POSCAR')
        if self._hp.mode == 'optimize':
            supercells = []
            for k in np.linspace(0.9, 1.1, 201):
                supercell = unitcell.copy()
                supercell.set_lattice(unitcell.lattice * k, scale_positions=True)
                supercells.append(supercell)
            self._atoms_objs = supercells
        elif self._hp.mode == 'phonon':
            unitcell.set_lattice(unitcell.lattice*scale, scale_positions=True)  # scaling
            unitcell = PhonopyAtoms(symbols=unitcell.get_chemical_symbols(),
                                    positions=unitcell.positions,
                                    numbers=unitcell.numbers,
                                    masses=unitcell.get_masses(),
                                    scaled_positions=unitcell.get_scaled_positions(),
                                    cell=unitcell.cell)
            phonon = Phonopy(unitcell,
                             dimension,
                             # primitive_matrix=primitive_matrix,
                             factor=VaspToCm)
            phonon.generate_displacements(distance=distance)
            supercells = phonon.get_supercells_with_displacements()
            self._phonopy = phonon
            self._atoms_objs = []
            for pa in supercells:
                atoms = Atoms(cell=pa.cell,
                              positions=pa.get_positions(),
                              numbers=pa.numbers,
                              masses=pa.masses)
                atoms.set_chemical_symbols(pa.get_chemical_symbols())
                self._atoms_objs.append(atoms)

        symbols = self._atoms_objs[0].get_chemical_symbols()
        composition = {'index': {k: set([i for i, s in enumerate(symbols) if s == k]) for k in set(symbols)},
                       'element': symbols}
        self._composition = DictAsAttributes(composition)
        self._nsample = len(supercells)
        self._natom = supercells[0].get_number_of_atoms()

        self._configure_mpi()
        Gs, dGs = self._make_input(save=save)
        self._datasets = (Gs, dGs)

    def reset_inputs(self, input, dinput):
        self._datasets = (input, dinput) + self._datasets[2:]

    def _configure_mpi(self):
        quo = self._nsample / mpi.size
        rem = self._nsample % mpi.size
        self._count = np.array([quo+1 if i < rem else quo
                                for i in xrange(mpi.size)], dtype=np.int32)
        self._disps = np.array([np.sum(self._count[:i])
                                for i in xrange(mpi.size)], dtype=np.int32)
        self._n = self._count[mpi.rank]  # the number of allocated samples in this node

    def _make_label(self, save):
        EF_file = path.join(self._data_dir, 'Energy_Force.npz')
        if path.exists(EF_file):
            ndarray = np.load(EF_file)
            Es = ndarray['E']
            Fs = ndarray['F']
        else:
            pprint('making {} ... '.format(EF_file), end='', flush=True)
            Es = np.empty((self._nsample, 1))
            Fs = np.empty((self._nsample, self._natom, 3))
            Es_send = np.array([data.cohesive_energy for data in self._atoms_objs]).reshape(-1, 1)
            Fs_send = np.array([data.force.T for data in self._atoms_objs]).reshape(-1, self._natom, 3)
            mpi.comm.Allgatherv((Es_send, (self._n), MPI.DOUBLE),
                                (Es, (self._count, self._disps), MPI.DOUBLE))
            num = self._natom * 3
            mpi.comm.Allgatherv((Fs_send, (self._n*num), MPI.DOUBLE),
                                (Fs, (self._count*num, self._disps*num), MPI.DOUBLE))
            if save:
                np.savez(EF_file, E=Es, F=Fs)
            pprint('done')
        Es = Es.astype(np.float32)
        Fs = Fs.astype(np.float32)
        return Es, Fs

    def _make_input(self, save):
        Gs = {}
        dGs = {}
        existing_keys = set()
        SF_file = path.join(self._data_dir, 'Symmetry_Function.npz')
        if save and path.exists(SF_file):
            ndarray = np.load(SF_file)
            if 'nsample' in ndarray and ndarray['nsample'] == self._nsample:
                for key, value in ndarray.iteritems():
                    params, GordG = path.split(key)
                    existing_keys.add(params)
                    if GordG == 'G':
                        Gs[key] = value
                    elif GordG == 'dG':
                        dGs[key] = value
            else:
                pprint('# of samples (or atoms) between Symmetry_Function.npz and structure.xyz doesn\'t match.\n'
                       're-calculate from structure.xyz.')
        new_keys = self._check_uncalculated_keys(existing_keys)

        for key, G_send, dG_send in self._calculate_symmetry_function(new_keys):
            G = np.empty((self._nsample,) + G_send.shape[1:])
            dG = np.empty((self._nsample,) + dG_send.shape[1:])
            num = G_send.size / self._n
            mpi.comm.Allgatherv((G_send, (G_send.size), MPI.DOUBLE),
                                (G, (self._count*num, self._disps*num), MPI.DOUBLE))
            num = dG_send.size / self._n
            mpi.comm.Allgatherv((dG_send, (dG_send.size), MPI.DOUBLE),
                                (dG, (self._count*num, self._disps*num), MPI.DOUBLE))
            Gs[path.join(key, 'G')] = G
            dGs[path.join(key, 'dG')] = dG

        if save and mpi.rank == 0:
            np.savez(SF_file, nsample=self._nsample, **{k: v for dic in [Gs, dGs] for k, v in dic.iteritems()})

        return np.concatenate([v for k, v in sorted(Gs.iteritems())], axis=2).astype(np.float32), \
            np.concatenate([v for k, v in sorted(dGs.iteritems())], axis=2).astype(np.float32)

    def _check_uncalculated_keys(self, existing_keys):
        all_keys = []
        for Rc in self._hp.Rc:
            key = path.join(*map(str, ['type1', Rc]))
            all_keys.append(key)
        for Rc, eta, Rs in product(self._hp.Rc, self._hp.eta, self._hp.Rs):
            key = path.join(*map(str, ['type2', Rc, eta, Rs]))
            all_keys.append(key)
        for Rc, eta, lambda_, zeta in product(self._hp.Rc, self._hp.eta, self._hp.lambda_, self._hp.zeta):
            key = path.join(*map(str, ['type4', Rc, eta, lambda_, zeta]))
            all_keys.append(key)
        new_keys = sorted(set(all_keys) - existing_keys)
        if new_keys:
            pprint('uncalculated symmetry function parameters are as follows:')
        for key in new_keys:
            pprint(key)
        return new_keys

    def _calculate_symmetry_function(self, keys):
        Gs = defaultdict(list)
        dGs = defaultdict(list)
        for at in self._atoms_objs:
            for key in keys:
                params = key.split('/')
                G, dG = getattr(self, '_'+params[0])(at, *map(float, params[1:]))
                Gs[key].append(G)
                dGs[key].append(dG)
        for key in keys:
            yield key, np.stack(Gs[key]), np.stack(dGs[key])

    def _type1(self, at, Rc):
        G = np.zeros((self._natom, 2))
        dG = np.zeros((self._natom, 2, self._natom, 3))
        for i, neighbour, homo_all, hetero_all, R, tanh, dR, _, _ in self._neighbour(at, Rc):
            g = tanh**3
            dg = -3./Rc * ((1.-tanh**2)*tanh**2)[:, None] * dR
            # G
            G[i, 0] = g[homo_all].sum()
            G[i, 1] = g[hetero_all].sum()
            # dG
            for j, (homo, indices) in enumerate(neighbour):
                if homo:
                    dG[i, 0, j] = dg.take(indices, 0).sum(0)
                else:
                    dG[i, 1, j] = dg.take(indices, 0).sum(0)
        return G, dG

    def _type2(self, at, Rc, eta, Rs):
        G = np.zeros((self._natom, 2))
        dG = np.zeros((self._natom, 2, self._natom, 3))
        for i, neighbour, homo_all, hetero_all, R, tanh, dR, _, _ in self._neighbour(at, Rc):
            g = np.exp(- eta * (R - Rs)**2) * tanh**3
            dg = (np.exp(- eta * (R - Rs)**2) * tanh**2 * (-2.*eta*(R-Rs)*tanh + 3./Rc*(tanh**2 - 1.0)))[:, None] * dR
            # G
            G[i, 0] = g[homo_all].sum()
            G[i, 1] = g[hetero_all].sum()
            # dG
            for j, (homo, indices) in enumerate(neighbour):
                if homo:
                    dG[i, 0, j] = dg.take(indices, 0).sum(0)
                else:
                    dG[i, 1, j] = dg.take(indices, 0).sum(0)
        return G, dG

    def _type4(self, at, Rc, eta, lambda_, zeta):
        G = np.zeros((self._natom, 3))
        dG = np.zeros((self._natom, 3, self._natom, 3))
        for i, neighbour, homo_all, hetero_all, R, tanh, dR, cos, dcos in self._neighbour(at, Rc):
            ang = 1. + lambda_ * cos
            rad1 = np.exp(-eta * R**2) * tanh**3
            rad2 = np.exp(-eta * R**2) * tanh**2 * (-2.*eta*R*tanh + 3./Rc*(tanh**2 - 1.0))
            ang[np.eye(len(R), dtype=bool)] = 0
            g = 2.**(1-zeta) * ang**zeta * rad1[:, None] * rad1[None, :]
            dg_radial_part = 2.**(1-zeta) * ang[:, :, None]**zeta * rad2[:, None, None] * rad1[None, :, None] * dR[:, None, :]
            dg_angular_part = zeta * lambda_ * 2.**(1-zeta) * ang[:, :, None]**(zeta-1) * rad1[:, None, None] * rad1[None, :, None] * dcos
            dg = dg_radial_part + dg_angular_part

            # G
            G[i, 0] = g.take(homo_all, 0).take(homo_all, 1).sum() / 2.0
            G[i, 1] = g.take(hetero_all, 0).take(hetero_all, 1).sum() / 2.0
            G[i, 2] = g.take(homo_all, 0).take(hetero_all, 1).sum()
            # dG
            for j, (homo, indices) in enumerate(neighbour):
                if homo:
                    dG[i, 0, j] = dg.take(indices, 0).take(homo_all, 1).sum((0, 1))
                    dG[i, 2, j] = dg.take(indices, 0).take(hetero_all, 1).sum((0, 1))
                else:
                    dG[i, 1, j] = dg.take(indices, 0).take(hetero_all, 1).sum((0, 1))
                    dG[i, 2, j] = dg.take(indices, 0).take(homo_all, 1).sum((0, 1))
        return G, dG

    def memorize_generator(f):
        cache = defaultdict(list)
        done = []

        def helper(self, at, Rc):
            if id(at) not in done:
                cache.clear()
                done.append(id(at))

            if Rc not in cache:
                for ret in f(self, at, Rc):
                    cache[Rc].append(ret)
                    yield cache[Rc][-1]
            else:
                for ret in cache[Rc]:
                    yield ret
        return helper

    @memorize_generator
    def _neighbour(self, at, Rc):
        at.set_cutoff(Rc)
        at.calc_connect()
        for i in xrange(self._natom):
            n_neighb = at.n_neighbours(i+1)
            if n_neighb == 0:
                continue

            r = np.zeros((n_neighb, 3))
            element = self._composition.element[i]
            neighbour = [(j_prime in self._composition.index[element], []) for j_prime in xrange(self._natom)]
            homo_all, hetero_all = [], []
            for j, neighb in enumerate(at.connect[i+1]):
                r[j] = neighb.diff
                neighbour[neighb.j-1][1].append(j)
                if neighbour[neighb.j-1][0]:
                    homo_all.append(j)
                else:
                    hetero_all.append(j)
            R = np.linalg.norm(r, axis=1)
            tanh = np.tanh(1-R/Rc)
            dR = r / R[:, None]
            cos = dR.dot(dR.T)
            # cosine(j - i - k) differentiate w.r.t. "j"
            # dcos = - rj * cos / Rj**2 + rk / Rj / Rk
            dcos = - r[:, None, :]/R[:, None, None]**2 * cos[:, :, None] \
                + r[None, :, :]/(R[:, None, None] * R[None, :, None])
            yield i, neighbour, homo_all, hetero_all, R, tanh, dR, cos, dcos


class DataGenerator(object):
    def __init__(self, hp, precond):
        self._hp = hp
        self._precond = precond
        self._data_dir = path.dirname(file_.xyz_file)
        self._config_type_file = path.join(self._data_dir, 'config_type.dill')
        if not path.exists(self._config_type_file):
            self._parse_xyzfile()

        with open(self._config_type_file, 'r') as f:
            config_type = dill.load(f)
        self._datasets = []
        self._elements = set()
        self._length = 0
        for type in file_.config:
            for config in filter(lambda config: match(type, config) or type == 'all', config_type):
                xyz_file = path.join(self._data_dir, config, 'structure.xyz')
                dataset = AtomicStructureDataset(self._hp)
                dataset.load_xyz(xyz_file)
                self._precond.decompose(dataset)
                self._datasets.append(dataset)
                self._elements.update(dataset.composition.element)
                self._length += len(dataset)

    def __iter__(self):
        """
        first iteration: training - validation set
        second iteration: config_type
        """
        if self._hp.mode == 'cv':
            splited = [[(train, val, d.config, d.composition)
                        for train, val in get_cross_validation_datasets_random(d, n_fold=self._hp.kfold)]
                       for d in self._datasets]
            for dataset in zip(*splited):
                yield dataset, self._elements
        elif self._hp.mode == 'training':
            splited = [split_dataset_random(d, len(d)*9/10) + (d.config, d.composition) for d in self._datasets]
            yield splited, self._elements

    def __len__(self):
        return self._length

    def _parse_xyzfile(self):
        if mpi.rank == 0:
            pprint('config_type.dill is not found.\nLoad all data from xyz file ... ', end='', flush=True)
            config_type = set()
            alldataset = defaultdict(list)
            for data in AtomsReader(file_.xyz_file):
                config = data.config_type
                config_type.add(config)
                alldataset[config].append(data)
            with open(self._config_type_file, 'w') as f:
                dill.dump(config_type, f)

            for config in config_type:
                composition = {'index': defaultdict(set), 'element': []}
                for i, atom in enumerate(alldataset[config][0]):
                    composition['index'][atom.symbol].add(i)
                    composition['element'].append(atom.symbol)
                composition = DictAsAttributes(composition)

                data_dir = path.join(self._data_dir, config)
                mkdir(data_dir)
                writer = AtomsWriter(path.join(data_dir, 'structure.xyz'))
                for data in alldataset[config]:
                    if data.cohesive_energy < 0.0:
                        writer.write(data)
                writer.close()
                with open(path.join(data_dir, 'composition.dill'), 'w') as f:
                    dill.dump(composition, f)
            pprint('done')

        mpi.comm.Barrier()
