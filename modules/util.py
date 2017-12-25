# -*- coding: utf-8 -*-

from config import mpi

from os import makedirs
from itertools import product
from collections import defaultdict
import random


def mpiprint(str):
    if mpi.rank == 0:
        print str


def mpisave(obj, *args):
    if mpi.rank == 0:
        obj.save(*args)


def mpimkdir(path):
    if mpi.rank == 0:
        makedirs(path)


def write(f, str):
    with open(f, 'a') as f:
        f.write(str)


class DictAsAttributes(dict):
    def __init__(self, dic):
        super(DictAsAttributes, self).__init__(self, **dic)

    def __dir__(self):
        return dir(super(DictAsAttributes, self))

    def __setstate__(self, state):
        self.__dict__.update(state)

    def __getattr__(self, name):
        # if name in dir(self):
        #     return self.name
        # elif name == '__setstate__':
        #     return self.__setstate__
        if name in dir(DictAsAttributes):
            return self.name

        value = self[name]
        if isinstance(value, defaultdict):
            pass
        elif isinstance(value, list):
            value = [DictAsAttributes(v) if isinstance(v, dict) else v for v in value]
        elif isinstance(value, dict):
            value = DictAsAttributes(value)
        return value

    def __setattr__(self, name, value):
        self[name] = value

    def __add__(self, other):
        new = self.copy()
        new.update(other)
        return DictAsAttributes(new)


class HyperParameter(object):
    def __init__(self, dic, random):
        self.hyperparameters = dic
        self.random = random
        self.indices = list(product(*[range(len(v)) for v in dic.values()]))

    def __iter__(self):
        if self.random:
            for i in range(self.random):
                yield DictAsAttributes({k: random.uniform(min(v), max(v)) if k != 'layer'
                                        else random.choice(v) for k, v in self.hyperparameters.items()})
        else:
            for i in range(len(self)):
                yield self.__getitem__(i)

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, i):
        return DictAsAttributes({k: v[i] for (k, v), i in zip(self.hyperparameters.items(), self.indices[i])})