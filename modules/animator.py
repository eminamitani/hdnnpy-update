from os import path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from config import hp
from config import file_


def update(i, ax, pred, true):
    if i != 0:
        ax.cla()

    ax.scatter(pred[i], true)
    ax.set_title('epochs={}'.format(i+1))
    min = np.min(true)
    max = np.max(true)
    ax.set_xlim([min, max])
    ax.set_ylim([min, max])


class Animator(object):
    def __init__(self, nsample):
        nfigure = 1 + 3 * hp.natom
        self.preds = np.empty((hp.nepoch, nsample, nfigure))
        self.true = np.empty((nsample, nfigure))

    def set_pred(self, m, E_pred, F_pred):
        # nepoch x nsample x 1+3*natom -> 1+3*natom x nepoch x nsample
        self.preds[m] = np.c_[E_pred, F_pred.reshape((-1, 3 * hp.natom))]

    def set_true(self, E_true, F_true):
        # nsample x 1+3*natom -> 1+3*natom x nsample
        self.true = np.c_[E_true, F_true.reshape((-1, 3 * hp.natom))]

    def save_fig(self):
        self.preds = self.preds.transpose(2, 0, 1)
        self.true = self.true.transpose(1, 0)

        for i, (pred, true) in enumerate(zip(self.preds, self.true)):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            anime = FuncAnimation(fig, update, fargs=(ax, pred, true), interval=250, frames=hp.nepoch)
            if i == 0:
                filename = 'energy.gif'
            elif i % 3 == 0:
                filename = 'force_{}x.gif'.format((i-1)/3+1)
            elif i % 3 == 1:
                filename = 'force_{}y.gif'.format((i-1)/3+1)
            elif i % 3 == 2:
                filename = 'force_{}z.gif'.format((i-1)/3+1)
            anime.save(path.join(file_.fig_dir, filename), writer='imagemagick')
            plt.close(fig)