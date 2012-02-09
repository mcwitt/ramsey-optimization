import numpy as np
import matplotlib.pyplot as plt
import IPython
from matplotlib.mlab import rec_groupby

def energyplot(fname, fmt='r-o', label=None):
    data = np.loadtxt(fname, usecols=[3,7], dtype=[('dmax','i4'),('emins','i4')])
    if not label: label = fname

    stats = (   ('emins', np.mean, 'emins_avg'),
                ('emins', lambda a: np.std(a)/np.sqrt(len(a)-1.), 'emins_err'))

    result = rec_groupby(data, groupby=('dmax',), stats=stats)
    plt.errorbar(result.dmax, result.emins_avg, result.emins_err, fmt=fmt, label=label)
    return result

if __name__=='__main__':
    import sys

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fmts = ['r-o', 'b-s', 'g-v', 'y-^', 'm-p', 'c-x', 'k-D', 'r-s']
    data = dict()

    for fname, fmt in zip(sys.argv[1:], fmts):
        data[fname] = energyplot(fname, fmt)

    plt.legend(loc='upper left')
    IPython.embed()
