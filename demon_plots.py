import numpy as np
import matplotlib.pyplot as plt
import IPython
from matplotlib.mlab import rec_groupby

def emin_stage_plot(data, fmt='r-o', label=None):

    err = lambda a: np.std(a)/np.sqrt(len(a)-1)
    r = rec_groupby(
        data, groupby=('dmax',),
        stats = [('emin', np.mean, 'emin_avg'), ('emin', err, 'emin_err')]
    )
    plt.errorbar(r.dmax, r.emin_avg, r.emin_err, fmt=fmt, label=label)
    plt.legend(loc='upper left')
    return r

def emin_run_hist(data, fmt='r-o', label=None):
    r = rec_groupby(data, groupby=('run',),
            stats = [('emin', np.min, 'emin_run')])

    r = r[:-1]
    
    counts, bins = np.histogram(r.emin_run,
            bins=range(min(r.emin_run), max(r.emin_run)))
    plt.plot(bins[:-1], counts, fmt, label=label)
    return bins[:-1], counts

if __name__=='__main__':
    import sys

    fmts = ['r-o', 'b-s', 'g-v', 'y-^', 'm-p', 'c-x', 'k-D', 'r-s']
    results = {}

    for fname, fmt in zip(sys.argv[1:], fmts):
        results[fname] = {}

        data = np.loadtxt(fname, usecols=[0,3,7],
                dtype=[ ('run', 'i4'), ('dmax','i4'), ('emin','i4') ])

        plt.figure(1)
        emin_stage = emin_stage_plot(data, fmt, fname)
        plt.figure(2)
        emin_run = emin_run_hist(data, fmt, fname)

        results[fname]['data'] = data
        results[fname]['emin_stage'] = emin_stage
        results[fname]['emin_run'] = emin_run

    IPython.embed()
