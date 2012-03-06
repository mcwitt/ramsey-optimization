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
    plt.xlabel('$D_{max}$')
    plt.ylabel(r'$\langle E \rangle$')
    plt.legend(loc='upper left')
    return r

def emin_run_hist(data, fmt='r-o', label=None):
    r = rec_groupby(data, groupby=('run',),
            stats = [('emin', np.min, 'emin_run')])

    r = r[:-1]
    
    counts, bins = np.histogram(r.emin_run,
            bins=range(min(r.emin_run), max(r.emin_run)))
    plt.plot(bins[:-1], counts, fmt, label=label)
    plt.xlabel('$E_{min}$')
    plt.ylabel('count')
    plt.legend(loc='upper left')
    return bins[:-1], counts

if __name__=='__main__':
    import sys

    plt.rc('legend', fontsize=10)
    markers = ['o', 's', 'v', '^', 'p', 'D']
    colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k']
    results = [] 
    files = sys.argv[1:]

    for i in xrange(len(files)):
        fname = files[i]

        data = np.loadtxt(fname, usecols=[0,2,6],
                dtype=[ ('run', 'i4'), ('dmax','i4'), ('emin','i4') ])

        fmt = '%s-%s' % (colors[i % len(colors)], markers[i % len(markers)])
        plt.figure(1)
        emin_stage = emin_stage_plot(data, fmt, fname)
        plt.figure(2)
        emin_run = emin_run_hist(data, fmt, fname)

        results.append(dict(
            fname = fname,
            data = data,
            emin_stage = emin_stage,
            emin_run = emin_run
            )
        )

    IPython.embed()
