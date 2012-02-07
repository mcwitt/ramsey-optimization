import numpy as np
import matplotlib.pyplot as plt
import IPython

def group(table, *groupby):
    '''Return a dictionary whose keys are values of the attributes specified
    in `groupby' and whose corresponding items are lists of rows with those
    attribute values.
    
    >>> groupby(array([[1,0,1],[0,1,0],[1,1,1]]), [0,2])
    {(1,1):[0,2], (0,0):[1]}'''

    rowd = dict()
    for i, row in enumerate(table):
        key = tuple([row[attr] for attr in groupby])
        rowd.setdefault(key, []).append(i)
    return rowd

if __name__=='__main__':
    import sys, os

    if len(sys.argv) != 2: sys.exit('Need to specify input file')

    fpath = sys.argv[1]
    fname, fext = os.path.splitext(fpath)
    if fext == '.npy': data = np.load(fpath)
    else:
        dt = [('igraph', 'u4'), ('nb', 'u1')]
        data = np.loadtxt(fpath, usecols=[0,2], dtype=dt)
        np.save(fname, data)
        
    groups = group(data, 'igraph')
    nb_min, nb_max = np.inf, 0
    freqs = []

    for igraph, rows in groups.items():
        nb = data[rows]['nb']
        nb_min = min(nb_min, min(nb))
        nb_max = max(nb_max, max(nb))

    x = np.arange(nb_min, nb_max+1, 1)
    bins = np.arange(nb_min-0.5, nb_max+1.5, 1)

    for igraph, rows in groups.items():
        nb = data[rows]['nb']
        freq, b = np.histogram(nb, bins)
        freqs.append(freq)

    freqs = np.array(freqs)

    freq_avg = np.average(freqs, axis=0)
    freq_err = np.std(freqs, axis=0)/np.sqrt(freqs.shape[0]-1.)

    plt.bar(x, freq_avg, yerr=freq_err, align='center', width=1, color='gray')
    plt.show()
