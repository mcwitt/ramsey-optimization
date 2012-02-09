import numpy as np
import matplotlib.pyplot as plt
import IPython

def nbhistplot(fname, fmt='r-o', label=None):
    data = np.loadtxt(fname, unpack=True)
    if not label: label = fname

    if len(data) == 3:
        nb, count, err = data
        plt.errorbar(nb, count, err, fmt=fmt, label=label)
        #plt.bar(nb, count, yerr=err, align='center', width=0.25, color=color)
        return nb, count, err
    else:
        nb, count = data
        plt.plot(nb, count, fmt, label=label)
        return nb, count

if __name__=='__main__':
    import sys

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fmts = ['r-o', 'b-s', 'g-v']
    data = dict()

    for fname, fmt in zip(sys.argv[1:], fmts):
        data[fname] = nbhistplot(fname, fmt)

    plt.legend(loc='upper left')
    IPython.embed()
