import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import rec_groupby

files = {
        'g45-24_300-0-0.003-1_all.log': (300, 0, 0.003, 1),
        'g45-24_300-0.2-0.003-1_all.log': (300, 0.2, 0.003, 1),
        'g45-24_300-0.6-0.003-1_all.log': (300, 0.6, 0.003, 1),
        'g45-24_300-0.6-0.001-1_all.log': (300, 0.6, 0.001, 1),
        'g45-24_300-0.6-0.01-1_all.log': (300, 0.6, 0.01, 1),
        }
genstep = 5

markers = ['o', 's', 'v', '^', 'p', 'D']
colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k']

plt.figure()
plt.gca().set_xscale('log')

for i, fname in enumerate(files.keys()):
    data = np.loadtxt(fname, usecols=[0,1], dtype=[('gen', int), ('emin', int)])
    err = lambda a: np.std(a)/np.sqrt(len(a)-1)

    r = rec_groupby(
            data, groupby=('gen',),
            stats = [('emin', np.mean, 'emin_avg'), ('emin', err, 'emin_err')]
        )

    np.sort(r, order=['gen'])
    gen = np.array([0, 5, 10, 30, 100, 300, 1000, 3000])
    emin_avg = r.emin_avg[gen/genstep]
    emin_err = r.emin_err[gen/genstep]

    kwargs = dict(
        label = r'$(%3d,\,%2.1f,\, %4.3f,\, %1d)$' % files[fname],
        marker = markers[i % len(markers)],
        color = colors[i % len(colors)]
    )
    plt.errorbar(gen, emin_avg, emin_err, **kwargs)
    plt.xlabel('generation')
    plt.ylabel(r'$\langle E_{min} \rangle$')

plt.legend(title='($n_{pop}$, $p_{cross}$, $p_{mutate}$, $n_{cp}$)', loc='lower left')
plt.show()
