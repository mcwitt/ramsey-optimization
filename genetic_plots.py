import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import rec_groupby

files = {
        'g45-24_200-0-0.003.log': (200, 0, 0.003),
        'g45-24_200-0.2-0.003.log': (200, 0.2, 0.003),
        'g45-24_200-0.6-0.003.log': (200, 0.6, 0.003),
        }

markers = ['o', 's', 'v', '^', 'p', 'D']
colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k']

for i, fname in enumerate(files.keys()):
    data = np.loadtxt(fname, usecols=[0,1], dtype=[('gen', int), ('emin', int)])
    err = lambda a: np.std(a)/np.sqrt(len(a)-1)

    r = rec_groupby(
            data, groupby=('gen',),
            stats = [('emin', np.mean, 'emin_avg'), ('emin', err, 'emin_err')]
        )

    kwargs = dict(
        label = r'$(%d,\,%.1g,\, %.1g)$' % files[fname],
        marker = markers[i % len(markers)],
        color = colors[i % len(colors)]
    )
    plt.errorbar(r.gen, r.emin_avg, r.emin_err, **kwargs)

plt.legend(title=r'$(N,\,p_{cross},\,p_{mutate})$', loc='upper right')
plt.show()
