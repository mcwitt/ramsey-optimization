import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import rec_groupby

files = [
        (1,   '4-6-35_1_all.log'),
        (1.1, '4-6-35_1.1_all.log'),
        (1.2, '4-6-35_1.2_all.log'),
        (1.3, '4-6-35_1.3_all.log'),
        (1.4, '4-6-35_1.4_all.log'),
        (1.5, '4-6-35_1.5_all.log'),
        (1.6, '4-6-35_1.6_all.log'),
        (1.7, '4-6-35_1.7_all.log'),
        (1.8, '4-6-35_1.8_all.log'),
        (1.9, '4-6-35_1.9_all.log'),
        (2.0, '4-6-35_1.9_all.log'),
        (2.1, '4-6-35_1.9_all.log')
        ]

markers = ['o', 's', 'v', '^', 'p', 'D']
colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k']

plt.figure(1)
plt.gca().set_xscale('log')

taus = []
emin_end_avg = []
emin_end_err = []

for i, f in enumerate(files):
    tau, fname = f
    data = np.loadtxt(fname, usecols=[0,1], dtype=[('nsweep', int), ('emin', int)])
    err = lambda a: np.std(a)/np.sqrt(len(a)-1)

    r = rec_groupby(
            data, groupby=('nsweep',),
            stats = [('emin', np.mean, 'emin_avg'), ('emin', err, 'emin_err')]
        )

    np.sort(r, order=['nsweep'])
    nsweep = np.array([1, 3, 10, 30, 100, 300, 999])
    emin_avg = r.emin_avg[nsweep]
    emin_err = r.emin_err[nsweep]

    kwargs = dict(
        label = r'$%.2f$' % tau,
        marker = markers[i % len(markers)],
        color = colors[i % len(colors)]
    )

    plt.errorbar(nsweep, emin_avg, emin_err, **kwargs)

    taus.append(tau)
    emin_end_avg.append(emin_avg[-1])
    emin_end_err.append(emin_err[-1])

plt.figure(1)
plt.text(0.5, 0.8, '(4, 6; 35)',
        horizontalalignment='center',
        transform=plt.gca().transAxes)
plt.xlabel('$N_{sweep}$')
plt.ylabel(r'$\langle E_{min} \rangle$')
plt.legend(title=r'$\tau$', loc='lower left')

plt.figure(2)
plt.errorbar(taus, emin_end_avg, emin_end_err, marker='o', label='(4, 6; 35)')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\langle E_{min} \rangle$ after 1000 sweeps')
plt.legend(title=r'$(r,\,s;\,N_v)$', loc='lower left')

plt.show()
