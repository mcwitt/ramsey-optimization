from pylab import *

files = [
        ['pt30_T2', 'pt31_T2', 'pt32_T2', 'pt33_T2'],
        ['sa30_1-0.1-10-10-5', 'sa31_1-0.1-10-50-5', 'sa32_1-0.1-10-100-15', 'sa33_1-0.2-300-1000-25'],
        ['de30_100-100-1', 'de31_100-100-1', 'de32_100-200-2', 'de33_1000-10000-12'],
        ['eo30_1.6-100000', 'eo31_1.6-100000', 'eo32_1.6-100000', 'eo33_1.6-100000']
        ]

methods = ['PT', 'SA', 'ADA', 'EO']
sizes = [30, 31, 32, 33]
symbols = ['s', 'd', 'o', '^']
colors = ['r', 'g', 'b', 'y']
rc('lines', markersize=7)

fig = figure()
ax = fig.add_subplot(111)

for method, mfiles, color, symbol in zip(methods, files, colors, symbols):
    nflip_avg = []
    nflip_err = []

    for size, f in zip(sizes, mfiles):
        data = np.loadtxt(f)
        nflip_avg.append(np.average(data))
        nflip_err.append(np.std(data) / sqrt(len(data) - 1))

    errorbar(sizes[:len(mfiles)], nflip_avg, nflip_err, label=method,
            color=color,
            marker=symbol,
            markeredgecolor=color,
            markerfacecolor='none',
            )

xlim(sizes[0]-0.5, sizes[-1]+0.5)
yscale('log')
xticks(sizes)
legend(loc='upper left')
text(0.8, 0.2, r'$N_{\mathrm{samp}}=96$', transform=ax.transAxes)
xlabel('$N_v$')
ylabel('$N_{\mathrm{flip}}/N$ to find ground state')
title('Performance of various searches on (4, 6; $N_v$)')
show()
