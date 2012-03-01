import numpy as np
import matplotlib.pyplot as plt
import sys

files = sys.argv[1:]

for i in xrange(len(files)):
    fname = files[i]
    data = np.loadtxt(fname, usecols=[0,1], dtype=[ ('gen', 'i4'), ('avg', 'f8') ])
    plt.plot(data['gen'], data['avg'])

plt.show()
