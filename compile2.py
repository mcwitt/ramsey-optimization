from scipy.misc import comb
import subprocess
import sys

if len(sys.argv) == 4:
    r, s, nv = (int(arg) for arg in sys.argv[1:4])
else:
    r   = int(input('r   = '))
    s   = int(input('s   = '))
    nv  = int(input('N_v = '))


nsgr = comb(nv, r, exact=True)
nsgs = comb(nv, s, exact=True)
nsgfer = comb(nv-2, r-2, exact=True)
nsgfes = comb(nv-2, s-2, exact=True)

cmd = "make ramsey2 NV={} R={} S={} NSGR={} NSGS={} NSGFER={} NSGFES={}".format(
            nv, r, s, nsgr, nsgs, nsgfer, nsgfes
        ).split()

p = subprocess.Popen(cmd)
p.wait()
