from choose import choose
import subprocess
import sys

if len(sys.argv) == 4:
    r, s, nv = (int(arg) for arg in sys.argv[1:4])
else:
    r   = int(input('r   = '))
    s   = int(input('s   = '))
    nv  = int(input('N_v = '))


ned = nv*(nv-1)/2
nsgr = choose(nv, r)
nsgs = choose(nv, s)
nsgfer = choose(nv-2, r-2)
nsgfes = choose(nv-2, s-2)

cmd = "make NV={} R={} S={} NED={} NSGR={} NSGS={} NSGFER={} NSGFES={}".format(
            nv, r, s, ned, nsgr, nsgs, nsgfer, nsgfes
        ).split()

p = subprocess.Popen(cmd)
p.wait()
