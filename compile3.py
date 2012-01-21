from scipy.misc import comb
import subprocess
import sys

if len(sys.argv) == 3:
    s, nv = int(sys.argv[1]), int(sys.argv[2])
else:
    s   = int(input('s   = '))
    nv  = int(input('N_v = '))


ned = nv*(nv-1)/2
neds = s*(s-1)/2
nsg = comb(nv, s, exact=True)
nsgfe = comb(nv-2, s-2, exact=True)

cmd = "make ramsey3 NV={} S={} NED={} NSG={} NSGFE={}".format(nv, s, ned, nsg, nsgfe).split()
p = subprocess.Popen(cmd)
p.wait()
