from choose import choose
import subprocess
import sys

if len(sys.argv) == 3:
    s, nv = int(sys.argv[1]), int(sys.argv[2])
else:
    s   = int(input('s   = '))
    nv  = int(input('N_v = '))


ned = nv*(nv-1)/2
nsg = choose(nv, s)
nsgfe = choose(nv-2, s-2)

cmd = "make ramsey3 NV={} S={} NED={} NSG={} NSGFE={}".format(nv, s, ned, nsg, nsgfe).split()
p = subprocess.Popen(cmd)
p.wait()
