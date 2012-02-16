#!/usr/bin/python

from choose import choose
import sys

outfile = 'defs.h'

if len(sys.argv) == 4:
    r, s, nv = (int(arg) for arg in sys.argv[1:4])
else:
    r   = int(input('r   = '))
    s   = int(input('s   = '))
    nv  = int(input('N_v = '))

if r > s: r, s = s, r

params = [
    ('R',       r),
    ('S',       s),
    ('NV',      nv),
    ('NED',     nv*(nv-1)/2),
    ('NSGR',    choose(nv, r)),
    ('NSGS',    choose(nv, s)),
    ('NSGFER',  choose(nv-2, r-2)),
    ('NSGFES',  choose(nv-2, s-2))
]

f = open(outfile, 'w')
f.write('/* This file is generated automatically by gendefs.py */\n')

for item in params:
    f.write('#define {:<8} {}\n'.format(*item))

f.close()
