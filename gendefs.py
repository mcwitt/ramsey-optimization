#!/usr/bin/python

from choose import choose
import sys

outfile = 'defs.h'

if len(sys.argv) == 4:
    r, s, nv = (int(arg) for arg in sys.argv[1:4])
else:
    r   = int(input('R  = '))
    s   = int(input('S  = '))
    nv  = int(input('NV = '))

if r > s: r, s = s, r

params = [
    ('R',       r,                  'red clique size'),
    ('S',       s,                  'blue clique size'),
    ('NV',      nv,                 'number of vertices'),
    ('NED',     nv*(nv-1)/2,        'number of edges = NV(NV-1)/2'),
    ('NSGR',    choose(nv, r),      'number of R-subgraphs = binomial(NV, R)'),
    ('NSGS',    choose(nv, s),      'number of S-subgraphs = binomial(NV, S)'),
    ('NSGFER',  choose(nv-2, r-2),  'number of R-subgraphs that use a given edge = binomial(NV-2, R-2)'),
    ('NSGFES',  choose(nv-2, s-2),  'number of S-subgraphs that use a given edge = binomial(NV-2, S-2)')
]

f = open(outfile, 'w')
f.write('/* This file is generated automatically by gendefs.py */\n')

for item in params:
    f.write('#define %-8s %-10d /* %s */\n' % item)
    #f.write('#define {:<8} {:<10} /* {} */\n'.format(*item))

f.close()
