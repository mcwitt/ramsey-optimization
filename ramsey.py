from numpy import *
from itertools import combinations

def clique_count(a, m):
    """Determines the number of m-cliques in the graph represented by the
    adjacency matrix a

    >>> clique_count(\
            array([[0,1,1,1],\
                   [1,0,1,1],\
                   [1,1,0,0],\
                   [1,1,0,0]]),\
            3)
    2
    
    >>> clique_count(\
            array([[0,1,0,0,1],\
                   [1,0,1,0,0],\
                   [0,1,0,1,0],\
                   [0,0,1,0,1],\
                   [1,0,0,1,0]]),\
            3)
    0
    """

    verts = range(len(a))
    count = 0

    # iterate over all ways of choosing m vertices
    for alpha in combinations(verts, m):
        isclique = True
        for j, k in combinations(alpha, 2):
            if not a[j,k]:
                isclique = False
                break
        if isclique:
            count += 1

    return count

def energy(a, n, m):
    return clique_count(a, m) + clique_count(a==False, n)

def draw(a, out='graph'):
    'Make an EPS drawing of the graph'
    import pyx
    c = pyx.canvas.canvas()
    N = len(a)
    verts = range(N)
    xy = lambda j: (5*sin(2*pi*j/N), 5*cos(2*pi*j/N))
    for j, k in combinations(verts, 2):
        color = pyx.color.rgb.red if a[j, k] else pyx.color.rgb.blue
        xj, yj = xy(j)
        xk, yk = xy(k)
        c.stroke(pyx.path.line(xj, yj, xk, yk), [color])
    for k in verts:
        x, y = xy(k)
        c.fill(pyx.path.circle(x, y, 0.2), [pyx.color.rgb.black])

    c.writeEPSfile(out)

def read(filename):
    with open(filename, 'r') as f:
        nv = int(f.readline())
        r = int(f.readline())
        s = int(f.readline())

        a = zeros((nv, nv))
        for k in xrange(nv):
            for j in xrange(k):
                a[j, k] = int(f.readline())
    return a, nv, r, s

if __name__=='__main__':
    import sys
    if len(sys.argv) == 2:
        a, nv, r, s = read(sys.argv[1])
        print energy(a, r, s)
        draw(a, 'graph.eps')
