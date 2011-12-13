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

def indset_count(a, n):
    '''Determines the number of n-independent sets in the graph represented by
    the adjacency matrix a'''

    return clique_count(a==False, n)

def energy(a, n, m):
    return clique_count(a, m) + indset_count(a, n)

def sweep(a, n, m, beta):
    'Sweep with brute-force energy delta calculation'
    verts = range(len(a))
    e0 = energy(a, n, m)
    for j, k in combinations(verts, 2):
        a[j,k] = not a[j,k]
        e1 = energy(a, n, m)
        de = e1 - e0
        if de < 0 or random.random() < exp(-beta*de):
            e0 = e1
        else:
            a[j,k] = not a[j,k]

def draw(a, out='graph'):
    'Make an EPS drawing of the graph'
    import pyx
    c = pyx.canvas.canvas()
    N = len(a)
    for k in xrange(N):
        for j in xrange(k):
            if a[j,k]: color=pyx.color.rgb.red
            else: color=pyx.color.rgb.blue
            coords = lambda j: (5*sin(2*pi*j/N), 5*cos(2*pi*j/N))
            xj, yj = coords(j)
            xk, yk = coords(k)
            c.stroke(pyx.path.line(xj, yj, xk, yk), [color])
    for k in xrange(N):
        x, y = coords(k)
        c.fill(pyx.path.circle(x, y, 0.2), [pyx.color.rgb.black])

    c.writeEPSfile(out)

