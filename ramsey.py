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

def restricted_clique_count(a, m, j, k):
    N = len(a)
    rverts = range(N)
    rverts.remove(j)
    rverts.remove(k)
    count = 0

    for alpha in combinations(rverts, m-2):
        alpha = list(alpha) + [j, k]
        alpha.sort()
        isclique = True
        for u, v in combinations(alpha, 2):
            if not a[u,v]:
                isclique = False
                break
        if isclique:
            count += 1

    return count

def energy(a, n, m):
    return clique_count(a, m) + clique_count(a==False, n)

def restricted_energy(a, n, m, j, k):
    return restricted_clique_count(a, m, j, k)\
            + restricted_clique_count(a==False, n, j, k)

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

def sweep2(a, n, m, beta):
    'Sweep with more efficient energy delta calculation'
    verts = range(len(a))
    for j, k in combinations(verts, 2):
        e0 = restricted_energy(a, n, m, j, k)
        a[j,k] = not a[j,k]
        e1 = restricted_energy(a, n, m, j, k)
        de = e1 - e0
        if de > 0 and random.random() > exp(-beta*de):
            a[j,k] = not a[j,k]

def draw(a, out='graph'):
    'Make an EPS drawing of the graph'
    import pyx
    c = pyx.canvas.canvas()
    N = len(a)
    verts = range(N)
    xy = lambda j: (5*sin(2*pi*j/N), 5*cos(2*pi*j/N))
    for j, k in combinations(verts, 2):
        color = pyx.color.rgb.red if a[j,k] else pyx.color.rgb.blue
        xj, yj = xy(j)
        xk, yk = xy(k)
        c.stroke(pyx.path.line(xj, yj, xk, yk), [color])
    for k in verts:
        x, y = xy(k)
        c.fill(pyx.path.circle(x, y, 0.2), [pyx.color.rgb.black])

    c.writeEPSfile(out)

