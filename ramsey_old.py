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
    'Count m-cliques involving vertices j, k'
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
    'Only compute terms in the energy involving vertices j, k'
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

def sweep3(a, n, m, beta):
    verts = range(len(a))
    for j, k in combinations(verts, 2):
        rverts = list(verts)
        rverts.remove(j)
        rverts.remove(k)
        for alpha in combinations(rverts, m)


def temper(copies, n, m, betas, nswaps):
    'Attempt parallel tempering swaps'
    
    for iT in xrange(1, len(betas)):
        a, b = copies[iT-1], copies[iT]
        ba, bb = betas[iT-1], betas[iT]
        Ea, Eb = energy(a,n,m), energy(b,n,m)
        logar = (bb-ba)*(Ea-Eb)

        if logar < 0 or random.random() < exp(logar):
            copies[iT-1], copies[iT] = b, a
            nswaps[iT] += 1

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

def test():
    temps = array(r_[0.01, 1:100.:11j])
    betas = 1./temps
    nswaps = zeros(12)
    copies = []
    for i in xrange(12):
        copies.append(random.random((17,17)) < 0.5)

    for i in xrange(1000):
        for copy, beta in zip(copies, betas):
            sweep2(copy, 4, 4, beta)
            print '%f %3d' % (1./beta, energy(copy, 4, 4))
        temper(copies, 4, 4, betas, nswaps)
        if energy(copies[0], 4, 4) == 0:
            break
        print nswaps/(i+1.)
        print '\n'
