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

def write_tex(a, out='graph.tex'):
    f = file(out, 'w')
    f.write(
        '\\documentclass{article}\n'
        '\\usepackage{tikz}\n'
        '\\begin{document}\n'
        '\\begin{tikzpicture}[scale=5]\n'
    )

    N = len(a)
    for j in range(N):
        x = cos(2*pi*j/N)
        y = sin(2*pi*j/N)
        f.write(r'\node (%d) at (%.2f,%.2f) {%d};' % (j, x, y, j+1))
        f.write('\n')

    for k in range(N):
        for j in range(k):
            if a[j,k]: color='red'
            else: color='blue'
            f.write(r'\draw[color=%s] (%d) -- (%d);' % (color, j, k))
            f.write('\n')

    f.write(
      '\\end{tikzpicture}\n'
      '\\end{document}\n'
    )
    f.close()

def draw(a, out='graph.ps'):
    import pygraphviz as pgv
    G = pgv.AGraph()
    N = len(a)
    for j in range(N):
        G.add_node(j)
    for k in range(N):
        for j in range(k):
            if a[j,k]:
                color='red'
            else:
                color='blue'

            G.add_edge(j, k, color=color)

    G.node_attr['size']=0
    G.node_attr['shape']='circle'
    G.draw(out, prog='circo')
