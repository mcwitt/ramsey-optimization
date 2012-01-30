from scipy import *
from matplotlib.pyplot import *
from itertools import combinations

def clique_count(a, s):
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
    for c in combinations(verts, s):
        isclique = True
        for j, k in combinations(c, 2):
            if not a[j, k]:
                isclique = False
                break
        if isclique:
            count += 1

    return count

def energy(a, r, s):
    return clique_count(a==False, r) + clique_count(a, s)

def random_graph(nv):
    return (rand(nv, nv) < 0.5).astype(int)

def draw(a, r, s):
    N = len(a)
    verts = range(N)
    xy = lambda j: (sin(2*pi*j/N), cos(2*pi*j/N))

    def add_graph(fig, row=1, col=1, index=1):
        ax = fig.add_subplot(row, col, index)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_aspect(1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        return ax

    def draw_verts(ax):
        for j in verts:
            c = Circle(xy(j), radius=0.05, color='k', zorder=10)
            ax.add_patch(c)

    def draw_edge(ax, j, k):
        color = 'b' if a[j, k] == 1 else 'r'
        xj, yj = xy(j)
        xk, yk = xy(k)
        l = Line2D([xj, xk], [yj, yk], color=color)
        ax.add_line(l)

    # draw all edges
    fig = figure()
    ax = add_graph(fig)
    draw_verts(ax)
    for j, k in combinations(verts, 2):
        draw_edge(ax, j, k)

    # draw cliques
    fig = figure()
    ncols = 3.
    nrows = int(ceil(energy(a, r, s) / ncols))
    i = 1
    for g, size in [(a, s), (a == False, r)]:
        for c in combinations(verts, size):
            isclique = True
            for j, k in combinations(c, 2):
                if not g[j, k]:
                    isclique = False
                    break
            if isclique:
                ax = add_graph(fig, nrows, ncols, i)
                draw_verts(ax)
                for j, k in combinations(c, 2):
                    draw_edge(ax, j, k)
                i += 1

    show()

def load_graph(filename):
    with open(filename, 'r') as f:
        nv = int(f.readline())
        a = empty((nv, nv))
        for k in xrange(nv):
            for j in xrange(k):
                a[j, k] = int(f.readline())
    return a

if __name__=='__main__':
    import sys
    if len(sys.argv) == 2:
        a = load_graph(sys.argv[1])
        r = input('r = ')
        s = input('s = ')
        print 'E = {}'.format(energy(a, r, s))
        #draw(a, r, s)
