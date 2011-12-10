from math import *
from itertools import combinations

def clique_count(N, g, m):
    '''Determines the number of m-cliques in the N-vertex graph represented by
    g, the length N(N-1)/2 binary string constructed from the adjacency matrix
    elements as follows:

    g = a_{2,1} ... a_{N,1} a_{3,2} ... a_{N,2} ... a_{N,N-1}'''

    assert(N*(N-1)/2 == len(g))

    g = list(g)
    count = 0

    # iterate over all ways of choosing m vertices
    for alpha in combinations(range(N), m):
        for pair in combinations(alpha, 2):





def indset_count(N, g, n):
    '''Determines the number of n-independent sets in g'''
    pass

def energy(g, n, m):
    pass
