#############################################################################
#
# In this file, we show that H_1(S) is a 4-dimensional irrep together with two 1-dimensional irreps
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

G = PermutationGroup([[(1,2,3,4),(5,6,7,8)], [(1,4,3,2),(5,6,7,8)], [(1,6),(2,7),(3,8),(4,5)]])
alpha = G.gens()[0]; r = G.gens()[1]; s = G.gens()[2]
x = s; y = s*r; z = alpha
number_of_edges = 6
hom = [x,x^(-1), x*z*y, (x*z*y)^(-1), y, y^(-1)]
gluing = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4}
edge_orientations = {0:1,1:-1,2:1,3:-1,4:1,5:-1}

base_surface = BaseSurface(number_of_edges, gluing, edge_orientations)
cover = Cover(base_surface, G, hom)

homology = HomologyGroup(cover)
deck_group_actions = {g:homology.action_of_deck_group_on_homology(g).matrix().transpose() for g in G}
# Sage assumes matrices act on the right, so we take the transpose to get an actual homomorphism
for g in G:
    for h in G:
        assert deck_group_actions[g]*deck_group_actions[h] == deck_group_actions[g*h]
assert homology.module.dimension() == 6


##### Building the irreps #####

big_rep = {alpha: matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]),
    r: matrix([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]]),
    s: matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])}

small_rep = {alpha: -1, r: 1, s: -1}

for i in range(0,2):
    for j in range(0,2):
        for k in range(0, 4):
            big_rep[alpha^i*s^j*r^k] = big_rep[alpha]^i*big_rep[s]^j*big_rep[r]^k
            small_rep[alpha^i*s^j*r^k] = small_rep[alpha]^i*small_rep[s]^j*small_rep[r]^k


##### Verifications #####

# First, we check that we have the right representation
for g in G.list():
    assert deck_group_actions[g].trace() == 2*small_rep[g] + big_rep[g].trace()

# Next, we can check that our big_rep is indeed a surjection onto M_2(Q(i))
mat_to_vec = lambda m : vector([m[i][j] for i in range(0,4) for j in range(0,4)])
vectors = [mat_to_vec(big_rep[g]) for g in G.list()]
assert dim(span(vectors, QQ)) == 8