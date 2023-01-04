#############################################################################
#
# In this file, we show that H_1(S) is a 8 dimensional?
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

G = PermutationGroup([[(1,2,3,4),(5,6,7,8)], [(1,4,3,2),(5,6,7,8)], [(1,6),(2,7),(3,8),(4,5)]])
alpha = G.gens()[0]; r = G.gens()[1]; s = G.gens()[2]
x = s; y = s*r; z = alpha
number_of_edges = 6
hom = [x, x^(-1), (x*z*y), (x*z*y)^(-1), z^2, z^(-2)]
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
print(homology.module.dimension())
assert homology.module.dimension() == 6
