#############################################################################
#
# In this file, we show that H_1(S) is 2 copies of the standard rep and that
# the action is arithmetic.
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

G = DihedralGroup(4)
r = G("(1,2,3,4)"); s = G("(1,4)(2,3)")
number_of_edges = 6
hom = [s, s, s*r, s*r, r^2, r^2]
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
assert homology.module.dimension() == 4

##### Building the irrep #####

rep = {r: matrix([[0, -1], [1, 0]]), s: matrix([[0,1],[1,0]])}
for i in range(0, 4):
    for j in range(0,2):
        rep[r^i*s^j] = rep[r]^i * rep[s]^j 

for g in G:
    assert deck_group_actions[g].trace() == 2*rep[g].trace()


##### Building the isomorphism between homology and rep #####

A = identity_matrix(4)
B = zero_matrix(4)

for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g])
    B = B + g_bad*A*(g_good.inverse())

# Checking that we indeed have an isomorphism
assert B.determinant() != 0
for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g])
    assert g_bad*B*(g_good.inverse()) == B


##### Lifting twists #####

twist_matrices = []

for power in range (1,4):
    for curve in liftable_curves(cover, power):
        T = action_of_twist_on_homology(cover, homology, curve, power)
        # Can verify that T commutes with deck group
        for g in G:
            assert T.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T.matrix().transpose()
        
        R = B.inverse()*T.matrix().transpose()*B
        # R will be a block matrix where each block is a multiple of the identity
        T_in_isotypic = matrix([[R[0][0], R[0][2]],[R[2][0],R[2][2]]])
        
        # Can verify that T has integer entries
        for i in range(0,2):
            for j in range(0,2):
                assert T_in_isotypic[i][j] in ZZ
        assert T_in_isotypic.determinant() == 1

        twist_matrices.append(T_in_isotypic)


##### Check finite index #####

print(is_finite_index(twist_matrices))