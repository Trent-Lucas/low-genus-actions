#############################################################################
#
# In this file, we show that H_1(S) is 1 copy of M_2(Q) and 2 copies of a 1-dimensional irrep,
# and we show the action is arithmetic
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

G = DihedralGroup(6)
r = G("(1,2,3,4,5,6)"); s = G("(1,6)(2,5)(3,4)")
number_of_edges = 6
hom = [s,s,s*r^2,s*r^2,r^3,r^3]
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


##### Building the irrep #####

rep = {r: matrix([[0,-1],[1,1]]), s: matrix([[0,1],[1,0]])}
small_rep = {r: -1, s: -1}
for i in range(0, 6):
    for j in range(0,2):
        rep[r^i*s^j] = rep[r]^i * rep[s]^j 
        small_rep[r^i*s^j] = small_rep[r]^i * small_rep[s]^j

for g in G:
    assert deck_group_actions[g].trace() == 2*rep[g].trace() + 2*small_rep[g]


##### Building the isomorphism between homology and rep #####

A = matrix([[2,0,0,-1,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
assert A.determinant() != 0
B = zero_matrix(6)

for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g], small_rep[g]*identity_matrix(2))
    B = B + g_bad*A*(g_good.inverse())

# Checking that we indeed have an isomorphism
assert B.determinant() != 0
for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g], small_rep[g]*identity_matrix(2))
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
        # The top-left 4x4 minor of R will be a block matrix where each block is a multiple of the identity
        T_in_isotypic = matrix([[R[0][0], R[0][2]],[R[2][0],R[2][2]]])
        assert T_in_isotypic.determinant() == 1
        
        # Can verify that T has integer entries
        int_entries = True
        for i in range(0,2):
            for j in range(0,2):
                if T_in_isotypic[i][j] not in ZZ:
                    int_entries = False
        if int_entries:
            twist_matrices.append(T_in_isotypic)


##### Check finite index #####

print(is_finite_index(twist_matrices))