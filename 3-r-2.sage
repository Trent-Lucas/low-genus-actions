#############################################################################
#
# In this file, we show that H_1(S) is 1 copy of M_2(Q) and 2 copies of a 1-dimensional irrep,
# and we show the action is arithmetic
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

import lifted_twist

G = DihedralGroup(4)
r = G("(1,2,3,4)"); s = G("(1,4)(2,3)")
number_of_edges = 8
hom = [s,s,s,s,s*r,s*r,s*r^(-1),s*r^(-1)]
gluing = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4, 6:7, 7:6}
edge_orientations = {0:1,1:-1,2:1,3:-1,4:1,5:-1, 6:1, 7:-1}

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

rep = {r: matrix([[0, -1], [1, 0]]), s: matrix([[0,1],[1,0]])}
small_rep = {r: 1, s: -1}
for i in range(0, 4):
    for j in range(0,2):
        rep[r^i*s^j] = rep[r]^i * rep[s]^j 
        small_rep[r^i*s^j] = small_rep[r]^i * small_rep[s]^j

for g in G:
    assert deck_group_actions[g].trace() == 2*rep[g].trace() + 2*small_rep[g]


##### Building the isomorphism between homology and rep #####

A = matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
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

curves = liftable_curves(cover, 1)
curves.append(lifted_twist.Curve(cover,[(7,0,1),(6,0,-1),(5,0,1),(4,0,-1),(3,0,1),(2,0,-1)]))
curves.append(lifted_twist.Curve(cover,[(1,0,1),(0,0,-1),(7,0,1),(6,0,-1),(5,0,1),(4,0,-1)]))
curves.append(lifted_twist.Curve(cover,[(3,0,1),(2,0,-1),(1,0,1),(0,0,-1),(7,0,1),(6,0,-1)]))
curves.append(lifted_twist.Curve(cover,[(5,0,1),(4,0,-1),(3,0,1),(2,0,-1),(1,0,1),(0,0,-1)]))

for curve in curves:
    T = action_of_twist_on_homology(cover, homology, curve, 2)
    # Can verify that T commutes with deck group
    for g in G:
        assert T.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T.matrix().transpose()
    
    R = B.inverse()*T.matrix().transpose()*B
    # The top-left 4x4 minor of R will be a block matrix where each block is a multiple of the identity
    T_in_isotypic = matrix([[R[0][0], R[0][2]],[R[2][0],R[2][2]]])
    
    # Can verify that T has integer entries
    for i in range(0,2):
        for j in range(0,2):
            assert T_in_isotypic[i][j] in ZZ
    assert T_in_isotypic.determinant() == 1

    twist_matrices.append(T_in_isotypic)


##### Check finite index #####

print(is_finite_index(twist_matrices))