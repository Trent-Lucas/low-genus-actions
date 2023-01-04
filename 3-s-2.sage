#############################################################################
#
# In this file, we show that H_1(S) is 1 copy of M_2(Q) and 2 copies of the trivial rep,
# and we show the action is arithmetic
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

import lifted_twist

G = DihedralGroup(4)
r = G("(1,2,3,4)"); s = G("(1,4)(2,3)")
number_of_edges = 6
hom = [s,s*r,s,s*r,r^2,r^2]
gluing = {0:2, 1:3, 2:0, 3:1, 4:5, 5:4}
edge_orientations = {0:1,1:1,2:-1,3:-1,4:1,5:-1}

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
small_rep = {r: 1, s: 1}
for i in range(0, 4):
    for j in range(0,2):
        rep[r^i*s^j] = rep[r]^i * rep[s]^j 
        small_rep[r^i*s^j] = small_rep[r]^i * small_rep[s]^j

for g in G:
    assert deck_group_actions[g].trace() == 2*rep[g].trace() + 2*small_rep[g]


##### Building the isomorphism between homology and rep #####

A = matrix([[1,0,0,0,0,1],[0,1,0,0,0,0],[0,0,1,0,0,1],[0,0,0,1,0,0],[0,0,0,0,1,0],[-1,0,0,0,0,1]])
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

curve1 = lifted_twist.Curve(cover, [(1,0,1),(3,0,-1)])
curve2 = lifted_twist.Curve(cover, [(0,0,1),(2,0,-1)])

T1 = action_of_twist_on_homology(cover, homology, curve1, 2)
T2 = action_of_twist_on_homology(cover, homology, curve2, 2)

for g in G:
    assert T1.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T1.matrix().transpose()
    assert T2.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T2.matrix().transpose()

R1 = B.inverse()*T1.matrix().transpose()*B
R2 = B.inverse()*T2.matrix().transpose()*B

T1_in_isotypic = matrix([[R1[0][0], R1[0][2]],[R1[2][0],R1[2][2]]])
T2_in_isotypic = matrix([[R2[0][0], R2[0][2]],[R2[2][0],R2[2][2]]])

twist_matrices.append(T1_in_isotypic)
twist_matrices.append(T2_in_isotypic)


##### Check finite index #####

print(is_finite_index(twist_matrices))