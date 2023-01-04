#############################################################################
#
# In this file, we show that H_1(S) is 2 copies of the trivial rep and
# two copies of Q(i), and that the action is arithmetic.
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

import lifted_twist

G = CyclicPermutationGroup(4)
x = G("(1,2,3,4)")
number_of_edges = 6
hom = [x, G.identity(), x^(-1), G.identity(), x^2, x^2]
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

rep = {x^i: matrix([[0,-1],[1,0]])^i for i in range(0,4)}

for g in G:
    assert deck_group_actions[g].trace() == 2*rep[g].trace() + 2


##### Building the isomorphism between homology and rep #####

A = matrix([[1,0,0,0,1,1],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
assert A.determinant() != 0
B = zero_matrix(6)

for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g], identity_matrix(2))
    B = B + g_bad*A*(g_good.inverse())

# Checking that we indeed have an isomorphism
assert B.determinant() != 0
for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g], identity_matrix(2))
    assert g_bad*B*(g_good.inverse()) == B


##### Lifting twists #####

# This case requires some particularly chosen twists

curve1 = lifted_twist.Curve(cover, [(1,0,1),(3,0,-1)])
curve2 = lifted_twist.Curve(cover, [(4,0,1),(5,0,-1),(0,0,1),(2,0,-1),(0,1,1),(2,1,-1),(1,0,1),(3,0,-1)])

T1 = action_of_twist_on_homology(cover, homology, curve1, 1)
T2 = action_of_twist_on_homology(cover, homology, curve2, 2)

for g in G:
    assert T1.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T1.matrix().transpose()
    assert T2.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T2.matrix().transpose()

R1 = B.inverse()*T1.matrix().transpose()*B
R2 = B.inverse()*T2.matrix().transpose()*B

# The upper-left 4x4 minor lies in the unitary group of the Reidemeister pairing, with complex numbers represented by 2x2 matrices
T1_in_isotypic = matrix([[R1[0][0]+I*R1[1][0], R1[0][2]+I*R1[1][2]],[R1[2][0]+I*R1[3][0],R1[2][2]+I*R1[3][2]]])
T2_in_isotypic = matrix([[R2[0][0]+I*R2[1][0], R2[0][2]+I*R2[1][2]],[R2[2][0]+I*R2[3][0],R2[2][2]+I*R2[3][2]]])


##### Putting the intersection pairing into standard form ######

# The following was used to see explicitly our homology basis
#for basis_vector in homology.module.basis():
#    edge_vector = homology.edge_chain_group.module.lift_map()(homology.module.lift_map()(basis_vector))
#    print([(homology.edge_chain_group.edge_module.label(homology.edge_chain_group.edge_module.module.basis()[i]), edge_vector[i]) for i in range(0, cover.deck_group.cardinality()*cover.number_of_edges) if edge_vector[i] != 0])

# The above gives the following intersection matrix in the given basis of H_1(S)
intersection_matrix = matrix([[0,0,0,0,-1,0],[0,0,1,0,1,-1],[0,-1,0,0,0,0],[0,0,0,0,0,-1],[1,-1,0,0,0,1],[0,1,0,1,-1,0]])
assert intersection_matrix.determinant() != 0
assert intersection_matrix.transpose() == -intersection_matrix
assert T1.matrix().transpose().transpose()*intersection_matrix*T1.matrix().transpose() == intersection_matrix
assert T2.matrix().transpose().transpose()*intersection_matrix*T2.matrix().transpose() == intersection_matrix

# Now, we change to our preferred basis and compute the Reidemeister pairing
S = B.transpose()*intersection_matrix*B

reidemeister_pairing = matrix([[2*(S[2*i][2*j] + I*S[2*i][2*j+1]) for j in range(0,2)] for i in range(0,2)])
assert reidemeister_pairing.determinant() != 0
assert reidemeister_pairing.transpose() == -1*reidemeister_pairing.conjugate()
assert T1_in_isotypic.transpose()*reidemeister_pairing*T1_in_isotypic.conjugate() == reidemeister_pairing
assert T2_in_isotypic.transpose()*reidemeister_pairing*T2_in_isotypic.conjugate() == reidemeister_pairing

# In our current basis, the Reidemeister pairing is the matrix [[16i, -8],[8,0]]
# We do one more change of basis to land in SL(2,Z)
P = matrix([[1,0],[-I,1]])
T1_in_isotypic_proper_basis = P.inverse()*T1_in_isotypic*P
T2_in_isotypic_proper_basis = P.inverse()*T2_in_isotypic*P
assert T1_in_isotypic_proper_basis.determinant() == 1
assert T2_in_isotypic_proper_basis.determinant() == 1
for i in range(0,2):
    for j in range(0,2):
        assert T1_in_isotypic_proper_basis[i][j] in ZZ
        assert T2_in_isotypic_proper_basis[i][j] in ZZ

twist_matrices = [matrix(ZZ,T1_in_isotypic_proper_basis), matrix(ZZ,T2_in_isotypic_proper_basis)]


##### Checking finite index #####

print(is_finite_index(twist_matrices))