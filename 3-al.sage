#############################################################################
#
# In this file, we show that H_1(S) is two copies of the standard rep tensor the sign rep, 
# and that the action is arithmetic.
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

G = SymmetricGroup(4)
x = G("(1,2)"); y = G("(2,3)"); z = G("(3,4)"); u = G("(1,3)"); v = G("(1,4)"); w = G("(2,4)")
number_of_edges = 6
hom = [x,x, y, y, y*x*z*y, y*x*z*y]
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

# Our isotypic component is 2 copies of the standard representation tensor the sign representation
# In this case, the isotypic component is the entire homology group
# We want an orthogonal representation with integer entries, so we build it by hand.

# First, we define the representation on the transpositions
rep = {x: (-1)*matrix([[0,1,0],[1,0,0],[0,0,1]]), y: (-1)*matrix([[0,0,-1],[0,1,0],[-1,0,0]]), z: (-1)*matrix([[0,-1,0],[-1,0,0],[0,0,1]]),
    u: (-1)*matrix([[1,0,0],[0,0,-1],[0,-1,0]]), v: (-1)*matrix([[0,0,1],[0,1,0],[1,0,0]]), w: (-1)*matrix([[1,0,0],[0,0,1],[0,1,0]])}

# Next, given an arbitrary permutation, we wish to decompose it into transpositions
def transposition_decomp(g):
    transposition_tuples = []
    for t in g.cycle_tuples():
        for i in range(0, len(t) - 1):
            transposition_tuples.append((t[i], t[i+1]))
    
    # We can check that our decomposition is correct
    product = G.identity()
    for transposition in transposition_tuples:
        product = G(transposition)*product
    assert product == g

    return transposition_tuples

# Now, we can build our representation
def get_rep(g):
    transpositions = transposition_decomp(g)
    image = matrix([[1,0,0],[0,1,0],[0,0,1]])
    for transposition in transpositions:
        image = rep[G(transposition)]*image
    return image

# To minimize computation, we store the representation in a dictionary
for g in G.list():
    rep[g] = get_rep(g)

# Checking that we have the correct representation
for g in G:
    assert deck_group_actions[g].trace() == 2*rep[g].trace()


##### Building isomorphism between homology and rep #####

A = matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
assert A.determinant() != 0
B = zero_matrix(6)

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
        T_in_isotypic = matrix([[R[0][0], R[0][3]],[R[3][0],R[3][3]]])
        
        # Can verify that T has integer entries
        for i in range(0,2):
            for j in range(0,2):
                assert T_in_isotypic[i][j] in ZZ
        assert T_in_isotypic.determinant() == 1

        twist_matrices.append(T_in_isotypic)


##### Check finite index #####

print(is_finite_index(twist_matrices))