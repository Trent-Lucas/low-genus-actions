from surfaces import *
from lifted_twist import *
from index import *

G = AlternatingGroup(4)
x = G("(1,2)(3,4)"); y = G("(1,2,3)")
number_of_edges = 6
hom = [x,x,x,x,y,y^(-1)]
gluing = {0:1,1:0,2:3,3:2,4:5,5:4}
edge_orientations = {0:1,1:-1,2:1,3:-1,4:1,5:-1}

base_surface = BaseSurface(number_of_edges, gluing, edge_orientations)
cover = Cover(base_surface, G, hom)

h = HomologyGroup(cover)
deck_group_actions = {g:h.action_of_deck_group_on_homology(g) for g in G}
assert h.module.dimension() == 6

##### Building the irrep #####

# Our isotypic component is 2 copies of the standard representation of A_4
# In this case, the isotypic component is the entire homology group

rep = SymmetricGroupRepresentation([3,1], "specht")

# Checking that we have the correct representation
for g in G:
    assert deck_group_actions[g].trace() == 2*rep(g).trace()

##### Building isomorphism between homology and rep #####

A = matrix([[1,0,0,0,1,1],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
assert A.determinant() != 0
B = zero_matrix(6)

for g in G:
    g_bad = deck_group_actions[g].matrix()
    g_good = block_diagonal_matrix(rep(g), rep(g))
    B = B + g_bad*A*(g_good.inverse())

# Checking that we indeed have an isomorphism
assert B.determinant() != 0
for g in G:
    g_bad = deck_group_actions[g].matrix()
    g_good = block_diagonal_matrix(rep(g), rep(g))
    assert g_bad*B*(g_good.inverse()) == B

##### Lifting twists #####

twist_matrices = []

for power in range (1,2):
    for curve in liftable_curves(cover, power):
        T = action_of_twist_on_homology(cover, h, curve, power)
        for g in G:
            assert T.matrix()*deck_group_actions[g].matrix() == deck_group_actions[g].matrix()*T.matrix()
        R = B.inverse()*T.matrix()*B
        T_in_isotypic = matrix([[R[0][0], R[0][3]],[R[3][0],R[3][3]]])
        for i in range(0,2):
            for j in range(0,2):
                assert T_in_isotypic[i][j] in ZZ
        assert T_in_isotypic.determinant() == 1

        twist_matrices.append(T_in_isotypic)
        #print(T_in_isotypic)

##### Check finite index #####

gap.eval("F:=FreeGroup(2)")
gap.eval("a:=F.1")
gap.eval("b:=F.2")
gap.eval("H:=Subgroup(F,[a^2,b])")
print(gap.eval("Index(F,H)"))
print(is_finite_index(twist_matrices) != "infinity")

