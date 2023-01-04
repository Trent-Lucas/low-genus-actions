

# This file was *autogenerated* from the file 3-i-1.sage
from sage.all_cmdline import *   # import sage library

_sage_const_4 = Integer(4); _sage_const_6 = Integer(6); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_3 = Integer(3); _sage_const_5 = Integer(5)#############################################################################
#
# In this file, we show that H_1(S) is 2 copies of the trivial rep and
# two copies of Q(i), and that the action is arithmetic.
#
#############################################################################

from surfaces import *
from lifted_twist import *
from index import *

import lifted_twist

G = CyclicPermutationGroup(_sage_const_4 )
x = G("(1,2,3,4)")
number_of_edges = _sage_const_6 
hom = [x, G.identity(), x**(-_sage_const_1 ), G.identity(), x**_sage_const_2 , x**_sage_const_2 ]
gluing = {_sage_const_0 :_sage_const_2 , _sage_const_1 :_sage_const_3 , _sage_const_2 :_sage_const_0 , _sage_const_3 :_sage_const_1 , _sage_const_4 :_sage_const_5 , _sage_const_5 :_sage_const_4 }
edge_orientations = {_sage_const_0 :_sage_const_1 ,_sage_const_1 :_sage_const_1 ,_sage_const_2 :-_sage_const_1 ,_sage_const_3 :-_sage_const_1 ,_sage_const_4 :_sage_const_1 ,_sage_const_5 :-_sage_const_1 }

base_surface = BaseSurface(number_of_edges, gluing, edge_orientations)
cover = Cover(base_surface, G, hom)

homology = HomologyGroup(cover)
deck_group_actions = {g:homology.action_of_deck_group_on_homology(g).matrix().transpose() for g in G}
# Sage assumes matrices act on the right, so we take the transpose to get an actual homomorphism
for g in G:
    for h in G:
        assert deck_group_actions[g]*deck_group_actions[h] == deck_group_actions[g*h]
assert homology.module.dimension() == _sage_const_6 


##### Building the irrep #####

rep = {x**i: matrix([[_sage_const_0 ,-_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]])**i for i in range(_sage_const_0 ,_sage_const_4 )}

for g in G:
    assert deck_group_actions[g].trace() == _sage_const_2 *rep[g].trace() + _sage_const_2 


##### Building the isomorphism between homology and rep #####

A = matrix([[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ]])
assert A.determinant() != _sage_const_0 
B = zero_matrix(_sage_const_6 )

for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g], identity_matrix(_sage_const_2 ))
    B = B + g_bad*A*(g_good.inverse())

# Checking that we indeed have an isomorphism
assert B.determinant() != _sage_const_0 
for g in G:
    g_bad = deck_group_actions[g]
    g_good = block_diagonal_matrix(rep[g], rep[g], identity_matrix(_sage_const_2 ))
    assert g_bad*B*(g_good.inverse()) == B


##### Lifting twists #####

# This case requires some particularly chosen twists

curve1 = lifted_twist.Curve(cover, [(_sage_const_1 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_3 ,_sage_const_0 ,-_sage_const_1 )])
curve2 = lifted_twist.Curve(cover, [(_sage_const_4 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_5 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_2 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_0 ,_sage_const_1 ,_sage_const_1 ),(_sage_const_2 ,_sage_const_1 ,-_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_3 ,_sage_const_0 ,-_sage_const_1 )])

T1 = action_of_twist_on_homology(cover, homology, curve1, _sage_const_1 )
T2 = action_of_twist_on_homology(cover, homology, curve2, _sage_const_2 )

for g in G:
    assert T1.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T1.matrix().transpose()
    assert T2.matrix().transpose()*deck_group_actions[g] == deck_group_actions[g]*T2.matrix().transpose()

R1 = B.inverse()*T1.matrix().transpose()*B
R2 = B.inverse()*T2.matrix().transpose()*B

# The upper-left 4x4 minor lies in SU(2,2), with complex numbers represented by 2x2 matrices
T1_in_isotypic = matrix([[R1[_sage_const_0 ][_sage_const_0 ]+I*R1[_sage_const_1 ][_sage_const_0 ], R1[_sage_const_0 ][_sage_const_2 ]+I*R1[_sage_const_1 ][_sage_const_2 ]],[R1[_sage_const_2 ][_sage_const_0 ]+I*R1[_sage_const_3 ][_sage_const_0 ],R1[_sage_const_2 ][_sage_const_2 ]+I*R1[_sage_const_3 ][_sage_const_2 ]]])
T2_in_isotypic = matrix([[R2[_sage_const_0 ][_sage_const_0 ]+I*R2[_sage_const_1 ][_sage_const_0 ], R2[_sage_const_0 ][_sage_const_2 ]+I*R2[_sage_const_1 ][_sage_const_2 ]],[R2[_sage_const_2 ][_sage_const_0 ]+I*R2[_sage_const_3 ][_sage_const_0 ],R2[_sage_const_2 ][_sage_const_2 ]+I*R2[_sage_const_3 ][_sage_const_2 ]]])
print(T1_in_isotypic)
print(T2_in_isotypic)

# The above matrices lie in SU(1,1), so we conjugate into SL(2,R)
C = matrix([[_sage_const_1 , -I],[_sage_const_1 ,I]])
twist_matrices = [C.inverse()*T1_in_isotypic*C, C.inverse()*T2_in_isotypic*C]
print(twist_matrices)


##### Check finite index #####

for basis_vector in homology.module.basis():
    edge_vector = homology.edge_chain_group.module.lift_map()(homology.module.lift_map()(basis_vector))
    print([(homology.edge_chain_group.edge_module.label(homology.edge_chain_group.edge_module.module.basis()[i]), edge_vector[i]) for i in range(_sage_const_0 , cover.deck_group.cardinality()*cover.number_of_edges) if edge_vector[i] != _sage_const_0 ])

#print(is_finite_index(twist_matrices))

