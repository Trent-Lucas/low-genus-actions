

# This file was *autogenerated* from the file testing.sage
from sage.all_cmdline import *   # import sage library

_sage_const_6 = Integer(6); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_5 = Integer(5)
from surfaces import *
from lifted_twist import *

G = DihedralGroup(_sage_const_6 )
number_of_edges = _sage_const_6 
r = G("(1,2,3,4,5,6)"); s = G("(1,6)(2,5)(3,4)")
hom = [s,s,s*r**_sage_const_2 ,s*r**_sage_const_2 ,r**_sage_const_3 ,r**_sage_const_3 ]
gluing = {_sage_const_0 :_sage_const_1 ,_sage_const_1 :_sage_const_0 ,_sage_const_2 :_sage_const_3 ,_sage_const_3 :_sage_const_2 ,_sage_const_4 :_sage_const_5 ,_sage_const_5 :_sage_const_4 }
edge_orientations = {_sage_const_0 :_sage_const_1 ,_sage_const_1 :-_sage_const_1 ,_sage_const_2 :_sage_const_1 ,_sage_const_3 :-_sage_const_1 ,_sage_const_4 :_sage_const_1 ,_sage_const_5 :-_sage_const_1 }

base_surface = BaseSurface(number_of_edges, gluing, edge_orientations)
cover = Cover(base_surface, G, hom)

h = HomologyGroup(cover)
print(h.module.dimension())
action_of_r = h.action_of_deck_group_on_homology(r)
action_of_s = h.action_of_deck_group_on_homology(s)

curve = Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_2 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_3 ,_sage_const_0 ,-_sage_const_1 )])
T = action_of_twist_on_chain_group(cover, h.edge_chain_group, curve, _sage_const_3 )
for v in h.boundary_kernel.basis():
    w = T(v)
    if w not in h.boundary_kernel:
        print([(v[i],h.edge_chain_group.edge_module.label(h.edge_chain_group.module.lift_map()(h.edge_chain_group.module.basis()[i]))) for i in range(_sage_const_0 , len(h.edge_chain_group.module.basis())) if v[i] != _sage_const_0 ])
        print([(w[i],h.edge_chain_group.edge_module.label(h.edge_chain_group.module.lift_map()(h.edge_chain_group.module.basis()[i]))) for i in range(_sage_const_0 , len(h.edge_chain_group.module.basis())) if w[i] != _sage_const_0 ])
        print("----")
#for edge in cover.edge_list:
    #print(edge)
    #print(action_of_twist_on_edge(cover, curve, edge, 1) == action_of_twist_on_edge(cover, curve, cover.gluing[edge], 1))

liftableTwists = []
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_2 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_3 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_4 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_5 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_3 ,_sage_const_1 ,_sage_const_1 ),(_sage_const_2 ,_sage_const_1 ,-_sage_const_1 )]),
        _sage_const_1 
    )
)
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_2 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_3 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_4 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_5 ,_sage_const_0 ,-_sage_const_1 )]),
        _sage_const_1 
    )
)
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_2 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_3 ,_sage_const_0 ,-_sage_const_1 )]),
        _sage_const_3 
    )
)
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_3 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_2 ,_sage_const_0 ,-_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_0 ,_sage_const_0 ,-_sage_const_1 )]),
        _sage_const_2 
    )
)

for pair in liftableTwists:
    T = action_of_twist_on_homology(cover, h, pair[_sage_const_0 ], pair[_sage_const_1 ])
    print("done")
    print(T.matrix()*action_of_r.matrix() == action_of_r.matrix()*T.matrix())
    print(T.matrix()*action_of_s.matrix() == action_of_s.matrix()*T.matrix())

T = action_of_twist_on_homology(cover, h, Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,-_sage_const_1 )]), _sage_const_1 )
print(T.matrix())
print("------------------------------------------")

print(is_liftable_twist(_sage_const_0 ,_sage_const_2 ,cover,_sage_const_1 ))
print([curve.downstairs_curve for curve in liftable_curves(cover, _sage_const_1 )])
print([curve.downstairs_curve for curve in liftable_curves(cover, _sage_const_2 )])
print([curve.downstairs_curve for curve in liftable_curves(cover, _sage_const_3 )])


print("------------------------------------------")

G = SymmetricGroup(_sage_const_4 )
x = G("(1,2)"); y = G("(2,3)"); z = G("(3,4)"); u = G("(1,3)"); v = G("(1,4)"); w = G("(2,4)")
number_of_edges = _sage_const_6 
hom = [x,x, y, y, y*x*z*y, y*x*z*y]
gluing = {_sage_const_0 :_sage_const_1 , _sage_const_1 :_sage_const_0 , _sage_const_2 :_sage_const_3 , _sage_const_3 :_sage_const_2 , _sage_const_4 :_sage_const_5 , _sage_const_5 :_sage_const_4 }
edge_orientations = {_sage_const_0 :_sage_const_1 ,_sage_const_1 :-_sage_const_1 ,_sage_const_2 :_sage_const_1 ,_sage_const_3 :-_sage_const_1 ,_sage_const_4 :_sage_const_1 ,_sage_const_5 :-_sage_const_1 }

base_surface = BaseSurface(number_of_edges, gluing, edge_orientations)
cover = Cover(base_surface, G, hom)

h = HomologyGroup(cover)
print(h.module.dimension())
action_of_x = h.action_of_deck_group_on_homology(x)
action_of_y = h.action_of_deck_group_on_homology(y)
action_of_z = h.action_of_deck_group_on_homology(z)

liftableTwists = []
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_1 ,_sage_const_0 ,-_sage_const_1 ), (_sage_const_2 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_3 ,_sage_const_0 ,-_sage_const_1 ), (_sage_const_4 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_5 ,_sage_const_0 ,-_sage_const_1 ), (_sage_const_3 ,_sage_const_1 ,_sage_const_1 ), (_sage_const_2 ,_sage_const_1 ,-_sage_const_1 )]),
        _sage_const_1 
    )
)
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_2 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_3 ,_sage_const_0 ,-_sage_const_1 ), (_sage_const_4 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_5 ,_sage_const_0 ,-_sage_const_1 )]),
        _sage_const_2 
    )
)
liftableTwists.append(
    (
        Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_1 ,_sage_const_0 ,-_sage_const_1 ), (_sage_const_2 ,_sage_const_0 ,_sage_const_1 ), (_sage_const_3 ,_sage_const_0 ,-_sage_const_1 )]),
        _sage_const_3 
    )
)

for pair in liftableTwists:
    T = action_of_twist_on_homology(cover, h, pair[_sage_const_0 ], pair[_sage_const_1 ])
    print("done")
    print(T.matrix()*action_of_x.matrix() == action_of_x.matrix()*T.matrix())
    print(T.matrix()*action_of_y.matrix() == action_of_y.matrix()*T.matrix())
    print(T.matrix()*action_of_z.matrix() == action_of_z.matrix()*T.matrix())

T = action_of_twist_on_homology(cover, h, Curve(cover, [(_sage_const_0 ,_sage_const_0 ,_sage_const_1 ),(_sage_const_1 ,_sage_const_0 ,-_sage_const_1 )]), _sage_const_1 )
print(T.matrix())

