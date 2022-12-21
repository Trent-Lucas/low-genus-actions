from surfaces import *
from lifted_twist import *

G = DihedralGroup(6)
number_of_edges = 6
r = G("(1,2,3,4,5,6)"); s = G("(1,6)(2,5)(3,4)")
hom = [s,s,s*r^2,s*r^2,r^3,r^3]
gluing = {0:1,1:0,2:3,3:2,4:5,5:4}
edge_orientations = {0:1,1:-1,2:1,3:-1,4:1,5:-1}

base_surface = BaseSurface(number_of_edges, gluing, edge_orientations)
cover = Cover(base_surface, G, hom)

h = HomologyGroup(cover)
print(h.module.dimension())
action_of_r = h.action_of_deck_group_on_homology(r)
action_of_s = h.action_of_deck_group_on_homology(s)

curve = Curve(cover, [(0,0,1),(1,0,-1),(2,0,1),(3,0,-1)])
T = action_of_twist_on_chain_group(cover, h.edge_chain_group, curve, 3)
for v in h.boundary_kernel.basis():
    w = T(v)
    if w not in h.boundary_kernel:
        print([(v[i],h.edge_chain_group.edge_module.label(h.edge_chain_group.module.lift_map()(h.edge_chain_group.module.basis()[i]))) for i in range(0, len(h.edge_chain_group.module.basis())) if v[i] != 0])
        print([(w[i],h.edge_chain_group.edge_module.label(h.edge_chain_group.module.lift_map()(h.edge_chain_group.module.basis()[i]))) for i in range(0, len(h.edge_chain_group.module.basis())) if w[i] != 0])
        print("----")
#for edge in cover.edge_list:
    #print(edge)
    #print(action_of_twist_on_edge(cover, curve, edge, 1) == action_of_twist_on_edge(cover, curve, cover.gluing[edge], 1))

liftableTwists = []
liftableTwists.append(
    (
        Curve(cover, [(0,0,1),(1,0,-1),(2,0,1),(3,0,-1),(4,0,1),(5,0,-1),(3,1,1),(2,1,-1)]),
        1
    )
)
liftableTwists.append(
    (
        Curve(cover, [(2,0,1),(3,0,-1),(4,0,1),(5,0,-1)]),
        1
    )
)
liftableTwists.append(
    (
        Curve(cover, [(0,0,1),(1,0,-1),(2,0,1),(3,0,-1)]),
        3
    )
)
liftableTwists.append(
    (
        Curve(cover, [(3,0,1),(2,0,-1),(1,0,1),(0,0,-1)]),
        2
    )
)

for pair in liftableTwists:
    T = action_of_twist_on_homology(cover, h, pair[0], pair[1])
    print("done")
    print(T.matrix()*action_of_r.matrix() == action_of_r.matrix()*T.matrix())
    print(T.matrix()*action_of_s.matrix() == action_of_s.matrix()*T.matrix())

print("------------------------------------------")

G = SymmetricGroup(4)
x = G("(1,2)"); y = G("(2,3)"); z = G("(3,4)"); u = G("(1,3)"); v = G("(1,4)"); w = G("(2,4)")
number_of_edges = 6
hom = [x,x, y, y, y*x*z*y, y*x*z*y]
gluing = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4}
edge_orientations = {0:1,1:-1,2:1,3:-1,4:1,5:-1}

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
        Curve(cover, [(0,0,1), (1,0,-1), (2,0,1), (3,0,-1), (4,0,1), (5,0,-1), (3,1,1), (2,1,-1)]),
        1
    )
)
liftableTwists.append(
    (
        Curve(cover, [(2,0,1), (3,0,-1), (4,0,1), (5,0,-1)]),
        2
    )
)
liftableTwists.append(
    (
        Curve(cover, [(0,0,1), (1,0,-1), (2,0,1), (3,0,-1)]),
        3
    )
)

for pair in liftableTwists:
    T = action_of_twist_on_homology(cover, h, pair[0], pair[1])
    print("done")
    print(T.matrix()*action_of_x.matrix() == action_of_x.matrix()*T.matrix())
    print(T.matrix()*action_of_y.matrix() == action_of_y.matrix()*T.matrix())
    print(T.matrix()*action_of_z.matrix() == action_of_z.matrix()*T.matrix())