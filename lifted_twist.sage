from surfaces import *

class Curve:
    """
    A class representing a curve on the base surface and its preimage on the cover.

    Attributes
    ----------
    downstairs_curve: list of 3-tuples (index, position, direction)
        The curve whose twist is being lifted, represented by a sequence of edges on the base surface.
        The index is an integer in 0, ..., cover.number_of_edges labeling the edges on the base surface.
        The position is an integer representing how close the curve is intersecting the edge towards the inital vertex,
        with 0 being the closest.
        The direction is +1 if the curve is entering the edge, -1 if the curve is exiting the edge.
        We assume that the direction of downstairs_curve[0] is +1.
    intersections: dict
        For each edge label i of the base surface, intersection_dict[i] is the list of indices of downstairs_curve
        corresponding to intersections with the ith edge, ordered by position
    length: int
        Length of downstairs_curve.  Represents the number of (doubly counted) edges traversed by each lift of downstairs_curve
    component_coset_rep: dict
        Given an index i of downstairs_curve and a sheet g of the cover, returns a representative of the component of the preimage of
        downstairs curve that lies on sheet g at index i.
    image_under_monodromy: cover.deck_group element
        Viewing downstairs curve as an element of pi_1(base_surface), returns the image of downstairs_curve under the monodromy homomorphism
    """
    def __init__(self, cover, downstairs_curve):
        """
        Parameters
        ----------
        cover: Cover
            The cover containing the curve
        downstairs_curve: list of 3-tuples (index, position, direction)
            The curve whose twist is being lifted, represented by a sequence of edges on the base surface.
            The index is an integer in 0, ..., cover.number_of_edges labeling the edges on the base surface.
            The position is an integer representing how close the curve is intersecting the edge towards the inital vertex,
            with 0 being the closest.
            The direction is +1 if the curve is entering the edge, -1 if the curve is exiting the edge.
            We assume that the direction of downstairs_curve[0] is +1
        """
        if downstairs_curve[0][2] != 1:
            raise ValueError("First entry of downstairs_curve must be an entering edge")

        intersections = {}
        for i in range(0, cover.number_of_edges):
            intersections_with_i = [index for index in range(0, len(downstairs_curve)) if downstairs_curve[index][0] == i]
            intersections_with_i.sort(key=(lambda index : downstairs_curve[index][1]))
            intersections[i] = intersections_with_i

        image_under_monodromy = cover.deck_group.identity()
        for i in range(0, len(downstairs_curve), 2):
            image_under_monodromy = image_under_monodromy * cover.monodromy[downstairs_curve[i][0]]
        
        component_coset_rep = {}
        transveral = [coset[0] for coset in cover.deck_group.cosets(cover.deck_group.subgroup([image_under_monodromy]), side='left')]
        for coset_rep in transveral:
            for j in range(0, image_under_monodromy.order()):
                component_coset_rep[(0, coset_rep*image_under_monodromy^j)] = coset_rep
                component_coset_rep[(len(downstairs_curve)-1, coset_rep*image_under_monodromy^j)] = coset_rep
        current_sheet = cover.monodromy[downstairs_curve[0][0]]

        for i in range(1, len(downstairs_curve)-1, 2):
            for coset_rep in transveral:
                for j in range(0, image_under_monodromy.order()):
                    component_coset_rep[(i, coset_rep*(image_under_monodromy^j)*current_sheet)] = coset_rep
                    component_coset_rep[(i+1, coset_rep*(image_under_monodromy^j)*current_sheet)] = coset_rep
            current_sheet = current_sheet * cover.monodromy[downstairs_curve[i][0]]
        


        self.downstairs_curve = downstairs_curve
        self.intersections = intersections
        self.length = len(downstairs_curve)
        self.image_under_monodromy = image_under_monodromy
        self.component_coset_rep = component_coset_rep


def lift(cover, curve, edge, index, power):
    """
    Returns a sequence of edges representing a lift of downstairs_curve, along with the last edge that the lift enters.

    Given an index of downstairs_curve, returns a pair (lift, entering_edge) were lift is a list of Edges and entering_edge is an Edge.
    Lift is the sequence of edges representing the lift of downstairs_curve at the point given by index on sheet edge.sheet.  The list lift
    will NOT contain edge, nor will it contain the copy of edge at the end of the lift.  entering_edge will be the copy of edge at the end of the lift
    if edge.orientation == 1, and it will be the edge glued to the copy of edge at the end of the lift if edge.orientation == -1.

    Parameters
    ----------
    cover: Cover
        The cover containing curve
    curve: Curve
        The curve being lifted
    edge: Edge
        The edge being acted on by the twist
    index: int
        An index of the list curve.downstairs_curve
    power: int
        The power of the twist being lifted
    """
    downstairs_curve = curve.downstairs_curve
    curve_length = curve.length

    # We assume that index corresponds to edge (in particlar, not the edge glued to edge)
    if edge.index != downstairs_curve[index][0]:
        raise ValueError("index does not match edge")
    
    lift = [] # list of edges to return
    direction = downstairs_curve[index][2] * edge.orientation # determines whether we follow the curve forwards or backwards.  recall that we are computing a left twist

    # because we twist to the left, the initial entering and exiting edges depend on the orientation of the edge
    if edge.orientation == 1:
        index = (index + 2*direction) % curve_length
        exiting_edge = cover.gluing[edge]
        entering_edge = cover.edges[(exiting_edge.sheet, downstairs_curve[index][0])]
    elif edge.orientation == -1:
        index = (index + direction) % curve_length
        exiting_edge = edge
        entering_edge = cover.edges[(edge.sheet, downstairs_curve[index][0])]
    else:
        raise ValueError("edge.orientation must be 1 or -1")
    lift = lift + cover.edges_between(exiting_edge, entering_edge)

    for i in range(0, (curve_length/2)*power - 1):
        index = (index + 2*direction) % curve_length
        exiting_edge = cover.gluing[entering_edge]
        entering_edge = cover.edges[(exiting_edge.sheet, downstairs_curve[index][0])]

        lift = lift + cover.edges_between(exiting_edge, entering_edge)
    
    # if edge is oriented clockwise, the last entering_edge should be a copy of edge on a different sheet 
    # if edge is oritned counterclockwise, the last entering_edge should glued to a copy of edge on a different sheet
    if edge.orientation == 1:
        assert entering_edge.index == edge.index
    else:
        assert cover.gluing[entering_edge].index == edge.index
    # note that this means in any case that the last entering_edge is oriented clockwise, and hence it is not included in the 
    # last call to edges_between

    return (lift, entering_edge)


def action_of_twist_on_edge(cover, curve, edge, power):
    """
    Computes the action of a power of a lifted twist on a single edge.

    Given an edge Edge on a cover Cover, returns a list of pairs (e, i) where e is an edge and i is an element of [1,-1].
    The list of pairs represents the image of edge under the action of the lifted twist around downstairs_curve; the element i is
    +1 or -1 depending on wheter the edge is traversed clockwise or counterclockwise respectively.

    Convention: all twists are assumed to be left twists.

    Parameters
    ----------
    cover: Cover
        The cover being acted on
    curve: Curve
        The curve whose twist is being lifted
    intersection_dict: dict
        For each edge label i of the base surface, intersection_dict[i] is the list of indices of downstairs_curve
        corresponding to intersections with the ith edge, ordered by position
        We pass this in as an argument to avoid recomputing it many times
    edge: Edge
        The edge on the cover being acted on.
    power: integer
        The power of the lifted twist to be computed.
    """
    image = [] # list of pairs (e,i) to be returned

    # First, we find all the indices of downstairs curve corresponding to intersections with edge
    intersections_with_edge = curve.intersections[edge.index]

    # Next, for each intersection, we add a lift of downstairs_curve to image
    current_edge = edge
    if edge.orientation == 1:
        lift_pair = ([], edge)
    else:
        lift_pair = ([], cover.gluing[edge])
    for index in intersections_with_edge:
        lift_pair = lift(cover, curve, current_edge, index, power)
        image = image + lift_pair[0]
        if edge.orientation == 1:
            current_edge = lift_pair[1]
        else:
            current_edge = cover.gluing[lift_pair[1]] # we want to always start on the same side of each copy of edge
    # Then, we have to add in the copy of edge at the end of the last lift
    image.append(lift_pair[1]) # We need to add the last entering edge of the last call of lift to ensure we're tranversing all edges clockwise

    # Finally, to find the deck transformation we post-compose with, we first find the closest intersection of
    # the pre-image of downstairs_curve with the initial vertex of edge in the clockwise direction
    # If the preimage intersects edge and edge is oriented clockwise, we just take the intersection with lowest position
    if len(curve.intersections[edge.index]) > 0 and edge.orientation == 1:
        closest_intersection = curve.intersections[edge.index][0]
    # Otherwise, we start searching for intersections clockwise from the initial vertex of edge
    else:
        edge_index = edge.index
        for i in range(1, cover.number_of_edges+1):
            j = (edge_index+i) % cover.number_of_edges
            if len(curve.intersections[j]) > 0:
                # If the closest edge with intersections is clockwise, we want the intersection with lowest position.
                # If it's oriented counterclockwise, we want the intersection with highest position
                closest_intersection = curve.intersections[j][0] if cover.edges[(edge.sheet,j)].orientation == 1 else curve.intersections[j][-1]
                break
    
    # If the intial vertex of edge lies to the left of the curve, i.e. the closest intersection is entering an edge, the deck transformation is trivial
    if curve.downstairs_curve[closest_intersection][2] == 1:
        deck_transformation = cover.deck_group.identity()
    # Otherwise, we find the deck transformation obtained by lifting downstairs curve to closest_intersection
    else:
        deck_transformation = curve.component_coset_rep[(closest_intersection, edge.sheet)]*curve.image_under_monodromy*curve.component_coset_rep[(closest_intersection, edge.sheet)]^(-1)

    if image[-1].orientation == -1:
        image[-1] = cover.gluing[image[-1]]
    return [cover.edges[((deck_transformation^power)*e.sheet, e.index)] for e in image]


def action_of_twist_on_edge_module(cover, edge_module, curve, power):
    action_on_edges = linear_transformation(
        edge_module.module, 
        edge_module.module,
        lambda e : edge_module.module.sum([edge_module.basis_labels[edge]*edge.orientation for edge in action_of_twist_on_edge(cover, curve, edge_module.label(e), power)])
    )
    return action_on_edges

def action_of_twist_on_chain_group(cover, edge_chain_group, curve, power):
    T = action_of_twist_on_edge_module(cover, edge_chain_group.edge_module, curve, power)
    action_on_chain_group = linear_transformation(
        edge_chain_group.module, 
        edge_chain_group.module,
        lambda e : edge_chain_group.module.quotient_map()(T(edge_chain_group.module.lift_map()(e)))
    )
    return action_on_chain_group

def action_of_twist_on_homology(cover, homology, curve, power):
    T = action_of_twist_on_chain_group(cover, homology.edge_chain_group, curve, power)
    action_on_homology = linear_transformation(
        homology.module, 
        homology.module,
        lambda x : homology.module.quotient_map()(T(homology.module.lift_map()(x)))
    )
    return action_on_homology