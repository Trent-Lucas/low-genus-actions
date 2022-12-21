

# This file was *autogenerated* from the file surfaces.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2)
class BaseSurface:
    """
    A class that assembles the data needed to describe a surface as polygon with edges identified.

    The edges are labeled sequentially by integers 0, .., number_of_edges-1.  The identifications are specified by a dictionary.
    Each edge is given an orientation recorded by +1 or -1 based on whether the edge is oriented clockwise or counterclockwise.

    Attributes
    ----------
    number_of_edges: int 
        the number of edges of the polygon
    gluing: dict 
        for each edge i, gluing[i] = j where j is the edge with which edge i is identified 
    edge_orientations: dict
        for each edge i, edge[i] = 1 if edge i is oriented clockwise, and edge[i] = -1 if edge i is oriented counterclockwise
    """

    def __init__(self, number_of_edges, gluing, edge_orientations):
        """
        Parameters
        ----------
        number_of_edges: int 
            the number of edges of the polygon
        gluing: dict 
            for each edge i, gluing[i] = j where j is the edge with which edge i is identified 
        edge_orientations: dict
            for each edge i, edge[i] = 1 if edge i is oriented clockwise, and edge[i] = -1 if edge i is oriented counterclockwise
        """

        self.number_of_edges = number_of_edges
        self.gluing = gluing
        self.edge_orientations = edge_orientations


class Edge:
    """
    A class for the edges of a cover of a BaseSurface S.

    An edge lies on a particular sheet, i.e. a copy of the polygon describing S labeled by an element of the deck group G.  Within the sheet,
    the edge is labeled by a number in 0, ..., number_of_edges-1.  Each edge is given an orientation recorded by +1 or -1 based on whether the edge is oriented clockwise or counterclockwise.

    Attributes
    ----------
    sheet: group element
        the group element labeling the sheet containing the edge
    index: int
        the label of the edge within its sheet
    orientation: int
        1 if the edge is oriented clockwise, -1 if the edge is oriented counterclockwise
    """

    def __init__(self, sheet, index, orientation):
        """
        Parameters
        ----------
        sheet: group element
            the group element labeling the sheet containing the edge
        index: int
            the label of the edge within its sheet
        orientation: int
            1 if the edge is oriented clockwise, -1 if the edge is oriented counterclockwise
        """

        self.sheet = sheet
        self.index = index
        self.orientation = orientation
    
    def __repr__(self):
        return f"({self.sheet},{self.index},{self.orientation})"


class Cover:
    """
    A class that assembles the data needed for a cover of a BaseSurface S with deck group G and a given monodromy homomorphism.

    The BaseSurface S is represented as a polygon P.  The Cover is obtained by taking |G| copies of P, labeled by the elements of P,
    and gluing the edges according to the monodromy.  Edge i on sheet g is glued to edge S.gluing[i] on sheet g*monodromy[i].

    Attributes
    ----------
    base_surface: BaseSurface:
        the base surface of the cover
    number_of_edges: int
        the number of edges on each sheet of the cover
    edges: dict
        a dictionary containing the edges of the cover.  the keys are (g,i) for g in G and i in 0, ..., number_of_edges-1
    edge_list: list
        the list of values of edges
    deck_group: Group
        the deck group of the cover
    monodromy: dict
        the homomorphism describing the cover. monodromy[i] is the element of G obtained by passing through edge i on P.
    gluing: dict
        for each edge e, contains a pair (e,e') where e' is the edge glued to e.

    Methods
    -------
    edges_between(edge1, edge2):
        Takes in two Edges on the same sheet and returns the list of sheets between their initial vertices, traveling clockwise from edge1 to edge2.
    """

    def __init__(self, base_surface, deck_group, monodromy):
        """
        Parameters
        ----------
        base_surface: BaseSurface
            the base surface of the cover
        deck_group: Group 
            the deck group of the cover
        monodromy: dict
            the homomorphism describing the cover. monodromy[i] is the element of G obtained by passing through edge i on P.
        """
        edges = {(g,i): Edge(g, i, base_surface.edge_orientations[i]) for g in deck_group.list() for i in range(_sage_const_0 , base_surface.number_of_edges)}
        edge_list = list(edges.values())
        gluing = {edge: edges[(edge.sheet*monodromy[edge.index], base_surface.gluing[edge.index])] for edge in edge_list}

        self.base_surface = base_surface
        self.deck_group = deck_group
        self.monodromy = monodromy
        self.number_of_edges = base_surface.number_of_edges
        self.edges = edges
        self.edge_list = edge_list
        self.gluing = gluing

    def edges_between(self, edge1, edge2):
        """Takes in two Edges on the same sheet and returns the list of sheets between their initial vertices, traveling clockwise from edge1 to edge2."""
        if edge1.sheet != edge2.sheet:
            raise ValueError("Edges must lie on the same sheet")

        number_of_edges_between = (edge2.index - edge1.index) % self.number_of_edges
        edge_sequence = []

        if edge1.orientation == _sage_const_1 :
            edge_sequence.append(edge1)
        for i in range(_sage_const_1 , number_of_edges_between):
            edge_sequence.append(self.edges[(edge1.sheet, (edge1.index + i) % self.number_of_edges)])
        if edge2.orientation == -_sage_const_1 :
            edge_sequence.append(edge2)

        return edge_sequence
        


class FreeModuleOnEdges:
    """
    A class for the free vector space spanned by the edges of the sheets of a cover.

    This class is a vector space over Q, togther with a dictionary labeling the basis vectors.  The basis vectors are labeled by edges of the cover.
    In particular, even if two edges are glued together, they represent distinct basis vectors.

    Attributes
    ----------
    module: VectorSpace
        a vector space over Q whose dimension is equal to the total number of edges across all the sheets of the cover.
    basis_labels: dict
        a dictionary assigning a basis vector of module to each edge of the cover.
    cover: Cover
        the cover
    
    Methods
    -------
    label(basis_vector):
        returns the edge labeling a basis vector
    action_of_deck_group_on_edge_module(g):
        for a group element g, returns a linear transformation of module given by the action of g
    """

    def __init__(self, cover):
        module = QQ**(cover.deck_group.cardinality()*cover.number_of_edges)
        basis_labels = {cover.edge_list[i]: module.basis()[i] for i in range(_sage_const_0 ,len(cover.edge_list))}
        
        self.module = module
        self.basis_labels = basis_labels
        self.cover = cover
    
    def label(self, basis_vector):
        """returns the edge labeling a basis vector of module"""
        return list(self.basis_labels.keys())[list(self.basis_labels.values()).index(basis_vector)]

    def action_of_deck_group_on_edge_module(self, g):
        """for a group element g, returns a linear transformation of module given by the action of g"""
        action_on_edges = linear_transformation(
            self.module, 
            self.module,
            lambda e : self.basis_labels[self.cover.edges[(g*self.label(e).sheet, self.label(e).index)]]
        )
        return action_on_edges


class ChainGroupEdges:
    """
    A class for the first cellular chain group of the cover.

    Attributes
    ----------
    edge_module: FreeModuleOnEdges
        the corresponding FreeModuleOnEdges
    module: VectorSpace
        the chain group

    Methods
    -------
    action_of_deck_group_on_edge_chain_group(g):
        for a group element g, returns a linear transformation of module given by the action of g
    """

    def __init__(self, cover):
        edge_module = FreeModuleOnEdges(cover)
        edge_relations = [edge_module.basis_labels[edge] - edge_module.basis_labels[cover.gluing[edge]] for edge in cover.edge_list]
        module = edge_module.module/(edge_module.module.span(edge_relations))

        self.edge_module = edge_module
        self.module = module
    
    def action_of_deck_group_on_edge_chain_group(self, g):
        action_on_edges = self.edge_module.action_of_deck_group_on_edge_module(g)
        action_on_chain_group = linear_transformation(
            self.module, 
            self.module,
            lambda e : self.module.quotient_map()(action_on_edges(self.module.lift_map()(e)))
        )
        return action_on_chain_group


class FreeModuleOnVertices:
    """
    A class for the free module on the vertices of the cover.

    A vertex is a pair (edge, i) where edge is an Edge and i is 0 or 1, representing the initial and terminal vertices respectively.
    This class is a vector space over Q, togther with a dictionary labeling the basis vectors.  The basis vectors are labeled by the vertices of each edge of the cover.
    In particular, even vertices of sequential edges or glued edges are considered distinct.

    Attributes
    ----------
    module: VectorSpace
        a vector space over Q whose dimension is equal to the total number of vertices across all the sheets of the cover.
    basis_labels: dict
        a dictionary assigning a basis vector of module to each vertex of the cover.
    
    Methods
    -------
    label(basis_vector):
        returns the vertex labeling a basis vector
    """

    def __init__(self, cover):
        module = free_module = QQ**(cover.deck_group.cardinality()*cover.number_of_edges*_sage_const_2 )
        vertex_list = [(edge, i) for edge in cover.edge_list for i in [_sage_const_0 ,_sage_const_1 ]]
        basis_labels = {vertex_list[i]: free_module.basis()[i] for i in range(_sage_const_0 ,len(vertex_list))}

        self.module = module
        self.basis_labels = basis_labels
    
    def label(self, basis_vector):
        """returns the vertex (edge, i) labeling a basis vector of module"""
        return list(self.basis_labels.keys())[list(self.basis_labels.values()).index(basis_vector)]



class ChainGroupVertices:
    """
    A class for the zeroth cellular chain group of the cover.

    Attributes
    ----------
    vertex_module: FreeModuleOnVertices
        the corresponding FreeModuleOnVertices
    module: VectorSpace
        the chain group
    """

    def __init__(self, cover):
        vertex_module = FreeModuleOnVertices(cover)

        vertex_relations = []
        for edge in cover.edge_list:
            # First, we identify the vertices along the gluing
            glued_edge = cover.gluing[edge]
            vertex_relations.append(vertex_module.basis_labels[(edge,_sage_const_0 )] - vertex_module.basis_labels[(glued_edge,_sage_const_0 )])
            vertex_relations.append(vertex_module.basis_labels[(edge,_sage_const_1 )] - vertex_module.basis_labels[(glued_edge,_sage_const_1 )])

            # Next, we identify the vertices of sequential edges
            next_edge = cover.edges[(edge.sheet, (edge.index + edge.orientation) % cover.number_of_edges)]
            previous_edge = cover.edges[(edge.sheet, (edge.index - edge.orientation) % cover.number_of_edges)]
            if edge.orientation == next_edge.orientation:
                vertex_relations.append(vertex_module.basis_labels[(edge,_sage_const_1 )] - vertex_module.basis_labels[(next_edge,_sage_const_0 )])
            else:
                vertex_relations.append(vertex_module.basis_labels[(edge,_sage_const_1 )] - vertex_module.basis_labels[(next_edge,_sage_const_1 )])
            if edge.orientation == previous_edge.orientation:
                vertex_relations.append(vertex_module.basis_labels[(edge,_sage_const_0 )] - vertex_module.basis_labels[(previous_edge,_sage_const_1 )])
            else:
                vertex_relations.append(vertex_module.basis_labels[(edge,_sage_const_0 )] - vertex_module.basis_labels[(previous_edge,_sage_const_0 )])

        module = vertex_module.module/(vertex_module.module.span(vertex_relations))
        
        self.vertex_module = vertex_module
        self.module = module


class HomologyGroup:
    """
    A class for the first homology group of the cover.

    Attributes
    ----------
    module: VectorSpace
        the homology group
    edge_chain_group: ChainGroupEdges
        the first cellular chain group of the cover
    vertex_chain_group: ChainGroupVertices
        the zeroth cellular chain group of the cover
    boundary_kernel: VectorSpace
        the kernel of the boundary map edge_chain_group -> vertex_chain_group
    boundary_map: linear_transformation
        the boundary map edge_chain_group -> vertex_chain_group
    
    Methods
    -------
    action_of_deck_group_on_homology(g):
        for a group element g, returns a linear transformation of module given by the action of g
    """

    def __init__(self, cover):
        edge_chain_group = ChainGroupEdges(cover)
        edge_module = edge_chain_group.edge_module
        vertex_chain_group = ChainGroupVertices(cover)
        vertex_module = vertex_chain_group.vertex_module
        
        # Our first goal is to build the kernel of the boundary map
        # First, we build the boundary map from the free module on edges to vertex chain group
        boundary_map_edge_module = linear_transformation(
            edge_module.module, 
            vertex_module.module, 
            lambda e : vertex_module.basis_labels[(edge_module.label(e),_sage_const_1 )] - vertex_module.basis_labels[(edge_module.label(e),_sage_const_0 )]
        )

        # Next, we descend this map to the chain groups
        boundary_map_chain_group = linear_transformation(
            edge_chain_group.module,
            vertex_chain_group.module,
            lambda x : vertex_chain_group.module.quotient_map()(boundary_map_edge_module(edge_chain_group.module.lift_map()(x)))
        )
        boundary_kernel = boundary_map_chain_group.kernel()

        # Next, we build the image of the boundary map B1
        face_relations = [
            sum(
                [edge_module.basis_labels[cover.edges[(g,i)]]*cover.edges[(g,i)].orientation for i in range(_sage_const_0 ,cover.number_of_edges)]
            ) 
            for g in cover.deck_group.list()
        ]
        boundary_image = edge_chain_group.module.span([edge_chain_group.module.quotient_map()(relation) for relation in face_relations])

        module = boundary_kernel/boundary_image

        self.edge_chain_group = edge_chain_group
        self.vertex_chain_group = vertex_chain_group
        self.module = module
        self.boundary_kernel = boundary_kernel
        self.boundary_map = boundary_map_chain_group
    
    def action_of_deck_group_on_homology(self, g):
        action_on_chain_group = self.edge_chain_group.action_of_deck_group_on_edge_chain_group(g)
        action_on_homology = linear_transformation(
            self.module, 
            self.module,
            lambda x : self.module.quotient_map()(action_on_chain_group(self.module.lift_map()(x)))
        )
        return action_on_homology


