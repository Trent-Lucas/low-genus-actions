#Input: a finite group G, an integer number_of_edges
#Attributes: module is a Q-vector space of dimension |G|*number_of_edges,
#            edges is a list of edges (g,i)
#            basis is a dictionary mapping edges (g,i) to a basis of module
class FreeModuleOnEdges:
    def __init__(self, G, number_of_edges):
        module = QQ^(G.cardinality()*number_of_edges)
        edges = [(g, i) for g in G.list() for i in range(0,number_of_edges)]
        basisLabels = {}
        for i in range(0, len(edges)):
            basisLabels[edges[i]] = module.basis()[i]
        
        self.module = module
        self.edges = edges
        self.basisLabels = basisLabels
    
    #For a basis vector v of module, returns the corresponding edge (g,i)
    def edge_label(self, v):
        return list(self.basisLabels.keys())[list(self.basisLabels.values()).index(v)]

    #Returns the action of a group element g on module
    def action_of_G_on_edge_module(self, g):
        action_on_edges = linear_transformation(
            self.module, 
            self.module,
            lambda e : self.basisLabels[(g*self.edge_label(e)[0], self.edge_label(e)[1])]
        )
        return action_on_edges

#Input: a finite group G, an integer number_of_edges
#Attributes: module is a Q-vector space of dimension |G|*number_of_edges*2,
#            vertices is a list of vertices ((g,i),j)
#            basis is a dictionary mapping edges ((g,i),j) to a basis of module
class FreeModuleOnVertices:
    def __init__(self, G, number_of_edges):
        module = QQ^(G.cardinality()*number_of_edges*2)
        vertices = [((g, i),j) for g in G.list() for i in range(0,number_of_edges) for j in [0,1]]
        basisLabels = {}
        for i in range(0, len(vertices)):
            basisLabels[vertices[i]] = module.basis()[i]

        self.module = module
        self.vertices = vertices
        self.basisLabels = basisLabels

#Input: 
#   G is a finite group,
#   number_of_edges is an integer representing the number of edges,
#   gluing is a dictionary representing which edges are glued together,
#   hom is a list of elements of G, where the ith entry is the deck transformation of
#       the ith edge
#Attributes:
#   edge_module is the free module on the edges 
#   module is the chain group C1
class ChainGroupEdges:
    def __init__(self, G, number_of_edges, gluing, hom):
        edge_module = FreeModuleOnEdges(G, number_of_edges)
        edgeRelations = []
        for g in G.list():
            for i in range(0,number_of_edges):
                edgeRelations.append(edge_module.basisLabels[(g,i)] - edge_module.basisLabels[(g*hom[i], gluing[i])])
        module = edge_module.module/(edge_module.module.span(edgeRelations))

        self.edge_module = edge_module
        self.module = module
    
    def action_of_G_on_edge_chain_group(self, g):
        action_on_edges = self.edge_module.action_of_G_on_edge_module(g)
        action_on_chain_group = linear_transformation(
            self.module, 
            self.module,
            lambda e : self.module.quotient_map()(action_on_edges(self.module.lift_map()(e)))
        )
        return action_on_chain_group

#Input: 
#   G is a finite group,
#   number_of_edges is an integer representing the number of edges,
#   gluing is a dictionary representing which edges are glued together,
#   hom is a list of elements of G, where the ith entry is the deck transformation of
#       the ith edge
#Attributes:
#   vertex_module is the free module on the vertices 
#   module is the chain group C0
class ChainGroupVertices:
    def __init__(self, G, number_of_edges, gluing, hom, orientation):
        vertex_module = FreeModuleOnVertices(G, number_of_edges)
        vertex_relations = []
        for g in G.list():
            for i in range(0,number_of_edges):
                vertex_relations.append(vertex_module.basisLabels[((g,i),0)] - vertex_module.basisLabels[((g*hom[i], gluing[i]),0)])
                vertex_relations.append(vertex_module.basisLabels[((g,i),1)] - vertex_module.basisLabels[((g*hom[i], gluing[i]),1)])
                edge_at_head = (i + orientation[i]) % number_of_edges
                edge_at_tail = (i - orientation[i]) % number_of_edges
                if orientation[i]*orientation[edge_at_head] == 1:
                    vertex_relations.append(vertex_module.basisLabels[((g,i),1)] - vertex_module.basisLabels[((g,edge_at_head),0)])
                else:
                    vertex_relations.append(vertex_module.basisLabels[((g,i),1)] - vertex_module.basisLabels[((g,edge_at_head),1)])
                if orientation[i]*orientation[edge_at_tail] == 1:
                    vertex_relations.append(vertex_module.basisLabels[((g,i),0)] - vertex_module.basisLabels[((g,edge_at_tail),1)])
                else:
                    vertex_relations.append(vertex_module.basisLabels[((g,i),0)] - vertex_module.basisLabels[((g,edge_at_tail),0)])
        module = vertex_module.module/(vertex_module.module.span(vertex_relations))

        self.vertex_module = vertex_module
        self.module = module

#Input: 
#   G is a finite group,
#   number_of_edges is an integer representing the number of edges,
#   gluing is a dictionary representing which edges are glued together,
#   hom is a list of elements of G, where the ith entry is the deck transformation of
#       the ith edge
#   edgeOrientations is a dictionary mapping i to [1,-1]
class HomologyGroup:
    def __init__(self, G, number_of_edges, gluing, hom, orientation):
        edge_chain_group = ChainGroupEdges(G, number_of_edges, gluing, hom)
        edge_module = edge_chain_group.edge_module
        vertex_chain_group = ChainGroupVertices(G, number_of_edges, gluing, hom, orientation)
        vertex_module = vertex_chain_group.vertex_module

        boundary_map_edge_module = linear_transformation(
            edge_module.module, 
            vertex_module.module, 
            lambda e : vertex_module.basisLabels[(edge_module.edge_label(e),1)] - vertex_module.basisLabels[(edge_module.edge_label(e),0)]
        )
        boundary_map_chain_group = linear_transformation(
            edge_chain_group.module,
            vertex_chain_group.module,
            lambda x : vertex_chain_group.module.quotient_map()(boundary_map_edge_module(edge_chain_group.module.lift_map()(x)))
        )
        Z1 = boundary_map_chain_group.kernel()

        faceRelations = [sum([edge_module.basisLabels[(g,i)]*orientation[i] for i in range(0,number_of_edges)]) for g in G.list()]
        B1 = edge_chain_group.module.span([edge_chain_group.module.quotient_map()(relation) for relation in faceRelations])

        module = Z1/B1

        self.edge_chain_group = edge_chain_group
        self.vertex_chain_group = vertex_chain_group
        self.module = module

    def action_of_G_on_homology(self, g):
        action_on_chain_group = self.edge_chain_group.action_of_G_on_edge_chain_group(g)
        action_on_homology = linear_transformation(
            self.module, 
            self.module,
            lambda x : self.module.quotient_map()(action_on_chain_group(self.module.lift_map()(x)))
        )
        return action_on_homology

