from homology import *

#Helper function: takes in two edges on the same sheet (g,i) and (g,j) and returns list of oriented edges (g, k, direction) between their inital vertices going clockwise
def edges_between(edge1, edge2, number_of_edges, edge_orientations):
    g = edge1[0]
    if edge2[0] != g:
        raise ValueError("Edges must lie on the same sheet")

    fullEdgesBetween = (edge2[1] - edge1[1]) % number_of_edges

    edges = []
    if edge_orientations[edge1[1]] == 1:
        edges.append((g,edge1[1],1))
    for i in range(edge1[1]+1, edge1[1] + fullEdgesBetween):
        edges.append((g, i % number_of_edges, 1))
    if edge_orientations[edge2[1]] == -1:
        edges.append((g,edge2[1],1))
    
    return edges


# Input:
#   G is a finite group
#   hom is a list of elements of G 
#   downstairs curve is a list of 3-tuples (i, position, orientation) where i is the edge label, indicates how close the intersecting arc is to 
#       the initial vertex (0 is closest), orientation is an element of [1,-1] indicating whether the arc is entering or exiting
#   e is an edge (g,i)
#   power is an integer
# Output: A list of tuples (g, i, direction) where direction is 1 if edge is traversed clockwise, -1 if edge is traversed counterclockwise
#
# Idea: 

def action_of_dehn_twist_on_edge(G, number_of_edges, hom, gluing, edge_orientations, downstairs_curve, e, power):
    liftToIdentitySheet = []
    currentSheet = G.identity()
    for edge in downstairs_curve:
        liftToIdentitySheet.append((currentSheet, edge[0], edge[1], edge[2]))
        if edge[2] == -1:
            continue
        else:
            currentSheet = currentSheet * hom[edge[0]]
    deckTransformation = currentSheet

    initialLift = [edge for edge in liftToIdentitySheet]
    for i in range(1, deckTransformation.order()):
        for edge in initialLift:
            liftToIdentitySheet.append((deckTransformation^(i) * edge[0], edge[1], edge[2], edge[3]))

    transveral = [coset[0] for coset in G.cosets(G.subgroup([deckTransformation]), side='left')]
    lifts = [[(r*edge[0], edge[1], edge[2], edge[3]) for edge in liftToIdentitySheet] for r in transveral]

    #intersections_with_e is a list of pairs (i,j) where lifts[i][j] == e
    intersections_with_e = []
    for i in range(0, len(lifts)):
        for j in range(0, len(lifts[i])):
            edge = lifts[i][j]
            if (edge[0], edge[1]) == (e[0],e[1]):
                intersections_with_e.append((i, j))
    #we sort intersections_with_e by position
    intersections_with_e.sort(key=(lambda pair : lifts[pair[0]][pair[1]][2]))

    if len(intersections_with_e) == 0:
        # To figure out "correcion" deck transformation, need to find out which component the initial vertex of e lies on 
        # Enough to figure out which lift is "closest" as you move clockwise around the boundary 
        
        if edge_orientations[(e[1] - 1) % number_of_edges] == 1:
            closestEdge = (e[0], (e[1] - 1) % number_of_edges, float('inf'), 1)
        else:
            closestEdge = (e[0], (e[1] - 1) % number_of_edges, 0, 1)

        for i in range(0, len(lifts)):
            for j in range(0, len(lifts[i])):
                edge = lifts[i][j]
                if edge[0] == e[0]:
                    if ((edge[1] - e[1]) % number_of_edges) < ((closestEdge[1] - e[1]) % number_of_edges):
                        closestPair = [i,j]
                        closestEdge = (edge[0],edge[1],edge[2],edge[3])
                    
                    if ((edge[1] - e[1]) % number_of_edges) == ((closestEdge[1] - e[1]) % number_of_edges):
                        if edge_orientations[edge[1]] == 1 and edge[3] <= closestEdge[3]:
                            closestPair = [i,j]
                            closestEdge = (edge[0], edge[1], edge[2], edge[3])
                        elif edge_orientations[edge[1]] == -1 and edge[3] >= closestEdge[3]:
                            closestPair = [i,j]
                            closestEdge = (edge[0], edge[1], edge[2], edge[3])

        if lifts[closestPair[0]][closestPair[1]][3] == 1:
            correctionDeckTransformation = G.identity()
        else:
            correctionDeckTransformation = transveral[closestPair[0]]*deckTransformation*transveral[closestPair[0]]^(-1)

        return [((correctionDeckTransformation^power)*e[0],e[1],edge_orientations[e[1]])]
    
    else:
        newEdges = []
        offset = 0 #keeps track of total deck transformation accumulated so far

        # To figure out "correction" deck transformation, need to determine the deck transformation for the component of initial vertex
        # Assumption: we fix the component to the left of lifts[0]
        # Given this, you just need to figure out if you're to the left or right of lifts[intersections_with_e[0][0]]
        firstPair = intersections_with_e[0]
        if edge_orientations[e[1]]*lifts[firstPair[0]][firstPair[1]][3] == 1:
            correctionDeckTransformation = G.identity()
        else:
            correctionDeckTransformation = transveral[firstPair[0]]*deckTransformation*transveral[firstPair[0]]^(-1)

        if edge_orientations[e[1]] == 1:
            for pair in intersections_with_e:
                givenLift = lifts[pair[0]]
                orientation = givenLift[pair[1]][3]
                index = (pair[1] + offset*len(initialLift)*orientation) % len(givenLift)
                currentSheet = givenLift[index][0]
                enteringEdge = (givenLift[index][0], givenLift[index][1])
                offset = offset + power

                for i in range(0, (len(initialLift)/2)*power):
                    currentSheet = currentSheet * hom[enteringEdge[1]]
                    index = (index + (2*orientation)) % len(givenLift)
                    exitingEdge = (currentSheet, gluing[enteringEdge[1]])
                    enteringEdge = (givenLift[index][0], givenLift[index][1])
                    
                    newEdges = newEdges + edges_between(exitingEdge, enteringEdge, number_of_edges, edge_orientations)
            newEdges.append((currentSheet, e[1], 1))

            return [((correctionDeckTransformation^power)*edge[0],edge[1],edge[2]) for edge in newEdges]

        if edge_orientations[e[1]] == -1:
            for pair in intersections_with_e:
                givenLift = lifts[pair[0]]
                orientation = -givenLift[pair[1]][3] #will pick up positive power of deck transform if givenLift is leaving this edge
                index = (pair[1] + offset*len(initialLift)*orientation) % len(givenLift)
                currentSheet = givenLift[index][0]
                exitingEdge = (givenLift[index][0], givenLift[index][1])
                offset = offset + power

                index = (index + orientation) % len(givenLift)
                enteringEdge = (givenLift[index][0],givenLift[index][1])
                newEdges = newEdges + edges_between(exitingEdge, enteringEdge, number_of_edges, edge_orientations)                

                for i in range(0, ((len(initialLift)/2)*power) - 1):
                    currentSheet = currentSheet * hom[enteringEdge[1]]
                    index = (index + (2*orientation)) % len(givenLift)
                    exitingEdge = (currentSheet, gluing[enteringEdge[1]])
                    enteringEdge = (givenLift[index][0], givenLift[index][1])
                    
                    newEdges = newEdges + edges_between(exitingEdge, enteringEdge, number_of_edges, edge_orientations)
            
            currentSheet = currentSheet * hom[enteringEdge[1]]
            newEdges.append((currentSheet, e[1], -1))

            return [((correctionDeckTransformation^power)*edge[0],edge[1], edge[2]) for edge in newEdges]


def action_of_dehn_twist_on_edge_module(G, number_of_edges, hom, gluing, downstairs_curve, edge_module, edge_orientations,power):
    edgeLabel = lambda e : list(edge_module.basisLabels.keys())[list(edge_module.basisLabels.values()).index(e)]
    actionOnEdges = linear_transformation(
        edge_module.module, 
        edge_module.module,
        lambda e : edge_module.module.sum([edge_module.basisLabels[(edge[0],edge[1])]*edge[2]*edge_orientations[edge[1]] for edge in action_of_dehn_twist_on_edge(G, number_of_edges, hom, gluing, edge_orientations, downstairs_curve, edgeLabel(e),power)])
    )
    return actionOnEdges


def action_of_dehn_twist_on_chain_group(G, number_of_edges, hom, gluing, downstairs_curve, edge_chain_group, edge_orientations,power):
    T = action_of_dehn_twist_on_edge_module(G, number_of_edges, hom, gluing, downstairs_curve, edge_chain_group.edge_module, edge_orientations,power)
    actionOnChainGroup = linear_transformation(
        edge_chain_group.module, 
        edge_chain_group.module,
        lambda e : edge_chain_group.module.quotient_map()(T(edge_chain_group.module.lift_map()(e)))
    )
    return actionOnChainGroup

def action_of_dehn_twist_on_homology(G, number_of_edges, hom, gluing, downstairs_curve, homology, edge_orientations,power):
    T = action_of_dehn_twist_on_chain_group(G, number_of_edges, hom, gluing, downstairs_curve, homology.edge_chain_group, edge_orientations,power)
    actionOnHomology = linear_transformation(
        homology.module, 
        homology.module,
        lambda x : homology.module.quotient_map()(T(homology.module.lift_map()(x)))
    )
    return actionOnHomology

###Testing###

#G = SymmetricGroup(4)
#x = G("(1,2)"); y = G("(2,3)"); z = G("(3,4)"); u = G("(1,3)"); v = G("(1,4)"); w = G("(2,4)")
#number_of_edges = 6
#hom = [x,x, y, y, y*x*z*y, y*x*z*y]
#gluing = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4}
#edge_orientations = {0:1,1:-1,2:1,3:-1,4:1,5:-1}

#print(action_of_dehn_twist_on_edge(G, number_of_edges, hom, gluing, edge_orientations, [(5,0,1), (4,0,-1), (1,0,1), (0,0,-1)], (G.identity(),4), 1))
#print(action_of_dehn_twist_on_edge(G, number_of_edges, hom, gluing, edge_orientations, [(5,0,1), (4,0,-1), (1,0,1), (0,0,-1)], (G.identity(),5), 1))

#h = HomologyGroup(G, number_of_edges, gluing, hom, edge_orientations)

#T = action_of_dehn_twist_on_homology(G, number_of_edges, hom, gluing, [(0,0,1), (1,0,-1), (2,0,1), (3,0,-1), (4,0,1), (5,0,-1), (3,1,1), (2,1,-1)], h, edge_orientations, 1)