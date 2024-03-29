import itertools

def entry_max(M):
    """ Given a matrix M, returns the max absolute value of its entries. """
    return max([abs(M[i][j]) for i in range(0,M.nrows()) for j in range(0,M.ncols())])


def schreier(gens, transversal):
    """ 
    Given a generating set of a subgroup H of SL(2,Z), returns a generating set for intersection with Gamma(2) 

    Parameters
    ----------
    gens: list of matrices 
        generators of subgroup
    transversal: dict
        given a matrix N in SL(2,2), transversal[N] is a representative of N in H
    """
    subgroup_gens = []
    for S in gens:
        for R in list(transversal.values()):
            prod = R*S
            prod_mod_2 = matrix(GF(2), prod)
            prod_mod_2.set_immutable()
            rep = transversal[prod_mod_2]
            new_gen = prod*(rep.inverse())
            subgroup_gens.append(new_gen)
    return subgroup_gens


def word_in_sanov_generators(mat):
    """ Given a matrix mat in the sanov subgroup, returns list of integers representing word in free generators. """

    # We decompose mat recursively, keeping track of the word we've accumulated
    def sanov_recursion(M, left_word, right_word):
        A = matrix([[1,2],[0,1]])
        B = matrix([[1,0],[2,1]])
        labels = {1 : A, -1 : A.inverse(), 2 : B, -2 : B.inverse()}
        entry_max_M = entry_max(M)

        # Base case: M is the identity
        if entry_max_M == 1:
            right_word.reverse()
            word = [-1*l for l in left_word] + [-1*r for r in right_word]
            
            # We can check that we have successfully decomposed mat
            assert prod([labels[i] for i in word])*identity_matrix(2) == mat
            
            return word

        else:
            for i in [1,-1,2,-2]:
                C = labels[i]
                if entry_max(C*M) < entry_max_M:
                    left_word.append(i)
                    result = sanov_recursion(C*M, left_word, right_word)
                    return result
                elif entry_max(M*C) < entry_max_M:
                    right_word.append(i)
                    result = sanov_recursion(M*C, left_word, right_word)
                    return result

    return sanov_recursion(mat, [], [])


def transversal(gens):
    """ Given a generating set of a subgroup H of SL(2,Z), returns a transversal of H cap Gamma(2) in H. """
    
    image_mod_2 = SL(2,2).subgroup([matrix(GF(2), gen) for gen in gens])
    immutable_image = [matrix(mat) for mat in image_mod_2]
    for mat in immutable_image:
        mat.set_immutable()

    # To build the transversal, we iterate over longer and longer words in the generators until we fill up the whole image
    transversal = {}
    reps_found = set()
    word_length = 1
    while len(reps_found) < image_mod_2.order():
        copy_of_gens = [gens for i in range(0, word_length)]
        words = itertools.product(*copy_of_gens)
        for word in words:
            group_element = product(word)
            group_element_mod_2 = matrix(GF(2), group_element)
            group_element_mod_2.set_immutable()

            if group_element_mod_2 not in reps_found:
                reps_found.add(group_element_mod_2)
                transversal[group_element_mod_2] = group_element
        word_length = word_length + 1

    return transversal


def is_finite_index(gens):
    """Given a list of generators of a subgroup Gamma of SL(2,Z), determines whether the index of Gamma is finite."""

    trans = transversal(gens)
    intersection_gens_unadjusted = schreier(gens, trans)
    
    # If necessary, can replace each gen with -1*gen to assume we're in the Sanov subgroup
    # The Sanov subgroup has diagonal matrices congruent to 1 mod 4
    intersection_gens = [gen if gen[0][0] % 4 == 1 else -1*gen for gen in intersection_gens_unadjusted]
    
    F = FreeGroup('a,b')
    free_subgroup_gens = []
    for gen in intersection_gens:
        free_subgroup_gens.append(F(word_in_sanov_generators(gen)))

    # We use GAP to check the index of the subroup
    gap.eval("F:=FreeGroup(2)")
    gap.eval("a:=F.1")
    gap.eval("b:=F.2")
    gap.eval("H:=Subgroup(F,[ %s ])" % ", ".join([str(w) for w in free_subgroup_gens if w != F.one()]))
    val = gap.eval("Index(F,H)")
    return (val != "infinity")