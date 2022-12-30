

# This file was *autogenerated* from the file index.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_3 = Integer(3)
import itertools

def entry_max(M):
    """ Given a matrix M, returns the max absolute value of its entries. """
    return max([abs(M[i][j]) for i in range(_sage_const_0 ,M.nrows()) for j in range(_sage_const_0 ,M.ncols())])


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
            prod_mod_2 = matrix(GF(_sage_const_2 ), prod)
            prod_mod_2.set_immutable()
            rep = transversal[prod_mod_2]
            new_gen = prod*(rep.inverse())
            subgroup_gens.append(new_gen)
    return subgroup_gens


def word_in_sanov_generators(mat):
    """ Given a matrix mat in the sanov subgroup, returns list of integers representing word in free generators. """

    # We decompose mat recursively, keeping track of the word we've accumulated
    def sanov_recursion(M, left_word, right_word):
        A = matrix([[_sage_const_1 ,_sage_const_2 ],[_sage_const_0 ,_sage_const_1 ]])
        B = matrix([[_sage_const_1 ,_sage_const_0 ],[_sage_const_2 ,_sage_const_1 ]])
        labels = {_sage_const_1  : A, -_sage_const_1  : A.inverse(), _sage_const_2  : B, -_sage_const_2  : B.inverse()}
        entry_max_M = entry_max(M)

        # Base case: M is the identity
        if entry_max_M == _sage_const_1 :
            right_word.reverse()
            word = [-_sage_const_1 *l for l in left_word] + [-_sage_const_1 *r for r in right_word]
            
            # We can check that we have successfully decomposed mat
            assert prod([labels[i] for i in word])*identity_matrix(_sage_const_2 ) == mat
            
            return word

        else:
            for i in [_sage_const_1 ,-_sage_const_1 ,_sage_const_2 ,-_sage_const_2 ]:
                C = labels[i]
                if entry_max(C*M) < entry_max_M:
                    left_word.append(i)
                    result = sanov_recursion(C*M, left_word, right_word)
                    return result
                if entry_max(M*C) < entry_max_M:
                    right_word.append(i)
                    result = sanov_recursion(M*C, left_word, right_word)
                    return result

    return sanov_recursion(mat, [], [])


def transversal(gens):
    """ Given a generating set of a subgroup H of SL(2,Z), returns a transversal of H cap Gamma(2) in H. """
    
    image_mod_2 = SL(_sage_const_2 ,_sage_const_2 ).subgroup([matrix(GF(_sage_const_2 ), gen) for gen in gens])
    immutable_image = [matrix(mat) for mat in image_mod_2]
    for mat in immutable_image:
        mat.set_immutable()

    # To build the transversal, we iterate over longer and longer words in the generators until we fill up the whole image
    transversal = {}
    reps_found = set()
    word_length = _sage_const_1 
    while len(reps_found) < image_mod_2.order():
        copy_of_gens = [gens for i in range(_sage_const_0 , word_length)]
        words = itertools.product(*copy_of_gens)
        for word in words:
            group_element = product(word)
            group_element_mod_2 = matrix(GF(_sage_const_2 ), group_element)
            group_element_mod_2.set_immutable()

            if group_element_mod_2 not in reps_found:
                reps_found.add(group_element_mod_2)
                transversal[group_element_mod_2] = group_element
        word_length = word_length + _sage_const_1 

    return transversal


def is_finite_index(gens):
    """Given a list of generators of a subgroup Gamma of SL(2,Z), determines whether the index of Gamma is finite."""

    trans = transversal(gens)
    intersection_gens = schreier(gens, trans)

    # If necessary, can replace each gen with -1*gen to assume we're in the Sanov subgroup
    # The Sanov subgroup has diagonal matrices congruent to 1 mod 4
    for gen in intersection_gens:
        if gen[_sage_const_0 ][_sage_const_0 ] % _sage_const_4  == _sage_const_3 :
            gen = -_sage_const_1 *gen
    
    F = FreeGroup('a,b')
    free_subgroup_gens = []
    for gen in intersection_gens:
        free_subgroup_gens.append(F(word_in_sanov_generators(gen)))

    # We use GAP to check the index of the subroup
    gap.eval("F:=FreeGroup(2)")
    gap.eval("a:=F.1")
    gap.eval("b:=F.2")
    gap.eval("H:=Subgroup(F,[ %s ])" % ", ".join([str(w) for w in free_subgroup_gens]))
    val = gap.eval("Index(F,H)")
    return (val != "infinity")

