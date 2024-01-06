import math
import numpy as np

def run_tests(force_cosheaf, constant_cosheaf, position_cosheaf):
    F = force_cosheaf
    J = constant_cosheaf
    P = position_cosheaf
    action = F.group_action

    print("Checking commutativity of cosheaf representation and boundary matrices\n")
    for generator_index in action.generator_indices():
        print("Group generator index:", generator_index)
        for dim in [1,2]:
            print("Force Cosheaf,    dim", dim, "error:", commutative_rep_bdd_test(F, dim, generator_index))
        for dim in [1,2]:
            print("Constant Cosheaf, dim", dim, "error:", commutative_rep_bdd_test(J, dim, generator_index))
        for dim in [1,2]:
            print("Position Cosheaf, dim", dim, "error:", commutative_rep_bdd_test(P, dim, generator_index))

    print("\n---------------------------------------\n")
    print("Check commutativity of cosheaf representation and cosheaf map matrices\n")
    for generator_index in action.generator_indices():
        print("Group generator index:", generator_index)
        for dim in [0,1,2]:
            print("Force Cosheaf,    dim", dim, "error:", commutative_rep_map_test(F, J, dim, generator_index))
        for dim in [0,1,2]:
            print("Constant Cosheaf, dim", dim, "error:", commutative_rep_map_test(J, P, dim, generator_index))

    print("\n---------------------------------------\n")
    print("Check commutativity of boundary and cosheaf map matrices\n")
    for dim in [1,2]:
        print("Force Cosheaf,    dim", dim, "error:", commutative_bdd_map_test(F, J, dim))
    for dim in [1,2]:
        print("Constant Cosheaf, dim", dim, "error:", commutative_bdd_map_test(J, P, dim))



def commutative_rep_bdd_test(Cosheaf, dim, index):
    #Check error in commutativity of representation matrices and boundary map of cosheaves
    A = np.dot(Cosheaf.representation(dim-1, index), Cosheaf.bdd(dim))
    B = np.dot(Cosheaf.bdd(dim), Cosheaf.representation(dim, index))
    c = np.absolute(A-B).round(3)
    if c.shape[0] > 0 and c.shape[1] > 0: #If the boundary matrix is non-trivial
        return c.max()
    else:
        return 0

def commutative_rep_map_test(Cosheaf_domain, Cosheaf_target, dim, index):
    #Check error in commutativity of representation matrices and cosheaf maps
    A = np.dot(Cosheaf_domain.to_next_cosheaf_map(dim), Cosheaf_domain.representation(dim, index))
    B = np.dot(Cosheaf_target.representation(dim, index), Cosheaf_domain.to_next_cosheaf_map(dim))
    c = np.absolute(A-B).round(3)
    if c.shape[0] > 0 and c.shape[1] > 0: #If the cosheaf map is non-trivial
        return c.max()
    else:
        return 0

def commutative_bdd_map_test(Cosheaf_domain, Cosheaf_target, dim):
    #Check error in commutativity of boundary maps and cosheaf maps
    A = np.dot(Cosheaf_domain.to_next_cosheaf_map(dim-1), Cosheaf_domain.bdd(dim))
    B = np.dot(Cosheaf_target.bdd(dim), Cosheaf_domain.to_next_cosheaf_map(dim))
    c = np.absolute(A-B).round(3)
    if c.shape[0] > 0 and c.shape[1] > 0: #If the cosheaf map and boundary map is non-trivial
        return c.max()
    else:
        return 0