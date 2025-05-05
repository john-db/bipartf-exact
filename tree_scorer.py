import numpy as np
from decimal import Decimal
import itertools

def comb(max, k):
    return [comb for comb in itertools.combinations(list(range(max)), k)]

def calc_prob_num(P, subtrees, order, clade, column): 
    
    cladeidx = None
    for i in range(len(subtrees)):
        if np.all(subtrees[i] == clade):
            cladeidx = i
    if cladeidx == None:
        return Decimal(0)
    else:
        subtrees_restrict = subtrees[:cladeidx] + subtrees[cladeidx+1:]
        P_restrict = np.delete(P, column, axis=1)
        col = P[:, column]
        colscore = 2 ** Decimal(np.dot(np.log2(col), clade) + np.dot(np.log2(1 - col), 1 - clade))
        s1 = calc_prob(P_restrict, subtrees, order)
        if sum(clade) == 1 or sum(clade) == len(subtrees[0]):
            return colscore * s1
        else:
            if order - 1 >= 0:
                s2 = calc_prob(P_restrict, subtrees_restrict, order - 1)
            else:
                s2 = calc_prob(P_restrict, subtrees_restrict, 0)
            return colscore * (s1 + s2)
            
#python exact_partf.py -i ./matrices/5x5_input.tsv -c "cell_2,cell_4,cell_5" -m "mut_1" -fp 0.5 -fn 0.5 -e True
def calc_prob(P, subtrees, order):
    n = len(subtrees[0])

    nontrivial = list(subtrees[n + 1:-1])

    logP = np.log2(P)
    logPC = np.log2(1 - P)
    p0 = 2 ** Decimal(prob_mat_mul_calc(logP, logPC, subtrees))
    correction = Decimal(0)
    for size in range(1, order + 1):
        temp = Decimal(0)
        for c in comb(len(nontrivial), size):
            restricted_subtrees = nontrivial[0:c[0]]
            for i in range(len(c)):
                if i < len(c) - 1:
                    restricted_subtrees += nontrivial[c[i] + 1 : c[i + 1]]
            restricted_subtrees += nontrivial[c[-1] + 1:]
            restricted_subtrees += subtrees[0:n+1]
            restricted_subtrees += [subtrees[-1]]
            temp += 2 ** Decimal(prob_mat_mul_calc(logP, logPC, np.array(restricted_subtrees)))
        if temp != 1:
            correction += (Decimal(-1) ** (size)) * temp
        else:
            print("???")
    return p0 + correction
            


def pf_cond_on_one_tree(P, subtrees, cond_c, cond_m):
    r"""
    Prob_{A\sim P}[\subtree(c, R, A)\cap A\in G| A\in T] in O(n^2).

    :param P:
    :param subtrees: cell lineage tree  n x (2n+1)
    :param cond_c: set of cells
    :param cond_m: one mutation
    :return: conditioned on the given tree what is probability of the given partition
        based on P numerator, denominator are returned separately here
    """

    denominator = Decimal(0)
    numerator = Decimal(0)
    col = P[:, cond_m]
    for v in subtrees:
        prob = Decimal(np.prod(col * v + (1 - col) * (1 - v)))
        denominator += prob
        if np.array_equal(v, cond_c):
            numerator = prob
    return numerator, denominator

def prob_mat_mul_calc(logP, logPC, subtrees):
    st = np.array(subtrees)

    res = np.exp2(np.matmul(st, logP) + np.matmul(1 - st, logPC))
    res = np.matmul(np.ones(res.shape[0]), res)
    res = np.matmul(np.ones(res.shape[0]), np.log2(res))

    return res

# def prob_mat_vec_calc(P, subtrees):
#     logP = np.log2(P)
#     logPC = np.log2(1 - P)
#     st = np.array(subtrees)
#     # res = np.exp2(np.matmul(st, lp))

#     return_value = Decimal(1)
#     for j in range(P.shape[1]):
#         # temp = Decimal(0)
#         colP = logP[:, j]
#         colPC = logPC[:, j]
#         res = np.exp2((np.matmul(st, colP) + np.matmul(1 - st, colPC)))
#         temp = Decimal(np.matmul(np.ones(res.shape[0]), res))
#         return_value *= temp
#     return return_value

# def cell_lineage_tree_prob(P, subtrees):
#     r"""
#     Calculate Prob_{A\sim P}[A\in T] in O(n m^2).

#     :param P:
#     :param subtrees: cell lineage tree
#     :return: Probability of this tree in the original distribution( based on P)
#     """

#     return_value = Decimal(1)
#     for j in range(P.shape[1]):
#         temp = Decimal(0)
#         col = P[:, j]
#         for v in subtrees:
#             temp += Decimal(np.prod(col * v + (1 - col) * (1 - v)))
#         return_value *= temp
#     return return_value
