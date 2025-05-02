import pandas as pd
from decimal import Decimal
import numpy as np
import argparse
#
from tree_scorer import calc_prob, pf_cond_on_one_tree, calc_prob_num
from utils import all_trees, newick_to_subtrees, binary_trees



def bipartf(P, cells, mutation, df, order, all_clts):
    my_cell = []
    for i in range(len(df.index)):
        if df.index[i] in cells:
            my_cell += [i]
    cond_c = np.zeros(P.shape[0], dtype=np.int8)
    cond_c[my_cell] = 1

    cond_m = np.where(df.columns == mutation)[0][0]

    n = len(df.index)

    numerator = Decimal(0)
    denominator = Decimal(0)

    treelist = None
    if all_clts == True:
        treelist = all_trees(n)
    else:
        treelist = binary_trees(n)

    if order == -1:
        for tree in treelist:
            sts = newick_to_subtrees(tree)
            nontrivial = sts[n + 1 : -1]
            
            order = len(nontrivial)
            p = calc_prob(P, sts, order)

            num = calc_prob_num(P, sts, order, cond_c, cond_m)
            numerator += Decimal(num)
            denominator += Decimal(p)
    elif order == 0:
        for tree in treelist:
            sts = newick_to_subtrees(tree)
            nontrivial = sts[n + 1 : -1]
            
            order = 0
            p = calc_prob(P, sts, order)

            cond = pf_cond_on_one_tree(P, sts, cond_c, cond_m)
            num = Decimal(cond[0] / cond[1]) * Decimal(p)
            numerator += Decimal(num)
            denominator += Decimal(p)
    else:
        raise Exception("don't use, not tested yet")
        for tree in treelist:
            sts = newick_to_subtrees(tree)
            nontrivial = sts[n + 1 : -1]
            
            order = min(len(nontrivial), order)
            p = calc_prob(P, sts, order)

            num = calc_prob_num(P, sts, order, cond_c, cond_m, len(nontrivial))
            numerator += Decimal(num)
            denominator += Decimal(p)
    return numerator, denominator, numerator / denominator

def exact_partf(P, cells, mutation, df):

    my_cell = []
    for i in range(len(df.index)):
        if df.index[i] in cells:
            my_cell += [i]
    cond_c = np.zeros(P.shape[0], dtype=np.int8)
    cond_c[my_cell] = 1

    cond_m = np.where(df.columns == mutation)[0][0]

    n = len(df.index)

    numerator = Decimal(0)
    denominator = Decimal(0)
    for tree in all_trees(n):
        sts = newick_to_subtrees(tree)
        nontrivial = sts[n + 1 : -1]
        
        order = len(nontrivial)
        p = calc_prob(P, sts, order)

        num = calc_prob_num(P, sts, order, cond_c, cond_m)
        numerator += Decimal(num)
        denominator += Decimal(p)

    return numerator, denominator, numerator / denominator

def our_partf(P, cells, mutation, df):

    my_cell = []
    for i in range(len(df.index)):
        if df.index[i] in cells:
            my_cell += [i]
    cond_c = np.zeros(P.shape[0], dtype=np.int8)
    cond_c[my_cell] = 1

    cond_m = np.where(df.columns == mutation)[0][0]

    n = len(df.index)

    numerator = Decimal(0)
    denominator = Decimal(0)
    for tree in binary_trees(n):
        sts = newick_to_subtrees(tree)
        
        order = 0
        p = calc_prob(P, sts, order)

        cond = pf_cond_on_one_tree(P, sts, cond_c, cond_m)
        num = Decimal(cond[0] / cond[1]) * Decimal(p)
        numerator += Decimal(num)
        denominator += Decimal(p)

    return numerator, denominator, numerator / denominator

def main(args):

    
    path = args.input_matrix
    cells = args.cells
    mutation = args.mutation
    alpha = args.alpha
    beta = args.beta

    df = pd.read_csv(path, sep=",", index_col=[0])

    I_mtr = df.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5

    pstr = None

    ret = bipartf(P, cells, mutation, df, order=args.order, all_clts=bool(args.all_trees))
    pstr = [str(ret[0]), str(ret[1]), str(float(ret[2]))]
    # print(args.all_trees)
    # if args.all_trees == 1:
    #     pstr = [str(exact_partf(P, cells, mutation, df)[0]), str(exact_partf(P, cells, mutation, df)[1]), str(float(exact_partf(P, cells, mutation, df)[2]))]
    # else:
    #     print("hi")
    #     pstr = [str(our_partf(P, cells, mutation, df)[0]), str(our_partf(P, cells, mutation, df)[1]), str(float(our_partf(P, cells, mutation, df)[2]))]

    output = list(map(lambda x: str(x), [args.input_matrix, args.cells, args.mutation, args.alpha, args.beta]))
    if args.all_trees == 0:
        output += ["binary"]
    else:
        output += ["all"]
    output += [str(args.order)]
    output += pstr
    print("\t".join(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run.py')

    parser.add_argument("-i", "--input_matrix", type=str,                                                        
                        help="Path to input genotype matrix where rows correspond to cells/sublines and columns correspond to mutations. See repo examples for formatting.", required=True)
    parser.add_argument("-c", "--cells", type=str,                                                        
                        help="List of cells (comma separated)", required=True)
    parser.add_argument("-m", "--mutation", type=str,                                                        
                        help="Name of the mutation (column) in the matrix", required=True)
    parser.add_argument("-fp", "--alpha", type=float,                                                        
                        help="False-positive rate (alpha in the paper)", required=True)
    parser.add_argument("-fn", "--beta", type=float,                                                        
                        help="False-negative rate (beta in the paper)", required=True)
    parser.add_argument("-a", "--all_trees", type=int, choices=[0,1],                                                    
                        help="1 for all trees, 0 for binary", required=True)
    parser.add_argument("-o", "--order", type=int,                                                    
                        help="Order of the approximation (0 for just X consistent with T, -1 for X yields T, i <= n of cells - 2^th order approximation)", required=True)
    

    main(parser.parse_args())