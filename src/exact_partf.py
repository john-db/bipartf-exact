import pandas as pd
from decimal import Decimal
import numpy as np
import argparse
import time
from tree_scorer import calc_prob, calc_prob_num

def bipartf(P, cells, mutation, df, trees):
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
    with open(trees) as f:
        for tree in f:
            tree = tree.strip()
            sts = np.zeros(shape=(len(tree) // n, n), dtype=np.int8)
            for i in range(len(tree)):
                sts[np.unravel_index(i, sts.shape)] = int(tree[i])
            nontrivial = sts[n + 1 : -1, :]
            
            order = len(nontrivial)
            p = calc_prob(P, sts, order)

            num = calc_prob_num(P, sts, order, cond_c, cond_m)
            numerator += Decimal(num)
            denominator += Decimal(p)

        

    
    return numerator, denominator, numerator / denominator

def main(args):

    
    path = args.input_matrix
    cells = args.cells
    mutation = args.mutation
    alpha = args.alpha
    beta = args.beta

    df = pd.read_csv(path, sep="\t", index_col=[0])

    I_mtr = df.values
    t1 = I_mtr * (1 - beta) / (alpha + 1 - beta)
    t2 = (1 - I_mtr) * beta / (beta + 1 - alpha)
    P = t1 + t2
    P[I_mtr == 3] = 0.5

    pstr = None

    start = time.time()
    ret = bipartf(P, cells, mutation, df, trees=args.trees)
    pstr = [str(ret[0]), str(ret[1]), str(float(ret[2]))]

    end = time.time()
    if args.cells == None or args.cells == "":
        cells_str = "Empty"
    else:
        cells_str = args.cells
    output = list(map(lambda x: str(x), [args.input_matrix, round(end-start, 5), cells_str, args.mutation, args.alpha, args.beta]))
    output += [args.trees]
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
    parser.add_argument("-t", "--trees", type=str,                                                    
                        help="Path to file containing trees", required=True)
    
    

    main(parser.parse_args())
