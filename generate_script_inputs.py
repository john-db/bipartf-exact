import pandas as pd

folder = ""

n_ls = [5,6,7,8]
m_ls = [8,10,12,15,20,50,100,200,500,1000]
fp_ls = [0.001]
fn_ls = [0.1,0.2,0.4]
na_ls = ["0","0.25","0.50"]
repl_ls = range(1,11)
for repl in repl_ls:
    for n in n_ls:
        for m in m_ls:
            for fp in fp_ls:
                for fn in fn_ls:
                    for na in na_ls:
                        params = [repl,n,m,fp,fn,na]
                        params = map(str, params)
                        sim = "simNo_{}-n_{}-m_{}-fp_{}-fn_{}-na_{}".format(*params)
                        sim += ".SC.after_noise"
                        
                        path = folder + sim
                        df = pd.read_csv(path, sep="\t", index_col=[0])
                        for j in min(m, 100):
                            mutation = "mut" + str(j)
                            col = df[mutation]
                            clade = [c for c in col.keys() if col[c] == 1]
                            print([sim, clade, mutation])