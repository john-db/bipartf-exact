import pandas as pd

folder = ""
folder_trees = ""

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
                        
                        path_trees = folder_trees  + str(n) + "_leaf_clade_matrices.lin"

                        path_noise = folder + sim + ".SC.after_noise"
                        df = pd.read_csv(path_noise, sep="\t", index_col=[0])
                        for j in min(m, 100):
                            mut = "mut" + str(j)
                            col = df[mut]
                            clade = [c for c in col.keys() if col[c] == 1]
                            clade = ",".join(clade)
                            print(" ".join([sim, clade, mut]))

                            cmd = "sbatch"
                            name = 'bpf-ex.' + '.'.join(params + [mut])
                            cmd += ' --job-name="' + name + '"'
                            cmd += ' --output="' + name + '.%j.out"'
                            cmd += ' --error="' + name + '.%j.err"'
                            cmd += ' --export=CLADE="' + clade + '"'
                            cmd += ',MUT="' + mut + '"'
                            cmd += ',ALPHA="' + str(fp) + '"'
                            cmd += ',BETA="' + str(fn) + '"'
                            cmd += ',INPUT="' + path_noise + '"'
                            cmd += ',TREES="' + path_trees + '"'
                            cmd += ',ORDER="' + "-1" + '"'
                            cmd += ',ALLTREES="' + "1" + '"'
                            cmd += ' bpf-exact-long-sims.sbatch'
                            print(cmd)