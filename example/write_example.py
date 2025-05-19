import numpy as np
lines = []

clades = [[1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1],
          [1,1,1,0,0], [0,0,0,1,1], [1,1,0,0,0], [1,1,1,1,1]]
matrix = []
n_muts = 1001
for i in range(1,n_muts):
    choice = int(np.random.rand() * len(clades))
    mut = clades[choice].copy()
    for j in range(len(mut)):
        if np.random.rand() < 0.25:
            mut[j] = 3
        else:
            if mut[j] == 0:
                if np.random.rand() < 0.01:
                    mut[j] = 1
            else:
                if np.random.rand() < 0.1:
                    mut[j] = 0
    matrix += [mut]

header = ["cell_id_x_mut_id"] + ["mut_" + str(i) for i in range(1,n_muts)]
matrix = np.array(matrix, dtype=int).T
print("\t".join(header))
for i in range(matrix.shape[0]):
    print("cell_" + str(i + 1), end="\t")
    print("\t".join(map(str, list(matrix[i]))))
