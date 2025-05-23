# bipartf-exact

All commands are given supposing that your shell is running in the top-level directory of this repo (i.e., the one containing the ```example/``` and ```src/``` directories).

To run the exact bipartion program we first generate a text file containing binary strings which represent trees. (the strings represent a linearized version of a matrix where each row corresponds to a (trivial or nontrivial) clade of the tree, where the ith entry of the row equals 1 if the ith cell is in the clade, and 0 otherwise.

To generate this file, please use gen_trees.py as such:
```python src/gen_trees.py -n <num_leaves> > <tree_filename>```
where ```<num_leaves>``` is the number of leaves that the tree should contain (i.e., the number of cells or rows in your input matrix).

Then to run the exact bipartition program: ```python src/exact_partf.py -i <path_to_input_matrix> -c <clade> -m <mutation> -fp <false_positive_rate> -fn <false_negative_rate> -t <tree_filename>```
where:

```<path_to_input_matrix>``` is the path to your input genotype matrix (see examples in example/ directory)

```<clade>``` is the clade of interest given as a comma-separated list of labels of the rows of the input matrix

```<mutation>``` label of the mutation of interest from input matrix

```<false_positive_rate>``` and ```<false_negative_rate>``` are the probablities of false-positives and false-negatives (expressed as real numbers between 0 and 1)

```<tree_filename>``` is the path to the file containing the trees (generated by gen_trees.py)

Here are some example commands that may be used to run the exact bipartition function on the data in ```example/```
5x5 example:
```
TEMP_NAME=$(mktemp); python src/gen_trees.py -n 5 > "$TEMP_NAME" && python src/exact_partf.py -i ./example/example-5x5.tsv -c "cell_1,cell_2" -m "mut_1" -fp 0.01 -fn 0.1 -t "$TEMP_NAME"
```

5x1000 example:
```
TEMP_NAME=$(mktemp); python src/gen_trees.py -n 5 > "$TEMP_NAME" && python src/exact_partf.py -i ./example/example-5x1000.tsv -c "cell_1,cell_2" -m "mut_1" -fp 0.01 -fn 0.1 -t "$TEMP_NAME"
```
(note that these commands are generating trees and storing them to a temporary file, and then running the exact bipartition function)
