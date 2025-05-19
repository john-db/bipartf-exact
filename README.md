# bipartf-exact

```
TEMP_NAME=$(mktemp); python src/gen_trees.py -n 5 > "$TEMP_NAME" && \
python src/exact_partf.py -i ./example/example-5x5.tsv -c "cell_1,cell_2" -m "mut_1" -fp 0.01 -fn 0.1 -t "$TEMP_NAME"
```
