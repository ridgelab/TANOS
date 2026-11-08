import sys
from collections import defaultdict

import dendropy


main_tree_path = sys.argv[1]
replicate_trees_path = sys.argv[2]
output_path = sys.argv[3]

tns = dendropy.TaxonNamespace()
main_tree = dendropy.Tree.get(
    path=main_tree_path,
    schema="newick",
    taxon_namespace=tns,
)
replicate_trees = dendropy.TreeList.get(
    path=replicate_trees_path,
    schema="newick",
    taxon_namespace=tns,
)

main_tree.encode_bipartitions()
support_counts = defaultdict(lambda: [0, 0])  # [count_present, count_evaluated]

for i, rep_tree in enumerate(replicate_trees):
    if not i % 10:
        print(f"Processing replicate {i}...", file=sys.stderr)

    rep_taxa = frozenset(leaf.taxon.label for leaf in rep_tree.leaf_node_iter())

    # Prune the main tree to the same taxon set as this replicate.
    pruned = main_tree.clone()
    pruned.prune_taxa_with_labels(
        [t.label for t in pruned.taxon_namespace if t.label not in rep_taxa]
    )
    pruned.encode_bipartitions()
    rep_tree.encode_bipartitions()

    rep_bips = {b.split_bitmask for b in rep_tree.bipartition_encoding}

    for bip in pruned.bipartition_encoding:
        key = bip.split_bitmask
        support_counts[key][1] += 1
        if key in rep_bips:
            support_counts[key][0] += 1

# Annotate main tree nodes with support values.
for node in main_tree.internal_nodes():
    if node.edge.bipartition is not None:
        key = node.edge.bipartition.split_bitmask
        if key in support_counts:
            present, evaluated = support_counts[key]
            support = round((present / evaluated) * 100, 1)
            node.label = str(support)
        else:
            node.label = "NA"

main_tree.write(
    path=output_path,
    schema="newick",
    suppress_internal_node_labels=False,
    suppress_rooting=True,
)

print(f"Tree with jackknife support written to {output_path}", file=sys.stderr)
