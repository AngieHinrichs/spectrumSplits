import os
import sys
import re
import argparse
import bte
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Process bootstrap spectra output files.")
    parser.add_argument('--spectrum_file', type=str, default="spectra_output.tsv", help="Spectrum file computed from the tree.")
    parser.add_argument("--input_tree", type=str, default="public-latest.all.masked.pb.gz", help="Input tree file (protobuf format)")
    parser.add_argument("--annotate_nodes_output_file", type=str, default="annotated_nodes_output.tsv", help="File to save annotated node data.")
    parser.add_argument("--metadata_output", type=str, default="metadata_output.tsv", help="File to save post-processed metadata.")
    parser.add_argument("--fold", action="store_true", help="Also write the six strand-folded classes (CT|GA and so on) produced by spectrumSplits.py --fold, after the twelve unfolded columns. If the spectrum file has no folded columns, each is computed as the sum of its pair.")
    return parser.parse_args()

### The twelve unfolded classes, and the six strand-folded pairs in the column
### order spectrumSplits.py --fold writes them.  Each folded label names its
### two members, e.g. "CT|GA" is CT plus GA.
MUTATION_TYPES = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]
FOLDED_TYPES = ["CA|GT","CG|GC","CT|GA","TA|AT","TC|AG","TG|AC"]

def folded_value( row, label ) :
    """The folded proportion for one split: read from the file if spectrumSplits
    wrote it, otherwise the sum of the pair (which is what it writes)."""
    if label in row :
        return row[ label ]
    return sum( float( row.get( member, 0 ) ) for member in label.split( "|" ) )

def import_tsv_to_dict(tsv_file):
    data = {}
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            # because bootstraps won't have representative nodes, skip that.
            row_dict = {header[i]: values[i] for i in range(len(header)-1)}
            data[row_dict[header[0]]] = row_dict
    return data

def get_spectrum_roots(tree, subset_nodes_set):
    annotations = {}

    def traverse_and_annotate(node, current_ancestor):
        # If the node is in the subset, update the current ancestor
        if node.id in subset_nodes_set:
            current_ancestor = node.id
        # Annotate the node and its tips with the current ancestor
        annotations[node.id] = current_ancestor

        # Recursively traverse children
        for child in node.children:
            traverse_and_annotate(child, current_ancestor)

    # Start the traversal from the root
    traverse_and_annotate(tree.root, None)

    return annotations

if __name__ == "__main__":
    # Take input data
    args = parse_args()
    spectra_data = import_tsv_to_dict(args.spectrum_file)
    splits = set(spectra_data.keys())
    tree = bte.MATree(args.input_tree)

    # Get the parent in the splits
    spectrum_roots = get_spectrum_roots(tree, splits)

    columns = MUTATION_TYPES + ( FOLDED_TYPES if args.fold else [] )
    if args.fold :
        first_row = next( iter( spectra_data.values() ) )
        missing = [ label for label in FOLDED_TYPES if label not in first_row ]
        if missing :
            print( f"Note: {args.spectrum_file} has no folded columns ({', '.join(missing)}); computing them as the sum of each pair.", file=sys.stderr )

    # Create a DataFrame to store the results
    result_data = { "NodeID": [], "SepctrumRoot": [] }
    for column in columns :
        result_data[ column ] = []

    # Populate the data
    for node, s in spectrum_roots.items():
        result_data["NodeID"].append(node)
        result_data["SepctrumRoot"].append(s)
        for mutation_type in MUTATION_TYPES :
            result_data[ mutation_type ].append( spectra_data[s][ mutation_type ] )
        if args.fold :
            for label in FOLDED_TYPES :
                result_data[ label ].append( folded_value( spectra_data[s], label ) )

    # Convert result data into a pandas DataFrame
    df = pd.DataFrame(result_data)

    # Save the annotated node data
    df.to_csv(args.annotate_nodes_output_file, sep='\t', index=False)

    # Post-process to rename 'NodeID' to 'strain' and save the output
    df.rename(columns={"NodeID": 'strain'}, inplace=True)  # rename the 'NodeID' column to 'strain'
    df.to_csv(args.metadata_output, sep='\t', index=False)  # save it as a .tsv extension
