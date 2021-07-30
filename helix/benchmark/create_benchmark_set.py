'''
Scan through the dataframe consisting of all the helices in PDB
interfaces. For each interface, see how many helices are part of the
other side of that interface. Then see how may of those have at least x
# of residues as part of the interface.

Usage:
    create_benchmark_set.py <dataframe_path> [options]

Options:
    --out=PATH, -o  Where to save final dataframe pickle  
    [default: interface_finder/benchmark_candidates.pkl]
    --interface_residues=NUM, -n  How many residues must be in the interface
    for a helix to be considered "flat"  [default: 5]
    --percent_interface=FLOAT, -p  Percentage of helix residues that
    must participate in interface for it to be considered "flat"
'''
import docopt
import sys
import pandas as pd
import pickle5 as pickle


def main():
    '''
    Here's what we gonna do.
    Scan through the dataframe.
    For each monomer (interface), see how many  helices.
    How many of those have at least, idk, 10 residues as part of the
    interface?
    Ok, so save those.
    '''
    args = docopt.docopt(__doc__)
    with open(args['<dataframe_path>'], 'rb') as f:
        df = pickle.load(f)
    # df = pd.read_pickle(args['<dataframe_path>'])
    df = df[df['interacting_chain'].apply(lambda x: x is not None)]
    df['target'] = df['interface'].apply(lambda x: x.split('_')[0])
    print
    df = df[df['interacting_length'].apply(lambda x: int(x) >=
        args['--interface_residues'])]
    df = df[df['chain'] != df['target']]
    df = df[df['interacting_chain'] == df['target']]
    groups = df.groupby(by=['name', 'target', 'chain'])
    final_df = pd.DataFrame()
    total_groups = groups.ngroups
    i = 0
    for name, group in groups:
        i += 1
        print('Processing group {} of {}'.format(i, total_groups),
                end='\r')
        group['n_helices'] = group.shape[0]
        final_df = pd.concat([final_df, group], ignore_index=True)

    final_df = final_df[final_df['n_helices'] > 2]

    final_df.to_pickle(args['--out'])


if __name__=='__main__':
    main()
