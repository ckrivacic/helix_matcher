import pandas as pd
import sys, os


def overlap(row1, row2):
    '''Check if there is overlap between 2 helices'''
    resis1 = row1['pdb_resis']
    resis2 = row2['pdb_resis']
    return any(x in max(resis1,resis2,key=len) for x in min(resis2,resis1,key=len))


def merge(row1, row2):
    '''Merge 2 rows and preserve interface score and other information'''
    keys = ['rosetta_resis', 'interface_resis', 'pdb_resis',
            'interface_scores']
    row1_indices = [i for i, resi in enumerate(row1['pdb_resis']) if
            resi not in row2['pdb_resis']]
    interface_indices = [i for i, resi in
            enumerate(row1['interface_resis']) if resi not in
            row2['interface_resis']]
    for key in keys:
        if key=='interface_resis' or key=='interface_scores':
            for idx in interface_indices:
                row2[key].append(row1[key][idx])
        else:
            for idx in row1_indices:
                row2[key].append(row1[key][idx])
    return row2


def merge_all(group, verbose=False):
    '''Recursively combine overlapping helices'''
    for idx1, row1 in group.iterrows():
        for idx2, row2 in group.iterrows():
            if idx1 != idx2:
                print('Checking {} vs. {}'.format(idx1, idx2))
            if overlap(row1, row2) and idx1 != idx2:
                if verbose:
                    print('Overlap found.')
                    print('Row1: {}'.format(row1['pdb_resis']))
                    print('Row2: {}'.format(row2['pdb_resis']))
                row2 = merge(row1, row2)
                group.loc[idx2] = row2
                if verbose:
                    print('Merged; new row2:')
                    print(group.loc[idx2]['pdb_resis'])
                group = group.drop(idx1)
                if verbose:
                    print(group)
                group = merge_all(group)
                return group
    return group



def main():
    try:
        df = pd.read_pickle(sys.argv[1])
    except:
        import pickle5 as pickle
        with open(sys.argv[1], 'rb') as f:
            df = pickle.load(f)
    groups = df.groupby(['name', 'target', 'chain'])
    final_df = pd.DataFrame()
    for name, group in groups:
        group = merge_all(group)
        group['n_helices'] = group.shape[0]
        final_df = pd.concat([final_df, group])
    
    final_df.to_pickle('final_consolidated.pkl')

if __name__=='__main__':
    main()
