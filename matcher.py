'''
Create bins or match a query protein.

Usage:
    matcher.py bin <helix_dataframe> [options]
    matcher.py match <pdb_folder> [options]

options:
    --local, -l  Run locally

    --tasks=NUM, -j  Run on the cluster using SGE. Argument should be # of
    tasks in total.

    --length, -e  Bin by length

    --verbose, -v  Verbose output

    --database=PATH, -d  Database of relative helix orientations  
    [default: database/]

    --out=PATH, -o  Where to save outputs  [default: .]

    --angstroms=NUM, -a  Binning option. How fine should the distance bins
    be?  [default: 2.5]
    --degrees=NUM, -g  Binning option. How fine should the angle bins be?
    [default: 15]
'''
import docopt
from copy import deepcopy
import pandas as pd
import numpy as np
import numeric
from itertools import product
import os, psutil, sys
import pickle
import subprocess
from scan_helices import final_vector
from pyrosetta import init, pose_from_file
import networkx as nx
# import graph_tool.all as gt
import collections


def plot_vectors(vectors, color='darkgray'):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for vector in vectors:
        x = [point[0] for point in vector]
        y = [point[1] for point in vector]
        z = [point[2] for point in vector]
        ax.plot(x, y, z, color=color, linewidth=4)
    plt.show()


def bin_array(array, bins):
    '''
    Digitize a numpy array.
    TO DO: Circularize the binning of angles somehow.
    Done?
    '''
    upper = bins[-1]
    lower = bins[0]
    array = numeric.wrap_angles(array, 0, upper, lower)
    inds = np.digitize(array, bins)
    binned = tuple([bins[inds[n]-1] for n in range(array.size)])
    return binned


def relative_position(row1, row2, vectortype='vector', clash=False,
        reverse=False):
    '''
    Gives the internal relative orientation of two lines, given their
    row from the pandas dataframe created in scan_helices.
    The relative orientation of two lines should be able to be described
    with just 4 parameters, since they are 2D objects in 3D space. If we
    have lines consisting of points [a,b] and [c,d], those parameters are:
    - The distance between their centroids
    - Angle abc
    - Angle bcd
    - Dihedral abcd
    '''

    out = {}
    norm_v1 = row1[vectortype]
    norm_v2 = row2[vectortype]
    out['dist'] = numeric.euclidean_distance(row1['centroid'],
            row2['centroid'])
    out['abc'] = numeric.angle(norm_v1[0], norm_v1[1], norm_v2[0])
    out['bcd'] = numeric.angle(norm_v1[1], norm_v2[0], norm_v2[1])
    out['dih'] = numeric.dihedral(norm_v1[0], norm_v1[1], norm_v2[0],
            norm_v2[1])
    # plot_vectors([norm_v1, norm_v2], color='black')
    if clash:
        centroid_vector1 = row1['centroid_vector']
        centroid_vector2 = row2['centroid_vector']
        out['cen1a'] = numeric.angle(
                norm_v1[0], centroid_vector1[0], centroid_vector1[1]
                )
        out['cen2a'] = numeric.angle(
                norm_v2[0], centroid_vector2[0], centroid_vector2[1]
                )
        out['cen1dih'] = numeric.wrap_angle(
                numeric.dihedral(
                    norm_v1[0], norm_v2[0],
                    centroid_vector1[0], centroid_vector1[1]
                    ), 180
                )
        out['cen2dih'] = numeric.dihedral(
                norm_v2[0], norm_v1[0],
                centroid_vector2[0], centroid_vector2[1]
                )
        if reverse:
            out['cen1a'] = numeric.wrap_angle(out['cen1a'], addition=180)
            out['cen2a'] = numeric.wrap_angle(out['cen2a'], addition=180)
            out['cen1dih'] = numeric.wrap_angle(out['cen1dih'], addition=180)
            out['cen2dih'] = numeric.wrap_angle(out['cen2dih'], addition=180)

    return out


class Match(object):
    '''
    Class to construct a potential match.
    '''
    def __init__(self, name, query_db, main_db, verbose=False):
        self.verbose = verbose
        self.name = name
        self.query = query_db
        self.db = main_db.xs(name, level='name')
        # self.graph = gt.Graph(directed=False)
        self.graph = nx.Graph()
        # Track helix pairs so we don't add them to the graph more than
        # once

    def max_subgraph(self):
        '''
        Finds dense subgraphs, which represent compatible sets of helix
        pairs between the query helices and the database PDB. The
        longest such subgraph represents the best overlay of the PDB
        with the set of query helices.
        '''
        max_subgraph_len = 0
        # for f in gt.max_cliques(self.graph):
        for f in nx.find_cliques(self.graph):
            if len(f) > max_subgraph_len:
                max_subgraph_len = len(f)

        print('Max number of matches:')
        print(max_subgraph_len)
        return max_subgraph_len


    def plot_graph(self):
        import matplotlib.pyplot as plt
        import graph_tool.draw as draw
        plt.subplot(111)
        # gt.remove_parallel_edges(self.graph)
        # pos = gt.fruchterman_reingold_layout(self.graph, n_iter=1000)
        # gt.graph_draw(self.graph, pos=pos)
        plt.show()

    def find_edges(self):
        '''
        Populate the graph with nodes and edges.
        Each node consists of a pair of indices, one from the main
        database and one from the query database. This pairing
        represents the case where the helix in the first index is
        overlaid on the helix of the second index. Edges represent
        compatibility between adjacent nodes.
        '''
        print('Finding edges')
        edges = []
        self.nodes = set()
        property_map = {}
        i = 0
        for doc in self.db.iterrows():
            if doc[0] in self.query.index:
                compatible_bins = self.query.xs(doc[0])
                # compatible_bins = self.query.find({'bin': doc['bin']})
                for result in compatible_bins.iterrows():
                    idx_pair1 = (doc[1]['idx1'], result[1]['idx1'])
                    idx_pair2 = (doc[1]['idx2'], result[1]['idx2'])
                    # Track which nodes have been sampled
                    if idx_pair1 not in self.nodes:
                        self.nodes.add(idx_pair1)
                        self.graph.add_node(idx_pair1)
                        # self.nodes[idx_pair1] = i
                        # property_map[i] = idx_pair1
                        i += 1
                        # self.nodes.append(idx_pair1)
                        # self.graph.add_node(idx_pair1)
                    if idx_pair2 not in self.nodes:
                        # self.nodes[idx_pair2] = i
                        # property_map[i] = idx_pair2
                        self.nodes.add(idx_pair2)
                        self.graph.add_node(idx_pair2)
                        i += 1
                        # self.nodes.append(idx_pair2)
                        # self.graph.add_node(idx_pair2)
                    self.graph.add_edge(idx_pair1, idx_pair2)
                    # print('Edge found:')
                    # print(idx_pair1)
                    # print(idx_pair2)
                    # edges.append((self.nodes[idx_pair1],
                        # self.nodes[idx_pair2]))
                # i += 2
        # nodes = set(self.nodes)
        # self.graph.add_edge(idx_pair1, idx_pair2)
        # print(nodes)
        # if self.verbose:
            # print('All edges:')
            # print(edges)
        # self.graph.add_edge_list(edges)

        # Add properties
        # prop_dict = self.graph.new_vertex_property('object')
        # for v in self.graph.vertices():
            # prop_dict[v] = {'query_idx':property_map[v][0],
                    # 'lookup_idx':property_map[v][1]}


class HelixBin(object):
    def __init__(self, helix_db, exposed_cutoff=0.3, length_cutoff=10.8,
            angstroms=2.5, degrees=15, clash_angle=None,
            verbose=False, start=None, stop=None):
        self.verbose = verbose
        self.df = helix_db
        self.df['idx'] = self.df.index

        # Binning parameters
        self.degrees = degrees
        self.angstroms = angstroms
        self.setup_bins()
        binned_name = 'bins_{}A_{}D'.format(self.angstroms,
                self.degrees)
        self.start = start
        self.stop = stop

        # Trimming dataframe
        if length_cutoff:
            self.df = self.df[self.df['length'] > length_cutoff]
        if exposed_cutoff:
            self.df = self.df[self.df['percent_exposed'] >
                    exposed_cutoff]
        # if 'normalized_vector' not in self.df.columns:
            # self.df['normalized_vector'] = self.df.apply(lambda x:
                    # final_vector(x['direction'], 1, x['centroid']), axis=1)

        self.clash_angle = clash_angle

    def setup_bins(self):
        nrbins = int(360//self.degrees) + 1
        self.rbins = np.linspace(-180, 180, nrbins)
        tstart = -10000
        tstop = 10000
        ntbins = int((tstop - tstart) // self.angstroms) + 1
        self.tbins = np.linspace(tstart, tstop, ntbins)

    def setup_clash_bins(self):
        start = -180
        bins = [start]
        while start < 180:
            start += (180 - self.clash_angle)
            bins.append(start)
        bins = np.array(bins)
        return bins, bins + (self.clash_angle / 2)

    def bin_clashes(self, cen_vector_angles):
        bins = self.setup_clash_bins()
        for angle in cen_vector_angles:
            upper = bins[-1]
            lower = bins[0]
            array = numeric.wrap_angles(cen_vector_angles, 0, upper, lower)
            arrays = [array]

            n = upper - 180
            overshot_angles = []
            for angle in array:
                if lower < angle < (lower + n):
                    overshot_angles.append(360 + angle)
                else:
                    overshot_angles.append(angle)
            if overshot_angles != array:
                arrays.append(overshot_angles)
                    
            binned_list = []
            for arr in arrays:
                inds = np.digitize(arr, bins)
                for bin_array in bins:
                    binned_list.append(tuple([bin_array[inds[n]-1] for n in
                        range(arr.size)]))
        return binned_list

    def bin_db(self, outdir=None, bin_length=False, clash_angle=None):
        '''
        Bin dataframes.
        Outdir: Where to save binned dataframe
        Bin length: Whether to bin by helix length
        Clash: Angle for clash filter bins
        '''

        from scipy.spatial.transform import Rotation as R
        import subprocess
        import time

        # Pandas indices are hash lookups and we can have multiple of
        # them, but they cannot be added piecewise. Therefore we will
        # create partial tables, then create the indices and save the
        # dataframes. Results will be saved in chunks.
        bins = pd.DataFrame(columns=['bin', 'name', 'idx1', 'idx2'])
        total_proteins = len(set(self.df['name']))
        interval = 500

        # i tracks # of names analyzed
        i = 0
        # saveno tracks how many dataframes have been saved.
        self.saveno = 1
        unsaved_docs = []
        start_time = time.time()

        def update(bins, start_time, unsaved_docs, interval, i,
                final=False):
            '''
            Helper function to add many PDBs to the dataframe at once
            and print out an estimated time to completion
            '''
            print('{} of {} PDBs processed so far.'.format(
                i, total_proteins))
            mem_used = psutil.Process(os.getpid()).memory_info().rss
            if self.verbose:
                print('Currently using {} GB of memory'.format(
                    mem_used * 10**-9
                    ))
            df_mem = bins.memory_usage(index=True, deep=True).sum()
            if self.verbose:
                print('Dataframe is using {} GB of memory'.format(
                    df_mem * 10**-9
                    ))
            elapsed = time.time() - start_time
            rate = interval / elapsed
            remaining = (total_proteins - i) / rate / 3600
            print('Analysis of 500 pdbs took {} seconds. Est. {} h remaining'.format(
                elapsed, remaining
                ))

            if len(unsaved_docs) > 0:
                if self.verbose:
                    print('Adding to dataframe...')
                bins = bins.append(unsaved_docs, ignore_index=True)
                if self.verbose:
                    print(bins)
            else:
                if self.verbose:
                    print('Nothing to update for this batch.')


            # Save when memory footprint of dataframe gets larger than 4
            # GB. This way each sub-dataframe can be read into memory.
            if outdir:
                if df_mem * 10**-9 > 4 or final:
                    bins.set_index(['bin', 'name'], inplace=True)
                    outfile = 'bins_{}A_{}D_{:04d}.pkl'.format(self.angstroms,
                            self.degrees, self.saveno)
                    out = os.path.join(outdir, outfile)
                    print('Saving current dataframe to {}'.format(out))
                    if not os.path.exists(outdir):
                        os.makedirs(outdir, exist_ok=True)
                    bins.to_pickle(out)
                    self.saveno += 1
                    if self.verbose:
                        print('Saved.')

                    # If saved to disk, return an empty dataframe.
                    return pd.DataFrame()

            elif final:
                bins.set_index(['bin', 'name'], inplace=True)
            # Return input dataframe if we have not saved it to disk.
            return bins

        '''
        Grouping dataframe
        '''
        groups = self.df.groupby(['name'])
        names = sorted(list(groups.groups.keys()))
        if self.start:
            names = names[self.start:]
        if self.stop:
            names = names[:self.stop]
        for name in names:
        # for name, group in df.groupby(['name']):
            group = groups.groups[name]
            i += 1

            for combination in product(self.df.loc[group].T.to_dict().values(),
                    repeat=2):
                if combination[0]['idx'] != combination[1]['idx']:
                    # Store index
                    idx1 = combination[0]['idx']
                    idx2 = combination[1]['idx']

                    # Get relative orientation
                    relative =\
                            relative_position(combination[0],
                                    combination[1], clash=clash_angle is not None)
                    dist = np.array([relative['dist']])
                    angles = np.array([relative['abc'], relative['bcd'], relative['dih']])
                    if clash_angle:
                        clash_angles = np.array([relative['cen1a'], relative['cen2a'],
                            relative['cen1dih'], relative['cen2dih']])
                        clashbin1, clashbin2 = self.bin_clashes(clash_angles)
                        clasha1 = set([clashbin1[0], clashbin2[0]])
                        clasha2 = set([clashbin1[1], clashbin2[1]])
                        clashdih1 = set([clashbin1[2], clashbin2[2]])
                        clashdih2 = set([clashbin1[3], clashbin2[3]])
                        clashes = [clasha1, clasha2, clashdih1,
                                clashdih2]

                    # Bin by length
                    lengths = np.array([combination[0]['length'],
                        combination[1]['length']])
                    lbin = bin_array(lengths, self.tbins)
                    lbin2 = bin_array(lengths, self.tbins +
                            (self.angstroms/2))
 
                    # Bin by distance and angle
                    rbin = bin_array(angles, self.rbins)
                    tbin = bin_array(dist, self.tbins)
                    rbin2 = bin_array(angles, 
                            self.rbins + (self.angstroms/2))
                    tbin2 = bin_array(dist, self.tbins +
                            (self.angstroms/2))


                    # Combine bins
                    x = set([tbin[0], tbin2[0]])
                    abc = set([rbin[0], rbin2[0]])
                    bcd = set([rbin[1], rbin2[1]])
                    dih = set([rbin[2], rbin2[2]])
                    lengths = set([lbin, lbin2])

                    # Add length to combined bins if option enabled;
                    # product bins to get all combinations of bins
                    all_bins = [x, abc, bcd, dih]
                    if bin_length:
                        all_bins.append(lengths)
                    if clash:
                        all_bins.extend(clashes)

                    # Iterate through combinations of bins and store
                    for bin_12 in all_bins:
                        bin_12 = ' '.join(map(str, bin_12))
                        doc = {
                                'bin':bin_12, 
                                'name': name,
                                'idx1':idx1,
                                'idx2':idx2
                        }
                        # if check_dups:
                            # if len(list(bins.find(doc))) == 0:
                                # unsaved_docs.append(doc)
                        # else:
                        unsaved_docs.append(doc)

            # Update dataframe ever <interval> PDBs
            if i%interval == 0:
                bins = update(bins, start_time, unsaved_docs, interval, i)
                start_time = time.time()
                unsaved_docs = []

        bins = update(bins, start_time, unsaved_docs, interval, i, final=True)

        return bins

class HelixLookup(object):
    '''
    Class to handle binning and matching of helix databases. This maybe
    should be two classes, one for binning and one for matching, but
    this is it for now.
    '''

    def __init__(self, lookup_folder, query, name='unknown',
            verbose=False):
        self.verbose = verbose
        self.lookup_folder = lookup_folder
        self.query = query
        self.name = name

    def score_match(self, list_of_index_pairs):
        """
        Idea (idk where else to put this):
            To get 3rd, 4th, etc. helices, do a reverse lookup. That is,
            for each bin in the FOUND PDB, look for matches in the QUERY
            pdb.
        """
        # TO DO: score clashes
        return

    def submit_local(self, outdir):
        import glob
        lookups = sorted(glob.glob(self.lookup_folder + '/*.pkl'))
        print(self.lookup_folder)
        print(lookups)
        i = 0
        for lookup in lookups:
            print('MATCHING AGAINST {}'.format(lookup))
            out = os.path.join(outdir, '{}_results_{:03d}.pkl'.format(
                self.name, i)
                )
            self.match(pd.read_pickle(lookup), out=out)
            i += 1

    def submit_cluster(self, outdir, total_tasks):
        import glob
        lookups = sorted(glob.glob(self.lookup_folder + '/*.pkl'))
        task = int(os.environ['SGE_TASK_ID']) - 1
        out = os.path.join(outdir, '{}_results_{:03d}.pkl'.format(self.name,
            task))
        print('Saving to {}'.format(out))
        # Warning: total_tasks must be a multiple of len(lookups) for
        # now.
        increment = total_tasks // len(lookups)
        print('Increment {}'.format(increment))
        lookups_idx = task//increment
        print('Reading database file # {}'.format(lookups_idx))

        lookup = pd.read_pickle(lookups[lookups_idx])
        num_rows = lookup.shape[0]
        row_increment = num_rows // increment
        rowstart = (task%increment) * row_increment
        rowend = rowstart + row_increment
        lookup = lookup.iloc[rowstart:rowend]
        print('Looking up rows {} through {}'.format(rowstart, rowend))
        print(lookup)
        self.match(lookup, out=out)

    def match(self, lookup, out=None):
        names = []

        # Pandas rewrite
        print('Starting forward search...')
        for _bin, group in self.query.groupby(level='bin'):
            if self.verbose:
                print('Searching bin {}'.format(_bin))
            if _bin in lookup.index:
                for result in lookup.xs(_bin, level='bin').iterrows():
                    # xs results in (index, row) tuples; db is indexed by
                    # name, so row[0] is the name.
                    if self.verbose:
                        print('Matched to pdb {}'.format(result[0]))
                    names.append(
                            result[0]
                            )

        print('Forward search done.')

        print('Original name list:')
        print(names)
        min_matches = 2
        names = [item for item, count in
                collections.Counter(names).items() if
                count >= min_matches]

        print('All matches:')
        print(names)
        print(len(names))

        results = []
        # TEMP

        # sys.exit()
        i = 0
        for name in names:
            i += 1
            result = {}
            result['name'] = name
            print('-------------------------------------------------')
            print('Name: {}'.format(name))
            match = Match(name, self.query, lookup, verbose=self.verbose)
            match.find_edges()
            result['matches'] = match.max_subgraph()
            result['graph'] = match.graph
            results.append(result)
            # match.plot_graph()
            # print('searching {}'.format(name))
            # for _bin in self.binned.find({'name': name[0]}):
                # if _bin['idx1'] == name[1]:
                    # print('-------')
                    # print(_bin)
                    # for doc in self.query_bins.find({'bin':_bin['bin']}):
                        # print('MATCH:')
                        # results[name].append((doc['idx1'], doc['idx2']))
                        # print(doc)

        df = pd.DataFrame(results)
        if out:
            df.to_pickle(out)

        return df
        # for key in results:
            # print('------------------RESULTS FOR {}----------------'.format(
                            # key
                        # ))
            # for pair in set(results[key]):
                # print(pair)
        # for key in results:
            # print('PDB {} had {} matching transformations'.format(
                # key, len(set(results[key]))
                # ))


def test():
    import scan_helices

    test_path = 'test_files/6r9d.cif'
    init()
    pose = pose_from_file(test_path).split_by_chain(2)
    print(pose.size())
    scanner = scan_helices.PoseScanner(pose)
    helices = scanner.scan_pose_helices()
    helices = pd.DataFrame(helices)
    print(helices)
    helices = helices[helices['percent_exposed'] > 0.3]
    print(helices)
    print(helices.shape)
    print(helices['name'])

    # lookup = HelixLookup(pd.read_pickle('dataframes/final.pkl'),
            # query_df=helices, query_name='6r9d')
    lookup = HelixLookup(pd.DataFrame(),
            query_df=helices, query_name='6r9d', angstroms=5,
            # degrees=15, reset_querydb=True, dbname='nr')
            degrees=30, reset_querydb=True, dbname='test_bins')
    lookup.match()

def test_rifdock():
    import scan_helices

    test_path = 'test_files/test_rifgen/cluster_representatives/matchme.pdb'
    init()
    pose = pose_from_file(test_path)
    print(pose.size())
    scanner = scan_helices.PoseScanner(pose)
    helices = scanner.scan_pose_helices(split_chains=False,
            name='rifdock_test')
    helices = pd.DataFrame(helices)
    helices.to_pickle('rifdock_helices.pkl')
    sys.exit()
    print(helices)
    # helices = helices[helices['percent_exposed'] > 0.3]
    print(helices)
    print(helices.shape)
    print(helices['name'])

    # lookup = HelixLookup(pd.read_pickle('dataframes/final.pkl'),
            # query_df=helices, query_name='6r9d')
    lookup = HelixLookup(pd.DataFrame(),
            query_df=helices, query_name='6r9d', angstroms=2.5,
            degrees=15, reset_querydb=True, dbname='nr')
            # degrees=30, reset_querydb=True, dbname='test_bins')
    lookup.match()

def make_hash_table():
    print('Loading database and setting up lookup object...')
    # length cutoff of 2 turns or 10.8 angstroms
    lookup = HelixLookup(pd.read_pickle('nr_dataframes/final.pkl'),
            exposed_cutoff=0.3, length_cutoff=10.8, angstroms=2.5,
            degrees=15, dbname='nr')
    print('Done.')
    # binned = lookup.bin_db(lookup.df)
    lookup.update_bin_db()
    # out = "binned_0p3/last.pkl"
    # with open(out, 'wb') as f:
        # pickle.dump(binned, f)

def make_test_hash_table():
    client = MongoClient()
    deg=15
    angstroms=2.5
    # client['test_bins']['bins_{}A_{}D'.format(angstroms, deg)].drop()
    lookup=HelixLookup(pd.read_pickle('out.pkl'), exposed_cutoff=0.3,
            length_cutoff=10.8, angstroms=angstroms, degrees=deg,
            dbname='test_bins')
    lookup.update_bin_db()


def main():
    args = docopt.docopt(__doc__)
    dbpath = os.path.join(
            args['--database'],
            "bins_{}A_{}D".format(
                float(args['--angstroms']),
                float(args['--degrees'])
                )
            )
    if args['bin']:
        lookup = HelixBin(pd.read_pickle(args['<helix_dataframe>']),
                exposed_cutoff=0.3, length_cutoff=10.8,
                angstroms=float(args['--angstroms']),
                degrees=float(args['--degrees']), 
                verbose=args['--verbose'])
        lookup.bin_db(outdir=dbpath, bin_length=args['--length'])
    if args['match']:
        import scan_helices
        # Import pdb
        pdbfolder = args['<pdb_folder>']
        init()


        helicepath = os.path.join(pdbfolder, 'query_helices.pkl')
        if os.path.exists(helicepath):
            helices = pd.read_pickle(helicepath)
        else:
            all_helices = []
            import glob
            gz = glob.glob(pdbfolder + '/*.pdb.gz')
            dotpdb = glob.glob(pdbfolder + '/*.pdb')
            gz.extend(dotpdb)
            pdbs = sorted(gz)
            for path in pdbs:
                pose = pose_from_file(path).split_by_chain(1)

                # Scan pdb helices
                scanner = scan_helices.PoseScanner(pose)
                helices = scanner.scan_pose_helices(name='query',
                        split_chains=False, path=path)
                all_helices.extend(helices)
            helices = pd.DataFrame(all_helices)
            helices.to_pickle(helicepath)
        print("HELICES")
        print(helices)
        print(helices['vector'])

        # Bin pdb helices
        query = HelixBin(helices, exposed_cutoff=0.3,
                length_cutoff=10.8, 
                angstroms=float(args['--angstroms']), 
                degrees=float(args['--degrees']),
                verbose=args['--verbose'])
        query_bins = query.bin_db(bin_length=args['--length'])
        print('QUERY BINS')
        print(query_bins)

        # Match
        # name = os.path.basename(path).split('.')[0]
        name = 'query'
        print('Database:')
        print(dbpath)
        matcher = HelixLookup(dbpath, query_bins, name=name,
                verbose=args['--verbose'])
        if args['--tasks']:
            matcher.submit_cluster(args['--out'], int(args['--tasks']))
        else:
            matcher.submit_local(args['--out'])


if __name__=='__main__':
    # test()
    # test_rifdock()
    # make_hash_table()
    # make_test_hash_table()
    main()
