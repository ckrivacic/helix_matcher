'''
Create bins or match a query protein.

Usage:
    matcher.py bin <helix_dataframe> [options]
    matcher.py bin_query <match_workspace> [options]
    matcher.py match <match_workspace> [options]

options:
    --local, -l  Run locally

    --tasks=NUM, -j  Run on the cluster using SGE. Argument should be # of
    tasks per dataframe.

    --length, -e  Bin by length

    --verbose, -v  Verbose output

    --database=PATH, -d  Database of relative helix orientations  
    [default: database/]

    --out=PATH, -o  Where to save outputs  [default: .]

    --angstroms=NUM, -a  Binning option. How fine should the distance bins
    be?  [default: 2.5]

    --degrees=NUM, -g  Binning option. How fine should the angle bins be?
    [default: 15]

    --settings=YML, -s  Provide a settings file.

    --scaffold=PDB  Only run matching for a given helix length/RIFDock
    scaffold.

    --min-dist=FLOAT  Minimum distance for two helices to have their relative orientations saved.
    Useful for query dataframes where you have lots of helices.

    --max-dist=FLOAT  Maximum distance before two helices do not have their relative orientations saved.

    --overwrite  Overwrite relative orientation dataframe even if another task has already made it

    --bin-tasks=INT  For binning query helices only. How many tasks to split the binning process into.

    --bin-task-id=INT  Task ID for binning (for local runs only)
'''
import collections
import os, psutil, sys
import pickle
import subprocess
import docopt
import numpy as np
import pandas as pd
import networkx as nx
from helix import workspace as ws
from helix.matching.scan_helices import final_vector
from helix.utils import numeric
from helix.utils import utils
from itertools import product
from pyrosetta import init, pose_from_file
# import graph_tool.all as gt


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
    '''
    inds = np.digitize(array, bins)
    binned = tuple([bins[inds[n]-1] for n in range(array.size)])
    return binned


def relative_position(row1, row2, vectortype='normalized_vector'):
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

    norm_v1 = row1[vectortype]
    norm_v2 = row2[vectortype]
    centroid_dist = numeric.euclidean_distance(row1['centroid'],
            row2['centroid'])
    abc = numeric.angle(norm_v1[0], norm_v1[1], norm_v2[0])
    bcd = numeric.angle(norm_v1[1], norm_v2[0], norm_v2[1])
    dihedral = numeric.dihedral(norm_v1[0], norm_v1[1], norm_v2[0],
            norm_v2[1])
    # plot_vectors([norm_v1, norm_v2], color='black')

    return centroid_dist, abc, bcd, dihedral


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
        # import graph_tool.draw as draw
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
            query_df=None, query_name=None, angstroms=2.5, degrees=15,
            verbose=False, start=None, stop=None, min_distance=None, max_distance=None,
            task=None):
        self.verbose = verbose
        self.df = helix_db
        self.df['idx'] = self.df.index
        if task:
            self.task = int(task)
        else:
            self.task = 1

        if min_distance:
            self.min_distance = float(min_distance)
        else:
            self.min_distance = min_distance
        if max_distance:
            self.max_distance = float(max_distance)
        else:
            self.max_distance = max_distance

        # Binning parameters
        self.degrees = degrees
        self.angstroms = angstroms
        self.setup_bins()
        binned_name = 'bins_{}A_{}D'.format(self.angstroms,
                self.degrees)


        # Trimming dataframe
        if length_cutoff:
            self.df = self.df[self.df['length'] > length_cutoff]
        if exposed_cutoff:
            self.df = self.df[self.df['percent_exposed'] >
                    exposed_cutoff]
        if 'normalized_vector' not in self.df.columns:
            self.df['normalized_vector'] = self.df.apply(lambda x:
                    final_vector(x['direction'], 1, x['centroid']), axis=1)

        if start:
            self.start = start
        else:
            self.start = 0
        if stop:
            self.stop = stop
        else:
            self.stop = self.df.shape[0]

    def setup_bins(self):
        nrbins = int(360//self.degrees) + 1
        self.rbins = np.linspace(-180, 180, nrbins)
        tstart = -10000
        tstop = 10000
        ntbins = int((tstop - tstart) // self.angstroms) + 1
        self.tbins = np.linspace(tstart, tstop, ntbins)

    def bin_db(self, outdir=None, bin_length=False):
        '''
        Bin dataframes.
        '''

        from scipy.spatial.transform import Rotation as R
        import subprocess
        import time

        # db = self.client[dbname]
        # bins = db['bins_{}A_{}D'.format(
            # self.angstroms, self.degrees
            # )]
        bins = pd.DataFrame(columns=['bin', 'name', 'idx1', 'idx2'])
        # Pandas indices are hash lookups and we can have multiple of
        # them, but they cannot be added piecewise. Therefore we will
        # create partial tables, then create the indices and save the
        # dataframes. Results will be saved in chunks.
        # bins.set_index(['bin', 'name'], inplace=True)
        total_proteins = len(set(self.df['name']))
        if total_proteins == 1:
            interval = 5000
        else:
            interval = 500

        # import shelve

        # binned = shelve.open('binned_0p3/hashtable', 'c', writeback=True)
        # i tracks # of names analyzed
        i = 0
        # j tracks # of rows analyzed
        j = 0
        # saveno tracks how many dataframes have been saved.
        self.saveno = 1
        unsaved_docs = []
        start_time = time.time()

        def update(bins, start_time, unsaved_docs, interval, i, j,
                final=False):
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
            if total_proteins == 1:
                total_combos = self.df.shape[0]**2
                remaining = (total_combos - j) / rate / 3600
                print('Analysis of {} rows took {:02f} seconds. Est. {:02f} h remaining'.format(
                    interval, elapsed, remaining
                ))
            else:
                remaining = (total_proteins - i) / rate / 3600
                print('Analysis of {} pdbs took {:02f} seconds. Est. {:02f} h remaining'.format(
                    interval, elapsed, remaining
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
                    if self.task:
                        outfile = 'bins_{}A_{}D_{:03d}_{:04d}.pkl'.format(self.angstroms,
                                self.degrees, self.task, self.saveno,)
                    else:
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

        groups = self.df.groupby(['name'])
        names = sorted(list(groups.groups.keys()))
        # if self.start:
        #     names = names[self.start:]
        # if self.stop:
        #     names = names[:self.stop]
        for name in names:
            group = groups.groups[name]
            group_df = self.df.loc[group]
            i += 1

            # for combination in product(self.df.loc[group].T.to_dict().values(),
            #         repeat=2):
            if total_proteins > 1:
                idx1_range = group_df.index
                idx2_range = group_df.index
            else:
                idx1_range = range(self.start, self.stop)
                idx2_range = range(0, group_df.shape[0])
            for idx1 in idx1_range:
                for idx2 in idx2_range:
                    # if self.verbose:
                    #     print('Indices are ({}, {})'.format(idx1, idx2))
                    #     print('Group:')
                    #     print(group_df)
                    # if combination[0]['idx'] == combination[1]['idx']:
                    j += 1
                    if idx1 == idx2:
                        continue
                    if idx1 not in self.df.index or idx2 not in self.df.index:
                        continue
                    # vector1 = combination[0]['vector']
                    # vector2 = combination[1]['vector']

                    # plot_vectors([vector1, vector2], color='purple')

                    # idx1 = combination[0]['idx']
                    # idx2 = combination[1]['idx']
                    # if self.verbose:
                        # print('------------------------------------')
                        # print(combination[0])
                        # print(combination[1])
                    # if total_proteins > 1:
                    #     combination = (group_df.iloc[idx1], group_df.iloc[idx2])
                    # else:
                    combination = (group_df.loc[idx1], group_df.loc[idx2])

                    dist, angle1, angle2, dihedral =\
                            relative_position(combination[0], combination[1])
                    dist = np.array([dist])
                    if self.min_distance:
                        if dist < self.min_distance:
                            # Add the indices to a list of indices we don't need to check? Hm.... Does that actually help?
                            # Would only save us a single distance calculation.
                            continue
                    if self.max_distance:
                        if dist > self.max_distance:
                            continue
                    angles = np.array([angle1, angle2, dihedral])

                    lengths = np.array([combination[0]['length'],
                        combination[1]['length']])
                    lbin = bin_array(lengths, self.tbins)
                    lbin2 = bin_array(lengths, self.tbins +
                            (self.angstroms/2))

                    rbin = bin_array(angles, self.rbins)
                    tbin = bin_array(dist, self.tbins)
                    rbin2 = bin_array(angles, self.rbins + (self.degrees/2))
                    tbin2 = bin_array(dist, self.tbins +
                            (self.angstroms/2))

                    x = [tbin[0], tbin2[0]]
                    abc = [rbin[0], rbin2[0]]
                    bcd = [rbin[1], rbin2[1]]
                    dih = [rbin[2], rbin2[2]]
                    lengths = [lbin, lbin2]

                    if bin_length:
                        all_bins = product(x, abc, bcd, dih, lengths)
                    else:
                        all_bins = product(x, abc, bcd, dih)

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
                    if total_proteins == 1:
                        if j%interval == 0:
                            bins = update(bins, start_time, unsaved_docs, interval, i, j)
                            start_time = time.time()
                            unsaved_docs = []

            if total_proteins > 1:
                if i%interval == 0:
                    bins = update(bins, start_time, unsaved_docs, interval, i, j)
                    start_time = time.time()
                    unsaved_docs = []

        bins = update(bins, start_time, unsaved_docs, interval, i, j, final=True)

        return bins

class HelixLookup(object):
    '''
    Class to handle binning and matching of helix databases. This maybe
    should be two classes, one for binning and one for matching, but
    this is it for now.
    '''

    def __init__(self, lookup_folder, query_list, name='unknown',
            verbose=False):
        self.verbose = verbose
        self.lookup_folder = lookup_folder
        self.query_list = query_list
        self.name = name

    def submit_local(self, outdir):
        import glob
        lookups = sorted(glob.glob(self.lookup_folder + '/*.pkl'))
        print(self.lookup_folder)
        print(lookups)
        i = 0
        os.makedirs(outdir, exist_ok=True)
        for lookup in lookups:
            print('MATCHING AGAINST {}'.format(lookup))
            for query in queries:
                print(f'USING QUERY DF {query}')
                self.query = utils.safe_load(query)
                out = os.path.join(outdir, '{}_results_{:03d}.pkl'.format(
                    self.name, i)
                    )
                self.match(utils.safe_load(lookup), out=out)
                i += 1

    def submit_cluster(self, outdir, tasks):
        # For now, "queries" is the list of query relative orientation dataframes.
        # For now just pass None to "query" when initiating this class if using submit_cluster.
        import glob
        lookups = sorted(glob.glob(self.lookup_folder + '/*.pkl'))
        total_tasks = tasks * len(lookups) * len(self.query_list)
        tasks_per_query = total_tasks // len(self.query_list)
        query_increment = tasks * len(lookups)

        task = int(os.environ['SGE_TASK_ID']) - 1
        self.query = utils.safe_load(queries[task // query_increment])
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, '{}_results_{:03d}.pkl'.format(self.name,
            task))
        print('Results will be saved to {}'.format(out))

        # Warning: total_tasks must be a multiple of len(lookups) for
        # now.
        increment = total_tasks // tasks_per_query //len(lookups)
        print('Increment {}'.format(increment))
        lookups_idx = task%increment
        print('Reading database file # {}'.format(lookups_idx))

        lookup = utils.safe_load(lookups[lookups_idx])
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

        # print('Original name list:')
        # print(names)
        min_matches = 2
        names = [item for item, count in
                collections.Counter(names).items() if
                count >= min_matches]

        print('All matches:')
        # print(names)
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
    # import scan_helices
    from helix.matchign import scan_helices

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
    from helix.matching import scan_helices

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
    lookup = HelixLookup(utils.safe_load('nr_dataframes/final.pkl'),
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
    lookup=HelixLookup(utils.safe_load('out.pkl'), exposed_cutoff=0.3,
            length_cutoff=10.8, angstroms=angstroms, degrees=deg,
            dbname='test_bins')
    lookup.update_bin_db()


def main():
    args = docopt.docopt(__doc__)
    print(args)

    if args['--settings']:
        # Deprecated; settings handled by submission command
        import yaml
        runtype = 'bin' if args['bin'] else 'match'
        settings = yaml.load(open(args['--settings'], 'r'))
        print(settings)
        for option in settings[runtype]:
            args[option] = settings[runtype][option]
        print(args)

    dbpath = os.path.join(
            args['--database'],
            "bins_{}A_{}D".format(
                float(args['--angstroms']),
                float(args['--degrees'])
                )
            )
    if args['bin']:
        lookup = HelixBin(utils.safe_load(args['<helix_dataframe>']),
                exposed_cutoff=0.3, length_cutoff=10.8,
                angstroms=float(args['--angstroms']),
                degrees=float(args['--degrees']),
                verbose=args['--verbose'], min_distance=args['--min-dist'],
                max_distance=args['--max-dist'],)
        lookup.bin_db(outdir=dbpath, bin_length=args['--length'])
    if args['match'] or args ['bin_query']:
        # import scan_helices
        from helix.matching import scan_helices
        workspace = ws.workspace_from_dir(args['<match_workspace>'])
        # Import pdb
        '''
        # This block is for working with CLUSTER folders. Uncomment if I decide to go back to using clusters.
        if args['--scaffold']:
            pdbfolders = [workspace.scaffold_clusters(args['--scaffold'])]
        else:
            pdbfolders = workspace.all_scaffold_clusters
        '''
        # This block is for NOT using clusters.
        pdbfolders = [workspace.cluster_outputs]

        init()

        if not args['--scaffold'] and \
                os.path.exists(workspace.all_scaffold_dataframe):
            all_helices = utils.safe_load(workspace.all_scaffold_dataframe)
        else:
            all_helices = []
            for pdbfolder in pdbfolders:
                # helicepath = os.path.join(pdbfolder, 'query_helices.pkl')
                helicepath = workspace.query_dataframe
                if os.path.exists(helicepath):
                    helices = utils.safe_load(helicepath)
                    all_helices.append(helices)
                else:
                    folder_helices = []
                    import glob
                    gz = glob.glob(pdbfolder + '/*.pdb.gz')
                    dotpdb = glob.glob(pdbfolder + '/*.pdb')
                    gz.extend(dotpdb)
                    pdbs = sorted(gz)
                    for path in pdbs:
                        # For PatcHMAN, second chain is the docked helix
                        pose = pose_from_file(path).split_by_chain(2)
                        path = os.path.relpath(path,
                                start=workspace.root_dir)

                        # Scan pdb helices
                        scanner = scan_helices.PoseScanner(pose)
                        helices = scanner.scan_pose_helices(name='query',
                                split_chains=False, path=path)
                        folder_helices.extend(helices)
                    helices = pd.DataFrame(folder_helices)
                    helices.to_pickle(helicepath)
                    all_helices.append(helices)
            all_helices = pd.concat(all_helices, ignore_index=True)
            if not args['--scaffold']:
                # Don't save to the all_scaffold path if not using all
                # scaffolds
                all_helices.to_pickle(workspace.all_scaffold_dataframe)

        print("HELICES")
        print(all_helices)

        if not os.path.exists(workspace.query_database_dir):
            os.makedirs(workspace.query_database_dir, exist_ok=True)
        if len(workspace.relative_orientation_dataframes) < 1 or args['--overwrite']:
            if 'SGE_TASK_ID' in os.environ:
                task_id = int(os.environ['SGE_TASK_ID'])
            else:
                if args['--bin-task-id']:
                    task_id = int(args['--bin-task-id'])
                else:
                    task_id = 1
            if args['--bin-tasks']:
                interval = (all_helices.shape[0] // int(args['--bin-tasks'])) + 1
                task = task_id - 1
                start = task * interval
                stop = start + interval
            else:
                task_id = 1
                start = 0
                stop = all_helices.shape[0]
            # Bin pdb helices
            query = HelixBin(all_helices, exposed_cutoff=0.3,
                    length_cutoff=10.8,
                    angstroms=float(args['--angstroms']),
                    degrees=float(args['--degrees']),
                    verbose=args['--verbose'],
                    min_distance=args['--min-dist'],
                    max_distance=args['--max-dist'],
                    start=start, stop=stop, task=task_id)
            query_bins = query.bin_db(bin_length=args['--length'], outdir = workspace.query_database_dir,)
            print('QUERY BINS')
            print(query_bins)
            # if not os.path.exists(workspace.relative_orientation_dataframe) and not args['--overwrite']:
            #     query_bins.to_pickle(workspace.relative_orientation_dataframe)
        # else:
        #     query_bins = workspace.relative_orientations

        if args['bin_query']:
            # Exit after binning
            sys.exit()
        query_bins = workspace.relative_orientation_dataframes

        # Match
        # name = os.path.basename(path).split('.')[0]
        name = 'query'
        print('Database:')
        print(dbpath)
        matcher = HelixLookup(dbpath, query_bins, name=name,
                verbose=args['--verbose'])
        if args['--local']:
            matcher.submit_local(workspace.output_dir)
        elif args['--tasks']:
            matcher.submit_cluster(workspace.output_dir, int(args['--tasks']))
        else:
            matcher.submit_cluster(workspace.output_dir, 1)


if __name__=='__main__':
    # test()
    # test_rifdock()
    # make_hash_table()
    # make_test_hash_table()
    main()
