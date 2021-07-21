'''
Aligns all RIFDock outputs to the input PDB, then clusters them based on
structure similarity. Representatives of each cluster are picked based
on their RIFDock score.

Usage:
    cluster_sge.py <rif_workspace> [options]

Options:
    --task=NUM, -t  Provide a task number. Each task is a different
    scaffold.
'''
import gzip
import docopt
import glob
import sys, os
import numpy as np
from shutil import copyfile
from helix import workspace as ws


class Design(object):
    def __init__(self, path):
        self.path = path
        self.backbone_coords = []
        self.score = None
        self.get_score()

    def read_coords(self, chain=None):
        with gzip.open(self.path, mode='rt') as f:
            lines = f.readlines()

        coords = []
        atoms = 'N', 'CA', 'C'

        for line in lines:
            if not line.startswith('ATOM'): continue

            atom_name = line[13:16].strip().upper()
            residue_id = int(line[22:26])
            chain_id = line[21]
            atom_coord = np.array((
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])))

            if atom_name in atoms:
                if chain:
                    if chain_id == chain:
                        coords.append(atom_coord)
                else:
                    coords.append(atom_coord)

        self.backbone_coords = np.array(coords)

    def get_score(self):
        # if len(glob.glob(os.path.dirname(self.path) + '/all.d*')) > 1:
            # print('WARNING: Multiple docking scorefiles detected.')
            # print('This may result in the wrong score being used.')
        scorefiles = sorted(glob.glob(os.path.join(
            os.path.dirname(self.path), 'all.dok*'
            )))
        # scorefile = os.path.join(os.path.dirname(self.path), 'all.dok')
        scorefile = scorefiles[-1]
        with open(scorefile, 'r') as f:
            for line in f:
                score = line[76:82]
                path = line[141:-1]

                if os.path.basename(path) == os.path.basename(self.path):
                    self.score = float(score)
                    return


class Cluster(object):
    def __init__(self, cluster):
        self.cluster = cluster
        self.designs = []
        return

    @property
    def rep(self):
        return self.find_rep()

    def add(self, design):
        self.designs.append(design)

    def find_rep(self):
        lowest_score = 9999
        for design in self.designs:
            if design.score < lowest_score:
                lowest_score = design.score
                rep = design

        return rep


class StructureCluster(object):
    def __init__(self, workspace, basename=None, threshold=None):
        self.threshold=threshold
        self.designs = []
        # Path sould be to parent rifdock directory.
        if basename:
            search = workspace.focus_dir + '/*/docked_full/{}*.pdb.gz'.format(basename)
        else:
            search = workspace.focus_dir + '/*/docked_full/*.pdb.gz'
        for path in glob.glob(
                search
                ):
            design = Design(path)
            design.read_coords(chain='A')
            self.designs.append(design)

        self.clusters = {}


    def cluster_coords(self, verbose=False):
        import scipy.spatial.distance as sp_dist
        import scipy.cluster.hierarchy as sp_clust
        from itertools import combinations

        num_designs = len(self.designs)
        if num_designs < 2: return

        # Calculate pariwise distance matrix

        dist_matrix = np.zeros((num_designs, num_designs))
        design_combos = combinations(enumerate(self.designs), 2)

        for (i, design_i), (j, design_j) in design_combos:
            dist_matrix[i,j] = self.calculate_rmsd(design_i,
                    design_j)
            dist_matrix[j,i] = dist_matrix[i,j]

        print(dist_matrix)



        # Cluster design such that no two designs in a cluster are
        # further apart than the threshold, in Angstroms.

        dist_vector = sp_dist.squareform(dist_matrix)
        mean_dist = np.mean(dist_vector)
        hierarchy = sp_clust.complete(dist_vector)
        print(self.threshold)
        print(mean_dist)
        print(self.threshold or mean_dist)
        clusters = sp_clust.fcluster(
                hierarchy, self.threshold or mean_dist, criterion='distance')
        print(clusters)
        # clusters = sp_clust.fcluster(
                # hierarchy, 20,
                # criterion='maxclust')

        for cluster, design in zip(clusters, self.designs):
            design.structure_cluster = cluster
            if cluster not in self.clusters:
                clst = Cluster(cluster)
                clst.add(design)
                self.clusters[cluster] = clst
            else:
                self.clusters[cluster].add(design)

        # Print some debugging information, if requested.

        if verbose == True:
            cluster_map = {}

            for cluster, design in zip(clusters, self.designs):
                cluster_map.setdefault(cluster, []).append(design)

            for cluster in sorted(set(clusters)):

                # Print out the designs contained in this cluster.

                print("Cluster {}:".format(cluster))
                for design in cluster_map[cluster]:
                    print(" ", design.path)
                print()

                # Print out pairwise distances for every cluster member.

                X = cluster_map[cluster]
                N = len(X)
                D = np.zeros((N, N))

                for i, design_i in enumerate(X):
                    for j, design_j in enumerate(X):
                        D[i,j] = self.calculate_rmsd(design_i, design_j)

                print(sp_dist.squareform(D))
                print()

                # Offer to display the cluster in pymol.

                command = ['pymol', '-qx', '-d', 'as cartoon']
                for design in cluster_map[cluster]:
                    command.append(design.path)

                # if raw_input("  View in pymol? [y/N] ") == 'y':
                    # import subprocess
                    # subprocess.check_output(command)

                print()

    def calculate_rmsd(self, design_1, design_2):
        # assert len(design_1.backbone_coords) and len(design_2.backbone_coords)
        difference = design_1.backbone_coords - design_2.backbone_coords
        num_atoms = design_1.backbone_coords.shape[0]

        return np.sqrt(np.sum(difference**2) / num_atoms)


if __name__=='__main__':
    # folders = sorted(glob.glob(sys.argv[1] + '/*_output'))
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<rif_workspace>'])
    assert(type(workspace).__name__ == 'RIFWorkspace')

    if 'SGE_TASK_ID' in os.environ:
        task = int(os.environ['SGE_TASK_ID']) - 1
    elif args['--task']:
        task = int(args['--task']) - 1
    else:
        task = 0
    # targets = workspace.targets
    helices = workspace.helices
    helix = helices[task]
    # folders = workspace.patches
    # workspace = ws.RIFWorkspace(workspace.root_dir, targets[task])

    # for helixlength in [3,4,6,8]:
    # for helix in workspace.helices:
    basename = workspace.basename(helix)
    clust = StructureCluster(workspace, basename=basename, threshold=10)
    print('Clustering...')
    clust.cluster_coords(verbose=True)
    clusters = clust.clusters

    outpath = os.path.join(workspace.focus_dir, 'cluster_representatives',
            basename)
    # if not os.path.exists(outpath):
        # print('PATH NO EXIST')
        # os.mkdir(outpath)
    os.makedirs(outpath, exist_ok=True)
    scores = open(os.path.join(
        outpath,
        '{}.scores'.format(basename)), 'w')

    overwrite = True
    # overwrite = False


    for clst in clusters:
        outfile = '{}_clst{}_rep.pdb.gz'.format(basename,
                clst)
        out = os.path.join(outpath, outfile)
        if overwrite:
            if os.path.exists(out):
                os.remove(out)
        print('CLUSTER {}'.format(clst))
        print('REP:')
        print(clusters[clst].rep.path)
        rep_dir = os.path.dirname(
                os.path.abspath(
                    clusters[clst].rep.path
                    )
                )
        relpath = os.path.relpath(
                rep_dir,
                outpath
                )
        # copyfile(clusters[clst].rep.path, os.path.join(outpath,
            # '{}turn_clst{}_rep.pdb.gz'.format(helixlength, clst)))
        os.symlink(
                os.path.join(relpath,
                    os.path.basename(clusters[clst].rep.path)),
                out
                )
        dokfile = sorted(glob.glob(rep_dir + '/*.dok*'))[-1]
        with open(dokfile, 'r') as f:
            for line in f:
                filename = os.path.basename(line.split(' ')[-1]).strip('\n')
                if filename == os.path.basename(
                                clusters[clst].rep.path):
                    # line.append(os.path.basename(out))
                    scores.write(line)

    scores.close()
