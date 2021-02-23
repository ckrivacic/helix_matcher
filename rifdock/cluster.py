import gzip
import glob
import sys, os
import numpy as np
from shutil import copyfile


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
        if len(glob.glob(os.path.dirname(self.path) + '/all.d*')) > 1:
            print('WARNING: Multiple docking scorefiles detected.')
            print('This may result in the wrong score being used.')
        scorefile = os.path.join(os.path.dirname(self.path), 'all.dok')
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
    def __init__(self, folder, threshold=None):
        self.threshold=threshold
        self.designs = []
        # Path sould be to parent rifdock directory.
        for path in glob.glob(folder + '/*/docked_full/*.pdb.gz'):
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



        # Cluster design such that no two designs in a cluster are
        # further apart than the threshold, in Angstroms.

        dist_vector = sp_dist.squareform(dist_matrix)
        mean_dist = np.mean(dist_vector)
        hierarchy = sp_clust.complete(dist_vector)
        print(self.threshold)
        print(mean_dist)
        print(self.threshold or mean_dist)
        sys.exit()
        clusters = sp_clust.fcluster(
                hierarchy, self.threshold or mean_dist, criterion='distance')
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
    clust = StructureCluster(sys.argv[1], threshold=10)
    clust.cluster_coords(verbose=True)
    clusters = clust.clusters

    outpath = os.path.join(sys.argv[1], 'cluster_representatives')
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    overwrite = True
    if overwrite:
        for f in glob.glob(sys.argv[1] +
                '/cluster_representatives/*.pdb.gz'):
            os.remove(f)


    for clst in clusters:
        print('CLUSTER {}'.format(clst))
        print('REP:')
        print(clusters[clst].rep.path)
        copyfile(clusters[clst].rep.path, os.path.join(outpath,
            'clst{}_rep.pdb.gz'.format(clst)))
