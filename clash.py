import alphashape
from shapely.ops import triangulate
import prody
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
# from descartes import PolygonPatch



def get_alphashape(pdb, chain='B'):
    '''
    Returns an AlphaShape object of a pdb file, outlining its general
    shape for use in a clash filter.
    '''
    atoms = prody.parsePDB(pdb, chain=chain)
    # For some reason there is a level which is not populated. This may
    # be for proteins w/ multiple chains.
    coords = atoms.getCoordsets()[0]
    # coords = [(0., 0.), (0., 1.), (1., 1.), (1., 0.), (0.5, 0.5)]

    alpha_shape = alphashape.alphashape(coords, 0.05)
    print(type(alpha_shape))

    fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    ax = Axes3D(fig)
    ax.scatter(*zip(*coords))
    # ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
    ax.plot_trisurf(*zip(*alpha_shape.vertices),
            triangles=alpha_shape.faces)
    plt.show()


def main():
    pdb = 'test_files/boundary/cluster_representatives/clst10_rep.pdb.gz'
    get_alphashape(pdb, chain='B')

main()
