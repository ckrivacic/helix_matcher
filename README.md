# Helical ELements Interface eXplorer

### Installation

Currently, HELIX relies on the C++ version of RIFDock to generate potential binding geometries for helices.
It also depends on PyRosetta.
Other dependencies can be found in the `helix.yml` conda environment.

Download the repo:
```
git clone https://github.com/ckrivacic/helix_matcher.git
cd helix_matcher
```

Install the conda environment (optional if you want to install dependencies yourself):
```
conda env create -f helix.yml
conda activate helix
```

Run setup.py:
```
python setup.py install
```

When prompted, you may want to download the default database. This consists of vectorized helix information for all proteins in the nonredundant PDB, as well as a database of their relative orientations. Otherwise, you will have to create a database yourself from a folder of PDBs using the `helix scan` and `helix bin` commands.


### Generating potential helix binding geometries

Helix matcher uses RIFDock to find potential binding geometries. Because RIFDock tends to converge on 
a single solution when docking single helices, we dock the helices to all possible surface patches of 
the target protein. 

To set up a RIFDock workspace, run `prep_rifgen.py`:

```bash
python3 rifgen.py <target_pdb> <output_folder>
```

This will create a sub-directory in output_folder for each surface patch on the target protein, with flags to be used with RIFGen.

Next, run RIFGen:

```bash
# You can skip the first command in the future if you add it to your .bashrc file.
export RIFGEN=/path/to/rifgen_binary
python3 rifgen.py <folder>
```

Where `<folder>` is the same as `<output_folder>` above.

You can also run RIFGEN on the cluster:
```bash
./submit_rifgen.sh
```
Currently you will need to edit this script to point to your python path, the `<folder>` used above, and to set the number of tasks. Streamlining this is on the to-do list.

Once RIFGen is complete, run RIFDock:

```bash
python3 run_rifdock.py <folder>
```


