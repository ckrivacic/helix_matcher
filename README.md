# Helical ELements Interface eXplorer

HELIX is a protocol meant to help determine an optimal scaffold to bind to a target protein. It does this by (1) determining the optimal positions of helical elements on the surface of the target protein, (2) finding scaffolds from a database that enable the largest number of these helices, and (3) filtering down to those with the best metrics (lowest clash score, largest # of helices supported, lowest RIFDock score, etc.).

### Installation

Currently, HELIX relies on the C++ version of RIFDock to generate potential binding geometries for helices.
It also depends on PyRosetta.
Other dependencies can be found in the `helix.yml` conda environment.

Download the repo:
```bash
git clone https://github.com/ckrivacic/helix_matcher.git
cd helix_matcher
```

Install the conda environment (optional if you want to install dependencies yourself):
```bash
conda env create -f helix.yml
conda activate helix
```

Run setup.py:
```bash
python setup.py install
```

When prompted, you may choose to download the default database. This consists of vectorized helix information for all proteins in the nonredundant PDB, as well as a database of their relative orientations. Otherwise, you will have to create a database yourself from a folder of PDBs using the `helix scan` and `helix bin` commands.

### HELIX commands

Once installed, HELIX can be run using commands. These take the form of `helix <command_name> [arguments]`. Commands that are integral parts of the pipeline are numbered, e.g. 01_prep_rifdock, 02_rifdock, etc. Note that you do not need to provide the full name of the command, so for example `helix 01 [arguments]` will run the `01_prep_rifdock` command.

For a full list of available commands, see `helix --help`.


### Setting up a HELIX workspace

HELIX helps track settings, targets, inputs and outputs for your project. Begin by setting up a HELIX workspace:
```bash
helix setup <workspace_name>
```
You will be asked for a path to your base RIFDock folder and the Python binary that contains your dependencies (if you used the conda environment provided with the repo, this will be in `/path/to/anacdona/envs/helix/bin/python`). You can skip these by defining the relevant environmental variables.
You will be asked for a database and one or more paths to target PDB file(s). You may ignore the database prompt if you downloaded the default database. For the target PDB, provide one or more paths, space-separated (wildcard extension is accepted as well), to the PDB files of the proteins you wish to design a binder for.

### Generating potential helix binding geometries

Helix matcher uses RIFDock to find potential binding geometries. Because RIFDock tends to converge on 
a single solution when docking single helices, we dock the helices to all possible surface patches of 
the target protein. 

To set up a RIFDock sub-workspace, run `01_prep_rifgen` command:

```bash
helix 01 <workspace_path>
```

This will create a sub-directory called `rifdock_outputs`, and within there another subdirectory for each target provided during setup. Inside these target folders, there will be a folder for each surface patch on the target protein, with flags to be used with RIFGen.

Next, run RIFDock:

```bash
helix 02 <workspace_path>```
```


