# Helical ELements Interface eXplorer

HELIX is a protocol meant to help determine an optimal scaffold to bind to a target protein. It does this by (1) determining the optimal positions of helical elements on the surface of the target protein, (2) finding scaffolds from a database that enable the largest number of these helices, and (3) filtering down to those with the best metrics (lowest clash score, largest # of helices supported, lowest RIFDock score, etc.).

### Installation

If you wish to use RIFDock to place fragments, HELIX relies on the C++ version of [RIFDock](https://github.com/rifdock/rifdock) to generate potential binding geometries for helices.

[PatchMAN](https://github.com/Alisa-Kh/PatchMAN) and its dependencies are required to place fragments based on tertiary motifs (recommended).

HELIX also requires [PyRosetta](https://www.pyrosetta.org/) and [klab](https://github.com/kortemme-lab/klab).
Other dependencies can be found in the `helix.yml` conda environment.

Download the repo:
```bash
git clone https://github.com/ckrivacic/helix_matcher.git
cd helix_matcher
```

Install the conda environment (optional if you want to install the dependencies not outlined above yourself):
```bash
conda env create -f helix.yml
conda activate helix
```

Run setup.py:
```bash
python setup.py install
```

When prompted, you may choose to download the default database. This consists of vectorized helix information for all proteins in the nonredundant PDB, as well as a database of their relative orientations. Otherwise, you will have to create a database yourself from a folder of PDBs using the `helix scan` and `helix bin` commands.


### Setting up a HELIX workspace

HELIX helps track settings, targets, inputs and outputs for your project. Begin by setting up a HELIX workspace:
```bash
helix setup <workspace_name>
```
You will be asked for a path to your base RIFDock folder and the Python binary that contains your dependencies (if you used the conda environment provided with the repo, this will be in `/path/to/anacdona/envs/helix/bin/python`). You can skip these by defining the relevant environmental variables.
You will be asked for a database and one or more paths to target PDB file(s). You may ignore the database prompt if you downloaded the default database. For the target PDB, provide one or more paths, space-separated (wildcard extension is accepted as well), to the PDB files of the proteins you wish to design a binder for.

### HELIX commands

Once installed, HELIX can be run using commands. These take the form of `helix <command_name> [arguments]`. Commands that are integral parts of the pipeline are numbered, e.g. 01_prep_rifdock, 02_rifdock, etc. Note that you do not need to provide the full name of the command, so for example `helix 01 [arguments]` will run the `01_prep_rifdock` command.

For a full list of available commands, see `helix --help`.

The recommended usage is as follows. For any command, pass the `--help` option to get more info.
```bash
helix setup <workspace_path>
# You will provide the input target and a few necessary paths to HELIX.
helix prep_patchman <workspace_path>
# Splits the target into surface patches in preparation for running PatchMAN.
helix patchman <workspace_path>
# Runs PatchMAN. Recommended that you run this on a cluster with SGE, but if you want to run it locally, do so by passing the -l option.
helix design_patchman <workspace_path>
# Designs the docked helices (ignoring anything <60% helical). Again, if running locally pass -l.
helix filter <workspace>
# Filters docked fragments. You may provide your own custom filters via a yaml-formated file (formatting described below) 
# by passing --filter <path/to/filter.yml>
helix bin_docked_helices <workspace_path>
# Calculates helix vectors and bins their relative orientations. Use -l to run locally.
helix 04_match <workspace_path>
# Matches the binned helices against your database. Use -l to run locally.
helix 05_score_matches <workspace_path>
# Scores the matches. Use -l to run locally.
helix export_matches <workspace_path>
# Exports matches as match-target complexes if they pass a clash score and RMSD (versus the docked fragments) threshold. This threshold 
# is currently hardcoded in commands/export_matches.py, where you can edit the values on lines 201-202 to get the desired number of 
# designable complexes.
helix 06_design_scaffolds <workspace_path>
# Designs the scaffolds. By default, transfers residues over from the docked fragments if they are within 1.5A and prevents repacking for 
# any that it is able to pack well in an initial round of design. Pass --special-rot to instead allow designs for these positions, but 
# with a score bonus for choosing the same rotamer as the docked fragment. Pass -l to run locally (not recommended generally, but 
# especially not here). By default creates 10 sequences per input, but you can change this by passing -n.
helix plot_designs <workspace> <plot_type> [options]
# Makes plots of your designs. Plot types implemented right now are "scatter" and "violin". If you make a scatterplot of just one target, 
# you can click on points to open a PyMOL session of the design. You can also open multiple designs at once by pressing <Shift>+A, then 
# clicking on a few points, then pressing <Shift>+C to open. If you click on a point where there are multiple designs, your terminal will 
# list the possible designs, and you can choose one by pressing the corresponding number (the plot needs to be the focus window) and then 
# pressing <Enter>.
```

At this stage, you should check your filtered output (located in `<worspace_path>/rifdock_outputs/<target_name>/filtered/`, to make sure you have enough helices (and not so many that matching will be intractable). 1-4 helices per protein residue seems to work well.

To provide HELIX a custom filter file, use the following formatting. These filters can be provided to the `filter`, `analyze_designs`, and `plot_designs` commands.
```
'threshold': {
  'metric_1': ['<', 30], # Any fragment with a value for this metric less than 30 passes filter
  'metric_2': ['>', 0.3], # Any fragment with a value for this metric greater than 0.3 passes filter
  }
'percentil': {
  'metric_3': ['<', 0.2], # Any fragment in the bottom 20th percentile for this metric passes filter
  'metric_4': ['>', 0.6], # Any fragment in the top 60th percentil for this metric passes filter
  }
```


