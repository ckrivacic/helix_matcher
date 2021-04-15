# Helix Matcher

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
python3 rifgen.py <folder> [options]
```

Where `<folder>` is the same as `<output_folder>` above.

Possible options are:
```python
--sge # If you want to use an SGE job distributor to run RIFGEN
--tasks=NUM # If using SGE, how many tasks do you want to split RIFGEN into?
```

Once RIFGen is complete, run RIFDock:

```bash
python3 run_rifdock.py <folder> [options]
```

Possible options are:
```python
--sge # Same as above; use this option if running on Wynton
--tasks=NUM # Same as above; how many tasks to split into if using the cluster
```

