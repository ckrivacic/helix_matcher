"""\
Query the user for all the input data needed for a design.  This includes a 
starting PDB file, the backbone regions that will be remodeled, the residues 
that will be allowed to design, and more.  A brief description of each field is 
given below.  This information is used to build a workspace for this design 
that will be used by the rest of the scripts in this pipeline.  

Usage:
    helix setup_workspace <workspace> [--remote] [--overwrite]

Options:
    --remote, -r
        Setup a link to a design directory on a remote machine, to help with 
        transferring data between a workstation and a cluster.  Note: the 
        remote and local design directories must have the same name.

    --overwrite, -o
        If a design with the given name already exists, remove it and replace 
        it with the new design created by this script.
"""

import os, re, shutil, subprocess, glob

def ensure_path_exists(path):
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise ValueError("'{0}' does not exist.".format(path))
    return path

class PythonPath:
    prompt = "Path to python binary where PyRosetta is installed: "
    description = """\
Python path: Path to your python binary. Make sure PyRosetta package is
installed. You can skip this step in the future by adding a ROSEASY_PYTHON
variable to your environment.
    """

    if 'ROSEASY_PYTHON' in os.environ:
        print('ROSEASY_PYTHON found in environment: {}'.format(
            os.environ.get('ROSEASY_PYTHON')))
        setting_env_var = os.environ.get('ROSEASY_PYTHON')

    @staticmethod
    def install(workspace, python_path):
        python_path = ensure_path_exists(python_path)

        os.symlink(python_path, workspace.python_path)

class RIFDock:
    prompt = "Path to python binary where RIFDock is installed: "
    description = """\
Python path: Path to folder where RIFDock package is
installed. You can skip this step in the future by adding a RIFDOCK 
variable to your environment.
    """

    if 'RIFDOCK' in os.environ:
        print("'rifdock' folder found in environment: {}".format(
            os.environ.get('RIFDOCK')))
        setting_env_var = os.environ.get('RIFDOCK')

    @staticmethod
    def install(workspace, rifdock_path):
        rifdock_bin = os.path.join(rifdock_path, 'build', 'apps',
                'rosetta', 'rif_dock_test')
        rifgen_bin = os.path.join(rifdock_path, 'build', 'apps',
                'rosetta', 'rifgen')
        rifdock_bin = ensure_path_exists(rifdock_bin)
        rifgen_bin = ensure_path_exists(rifgen_bin)

        os.symlink(rifdock_bin, workspace.rifdock)
        os.symlink(rifgen_bin, workspace.rifgen)


class RosettaDir:
    prompt = "Path to rosetta: "
    description = """\
    Rosetta checkout: Path to the main directory of a Rosetta source code checkout.  
    This is the directory called 'main' in a normal rosetta checkout.  Rosetta is 
    used both locally and on the cluster, but the path you specify here probably 
    won't apply to both machines.  You can manually correct the path by changing 
    the symlink called 'rosetta' in the workspace directory."""

    if 'ROSETTA_DIR' in os.environ:
        setting_env_var = os.environ.get('ROSETTA_DIR')

    @staticmethod
    def install(workspace, rosetta_dir):
        rosetta_dir = ensure_path_exists(rosetta_dir)
        rosetta_subdirs = [
                os.path.join(rosetta_dir, 'database'),
                os.path.join(rosetta_dir, 'tests'),
                os.path.join(rosetta_dir, 'source'),
                os.path.join(rosetta_dir, 'source', 'bin'),
        ]
        rosetta_subdirs_exist = list(map(os.path.exists, rosetta_subdirs))

        if not all(rosetta_subdirs_exist):
            message = [
                    "'{0}' does not appear to be the main rosetta directory.".format(rosetta_dir),
                    "The following subdirectories are missing:"
            ]
            for path in rosetta_subdirs:
                if not os.path.exists(path):
                    message.append('    ' + path)
            raise ValueError('\n'.join(message))

        os.symlink(rosetta_dir, workspace.rosetta_dir)


class InputPdb:
    prompt = "Path to the input PDB file: "
    description = """\
Input PDB file: A structure containing the functional groups to be positioned.  
This file should already be parse-able by rosetta, which often means it must be 
stripped of waters and extraneous ligands."""

    @staticmethod
    def install(workspace, pdb_path):
        pdbs = pdb_path.split(' ')
        globlist = []
        for pdb in pdbs:
            globlist.extend(glob.glob(pdb))
        os.makedirs(workspace.target_dir, exist_ok=True)
        for f in globlist:
            f = ensure_path_exists(f)
            destination = os.path.join(workspace.target_dir,
                    os.path.basename(f))
            if f.endswith('.pdb.gz'):
                shutil.copyfile(f, destination)
            elif pdb_path.endswith('.pdb'):
                destination += '.gz'
                subprocess.call('gzip -c {0} > {1}'.format(
                        f, destination), shell=True)
            else:
                raise ValueError("'{0}' is not a PDB file.".format(pdb_path))

class Database:
    description="""\
Symlinks to the default database. Users can run the bin command to
create a new database directory in project_params. Database should be 
downloaded from <url> and placed in the folder's install directory.
(Will make this automatic during setup.py)"""
    
    @staticmethod
    def install(workspace):
        default_dbpath = os.path.join(
                os.path.realpath(__file__),
                '..', 'database'
                )
        os.symlink(default_dbpath, os.path.join(
            workspace.standard_params_dir, 'database'
            ))

class Helices:
    prompt = None
    description = """\
Installs helices to be docked to the target PDBs."""

    @staticmethod
    def install(workspace):
        helicepaths = glob.glob(
                os.path.join(
                    os.path.dirname(__file__), '..', 'helices', '*.pdb'
                    )
                )
        os.makedirs(workspace.helix_dir, exist_ok=True)
        for f in helicepaths:
            f = ensure_path_exists(f)
            workspace_path = os.path.join(workspace.helix_dir,
                    os.path.basename(f))
            shutil.copyfile(f, workspace_path)


class DefaultScripts:
    prompt = None
    Description = """\
Installing default scripts."""
    
    @staticmethod
    def install(workspace):
        script_dir = os.path.join(os.path.dirname(__file__), '..',
                'standard_params')
        python = glob.glob(script_dir + '/*.py')
        yaml = glob.glob(script_dir + '/*.yml')
        wts = glob.glob(script_dir + '/*.wts')
        sho = glob.glob(script_dir + '/*.sho')
        pickle = glob.glob(script_dir + '/*.pkl')
        for script in python + yaml + wts + pickle:
            script_path = os.path.join(script_dir, script)
            workspace_path = os.path.join(workspace.standard_params_dir,
                    os.path.basename(script))
            shutil.copyfile(script_path, workspace_path)
        for script in sho:
            script_path = os.path.join(script_dir, script)
            workspace_path = os.path.join(workspace.root_dir,
                    os.path.basename(script))
            shutil.copyfile(script_path, workspace_path)


class RsyncUrl:
    prompt = "Path to project on remote host: "
    description = """\
Rsync URL: An ssh-style path to the directory that contains (i.e. is one level 
above) the remote workspace.  This workspace must have the same name as the 
remote one.  For example, to link to "~/path/to/my_design" on chef, name this 
workspace "my_design" and set its rsync URL to "chef:path/to"."""

    @staticmethod
    def install(workspace, rsync_url):
        with open(workspace.rsync_url_path, 'w') as file:
            file.write(rsync_url.strip() + '\n')



from klab import scripting
import docopt
from helix import workspace as ws

@scripting.catch_and_print_errors()
def main():
    arguments = docopt.docopt(__doc__)
    workspace = ws.Workspace(arguments['<workspace>'])

    # Make a new workspace directory.

    if workspace.incompatible_with_fragments_script:
        scripting.print_error_and_die("""\
Illegal character(s) found in workspace path:

  {}

The full path to a workspace must contain only characters that are alphanumeric
or '.' or '_'.  The reason for this ridiculous rule is the fragment generation
script, which will silently fail if the full path to its input file contains 
any characters but those.""", workspace.abs_root_dir)

    if workspace.exists():
        if arguments['--overwrite']:
            shutil.rmtree(workspace.root_dir)
        else:
            scripting.print_error_and_die("""\
Design '{0}' already exists.  Use '-o' to overwrite.""", workspace.root_dir)

    workspace.make_dirs()

    # Decide which settings to ask for.

    if arguments['--remote']:
        installers = (
                # RosettaDir,
                RsyncUrl,
                PythonPath,
        )
    else:
        installers = (
                # RosettaDir,
                InputPdb,
                PythonPath,
                DefaultScripts,
                Helices,
                RIFDock,
                Database,
                # LoopsFile,
                # Resfile,
                # ParamsFile,
        )

    # Get the necessary settings from the user and use them to fill in the 
    # workspace.

    print("Please provide the following pieces of information:")
    print()

    scripting.use_path_completion()

    for installer in installers:

        # If the installer doesn't have a prompt, just install it without 
        # asking any questions.

        if installer.prompt is None:
            installer.install(workspace)
            continue

        # Check if an environment variable is defined for the installer
        if hasattr(installer, 'setting_env_var'):
            installer.install(workspace, installer.setting_env_var)
            continue

        # Otherwise, print a description of the setting being installed and 
        # prompt the user for a value.


        print(installer.description)
        print()

        while True:
            try:
                setting = input(installer.prompt)
                installer.install(workspace, setting)
            except (ValueError, IOError) as problem:
                print(problem)
                continue
            except (KeyboardInterrupt, EOFError):
                shutil.rmtree(workspace.root_dir)
                scripting.print_error_and_die("\nReceived exit command, no workspace created.")
            else:
                break

        print()

    # If we made a link to a remote workspace, immediately try to synchronize 
    # with it.  Rsync will say whether or not it succeeded.  Otherwise just 
    # print a success message.

    if arguments['--remote']:
        ws.fetch_data(workspace.root_dir)
    else:
        print("Setup successful for design '{0}'.".format(workspace.root_dir))

if __name__=='__main__':
    main()
