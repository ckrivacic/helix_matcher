#!/usr/bin/env python2

"""\
This module defines the Workspace classes that are central to every script.
The role of these classes is to provide paths to all the data files used in any
part of the pipeline and to hide the organization of the directories containing
those files.  The base Workspace class deals with files in the root directory
of a design.  It's subclasses deal with file in the different subdirectories of
the design, each of which is related to a cluster job.
"""

import os, re, glob, json, pickle, sys
from klab import scripting
from pprint import pprint
import yaml
# from roseasy.standard_params import *

class Workspace(object):
    """
    Provide paths to every file used in the design pipeline.

    Each workspace object is responsible for returning paths to files that are
    relevant to a particular stage of the design pipeline.  These files are
    organized hierarchically: files that are relevant to many parts of the
    pipeline are stored in the root design directory while files that are
    relevant to specific stages are stored in subdirectories.  You can think of
    each workspace class as representing a different directory.

    The Workspace class itself represents the root directory, but it is also
    the superclass from which all of the other workspace derive.  The reason
    for this is that the root workspace knows where all the shared parameter
    files are located, and this information is needed in every workspace.

    When modifying or inheriting from this class, keep in mind two things.
    First, workspace objects should do little more than return paths to files.
    There are a few convenience functions that clear directories and things
    like that, but these are the exception rather than the rule.  Second, use
    the @property decorator very liberally to keep the code that uses this API
    succinct and easy to read.
    """

    def __init__(self, root):
        self._root_dir = os.path.abspath(root)

    def job_info_path(self, job_id):
        return os.path.join(self.focus_dir, '{0}.json'.format(job_id))

    @classmethod
    def from_directory(cls, directory):
        # Force subclasses to reimplement this method
        if cls != Workspace:
            raise NotImplementedError

        return Workspace(directory)

    @property
    def incompatible_with_fragments_script(self):
        return re.search('[^a-zA-Z0-9_/.]', self.abs_root_dir)

    @property
    def settings_filename(self):
        return 'settings.yml'

    @property
    def settings_path(self):
        return self.find_path(self.settings_filename)

    @property
    def log_dir(self):
        return os.path.join(self.project_params_dir, 'logs')

    @property
    def settings(self):
        settings_paths =\
                self.find_all_paths(self.settings_filename)
        these_settings = None
        settings = {}
        # Open settings files in order from lowest to highest
        # priority
        for f in reversed(settings_paths):
            with open(f, 'r') as stream:
                try:
                    these_settings = yaml.safe_load(stream)
                except yaml.YAMLError as exec:
                    print(exec)
            # Update settings in reverse priority order
            if these_settings:
                for category in these_settings:
                    # Some settings may be nested dictionaries, others
                    # not. For nested dictionaries, update items
                    # individually.
                    if category not in settings:
                        settings[category] = these_settings[category]
                    if isinstance(these_settings[category], dict):
                        for key in these_settings[category]:
                            settings[category][key] =\
                                    these_settings[category][key]
        return settings

    @property
    def database_path(self):
        return self.find_path('database')

    @property
    def dataframe_path(self):
        dbpath = self.settings['match']['--database']
        picklepath = glob.glob(os.path.join(dbpath, 'helixdf*.pkl'))
        if len(picklepath) > 1:
            print("Multiple helix dataframes found in {0}! Please consolodate them "\
                    "or remove extraneous pkl files beginning with "\
                    "'helixdf'.".format(database))
            sys.exit()
        elif len(picklepath) == 0:
            print("No helix dataframe ('helixdf*.pkl') found in provided "\
                    "database path ({0}). Please runs 'helix scan' command "\
                    "before binning.")
        else:
            helixdf = picklepath[0]
        return helixdf

    def clear_database(self):
        scripting.clear_directory(os.path.join(
            self.project_params_dir, 'database'
            ))

    def is_default_database(self, dbpath):
        '''Checks if a database path is the default one by
        checking whether it is in standard_params or not'''
        return os.path.dirname(os.path.abspath(dbpath)) == \
                os.path.abspath(self.standard_params_dir)

    @property
    def parent_dir(self):
        return os.path.dirname(self.root_dir)

    @property
    def root_dir(self):
        return self._root_dir

    @property
    def abs_root_dir(self):
        return os.path.abspath(self.root_dir)

    def basename(self, target):
        if os.path.isfile(target):
            if target.endswith('.pdb.gz'):
                return os.path.basename(target)[:-len('.pdb.gz')]
            elif target.endswith('.pdb'):
                return os.path.basename(target)[:-len('.pdb')]
            else:
                return os.path.basename(target).split('.')[0]
        elif os.path.isdir(target):
            return os.path.basename(target)

    @property
    def focus_name(self):
        return ''

    @property
    def rifdock_outdir(self):
        return os.path.join(self.root_dir, 'rifdock_outputs')
    
    def target_rifdock_path(self, target):
        dirname = self.basename(target)
        return os.path.join(self.rifdock_outdir, dirname)

    def target_clusters(self, target):
        dirname = self.basename(target)
        return os.path.join(self.target_rifdock_path(target),
                'cluster_representatives')

    def target_path_from_name(self, name):
        return os.path.join(self.target_dir, '{}.pdb.gz'.format(name))

    @property
    def all_rifdock_workspaces(self):
        all_target_paths = []
        for target in self.targets:
            all_target_paths.append(self.target_rifdock_path(target))
        return all_target_paths

    @property
    def match_outdir(self):
        return os.path.join(self.root_dir, 'match_outputs')

    def target_match_path(self, target):
        out_folder = os.path.join(self.match_out_dir,
                self.basename(target))

    @property
    def slurm_custom_jobno(self):
        """
        Return the number associated with the latest slurm command file. 
        Slurm command files should have the format slurm_<jobnumber>, where 
        jobnumber increments by 1 for each job submitted in a particular 
        workspace. This function simply returns the latest <jobnumber>.
        """
        latest = 0
        for f in glob.glob(self.focus_dir + '/slurm_*'):
            f = os.path.basename(f)
            num = int(f.split('_')[1])
            if num > latest:
                latest = num

        return latest

    @property
    def user_email(self):
        '''
        Slurm jobs can accept an email to be submitted with a job. This can
        be provided in either standard_params or project_params.
        '''
        path = self.find_path('email')
        if not os.path.exists(path):
            print('WARNING: Slurm requires an email to submit a job.')
            print('Please provide this via a file named \'email\' in project_params.')
            raise PathNotFound(path)
        with open(path, 'r') as f:
            email = f.read().strip('\n')
        return email

    @property
    def slurm_submit_file(self):
        jobno = self.slurm_custom_jobno
        return os.path.join(self.focus_dir,
                'submit_{}.slurm'.format(jobno))

    @property
    def slurm_cmd_file(self):
        """
        Returns the latest slurm command file.
        """
        jobno = self.slurm_custom_jobno
        return os.path.join(self.focus_dir, 'slurm_{}'.format(jobno))

    @property
    def focus_dir(self):
        """
        The particular directory managed by this class.  This is meant to be
        overridden in subclasses.
        """
        return self.root_dir

    @property
    def standard_params_dir(self):
        return os.path.join(self.root_dir, 'standard_params')

    @property
    def project_params_dir(self):
        return os.path.join(self.root_dir, 'project_params')

    @property
    def io_dirs(self):
        return []

    @property
    def python_path(self):
        return self.find_path('python', self.root_dir)

    @property
    def rosetta_dir(self):
        return self.find_path('rosetta', self.root_dir)

    @property
    def rosetta_scripts_path(self):
        pattern = self.rosetta_subpath('source', 'bin', 'rosetta_scripts*')
        executables = glob.glob(pattern)

        # Sometimes dead symlinks end up in the `bin/` directory, so explicitly
        # ignore those.

        executables = [x for x in executables if os.path.exists(x)]

        # Print a (hopefully) helpful error message if no ``rosetta_scripts``
        # executables are found.

        if len(executables) == 0:
            raise PipelineError("""\
No RosettaScripts executable found.

Expected to find a file matching '{0}'.  Did you forget to compile rosetta?
""".format(pattern))

        # Sort the ``rosetta_scripts`` executables such that those containing
        # the word 'release' end up at the front of the list, those containing
        # 'debug' end up at the back, and shorter file names (which have fewer
        # weird compilation options) end up in front of longer ones.  We'll
        # ultimately pick the first path in the list, so we're doing our best
        # to use a basic release mode executable.

        executables.sort(key=lambda x: len(x))
        executables.sort(key=lambda x: 'debug' in x)
        executables.sort(key=lambda x: 'release' not in x)

        return executables[0]

    @property
    def rosetta_database_path(self):
        return self.rosetta_subpath('database')

    def rosetta_vall_path(self, test_run=False):
        return os.path.join(self.rosetta_database_path, 'sampling',
                'small.vall.gz' if test_run else 'vall.jul19.2011')

    def rosetta_subpath(self, *subpaths):
        return os.path.join(self.rosetta_dir, *subpaths)

    @property
    def input_pdb_path(self):
        return self.find_path('input.pdb.gz')

    @property
    def target_dir(self):
        return os.path.join(self.root_dir, 'targets')

    @property
    def targets(self):
        return sorted(glob.glob(os.path.join(self.target_dir,
            '*.pdb.gz')))

    @property
    def helix_dir(self):
        return os.path.join(self.root_dir, 'helices')

    @property
    def scaffolds(self):
        return sorted(glob.glob(os.path.join(self.helix_dir,
            '*.pdb')) + glob.glob(os.path.join(self.helix_dir,
                '*.pdb.gz')))

    @property
    def flags_path(self):
        return self.find_path('flags')

    @property
    def rsync_url_path(self):
        return self.find_path('rsync_url')

    @property
    def rsync_url(self):
        if not os.path.exists(self.rsync_url_path):
            raise UnspecifiedRemoteHost()
        with open(self.rsync_url_path) as file:
            return file.read().strip()

    @property
    def rsync_recursive_flag(self):
        return False

    @property
    def rsync_include_patterns(self):
        return []

    @property
    def rsync_exclude_patterns(self):
        return ['rosetta', 'rsync_url']

    @property
    def preferred_install_dir(self):
        return os.path.join(self.project_params_dir, self.focus_name)

    @property
    def find_path_dirs(self):
        # The search path has to be a little different for the root directory, 
        # otherwise you end up with some weird behavior dues to the focus 
        # directory (which is the first place to look for files) being the same 
        # as the root directory (which is the last place to look for files).

        if self.focus_dir == self.root_dir:
            return [
                    self.root_dir,
                    self.project_params_dir,
                    self.standard_params_dir,
            ]
        else:
            return [
                    self.focus_dir,
                    os.path.join(self.root_dir, self.focus_name),
                    os.path.join(self.project_params_dir, self.focus_name),
                    self.root_dir,
                    self.project_params_dir,
                    os.path.join(self.standard_params_dir, self.focus_name),
                    self.standard_params_dir,
            ]

    def find_path(self, basename, install_dir=None):
        """
        Look in a few places for a file with the given name.  If a custom
        version of the file is found in the directory being managed by
        this workspace, return it.  Otherwise look in the custom and default 
        input directories in the root directory, and then finally in the root 
        directory itself.

        This function makes it easy to provide custom parameters to any stage
        to the design pipeline.  Just place the file with the custom parameters
        in a directory associated with that stage.
        """

        # Look for the file we were asked for.
        for dir in self.find_path_dirs:
            path = os.path.join(dir, basename)
            if os.path.exists(path):
                return path

        # If we didn't find the file, return the path to where we'd like it to 
        # be installed.
        return os.path.join(install_dir or self.preferred_install_dir, basename)

    def find_all_paths(self, basename):
        """
        Looks in a few places for any files with a given name or
        pattern and returns them as a list.
        """

        # Look for the file in standard folders
        hits = []
        for directory in self.find_path_dirs:
            paths = glob.glob(os.path.join(directory, basename))
            hits.extend([os.path.abspath(path) for path in paths])

        return hits

    def check_paths(self):
        required_paths = []
                # self.input_pdb_path]
        '''
        self.loops_path,
        self.resfile_path,
        self.restraints_path,
        self.build_script_path,
        self.design_script_path,
        self.validate_script_path,
        self.flags_path,
        '''
        for path in required_paths:
            if not os.path.exists(path):
                raise PathNotFound(path)
    
    def make_dirs(self):
        scripting.mkdir(self.focus_dir)
        scripting.mkdir(self.standard_params_dir)
        scripting.mkdir(self.preferred_install_dir)

        pickle_path = os.path.join(self.focus_dir, 'workspace.pkl')
        print(self.__class__)
        with open(pickle_path, 'wb') as file:
            pickle.dump(self.__class__, file)

    def cd(self, *subpaths):
        """
        Change the current working directory and update all the paths in the
        workspace.  This is useful for commands that have to be run from a
        certain directory.
        """
        target = os.path.join(*subpaths)
        os.chdir(target)

    def cd_to_root(self):
        self.cd(self.root_dir)

    def exists(self):
        return os.path.exists(self.focus_dir)

    # def patch_folder(self, patch):
        # return os.path.join(self.root_dir, 'rifdock',
                # 'patch_{}'.format(self.patchno))

    def rifgen_flags_template(self):
        return self.find_path('rifgen_flags')

    @property
    def rifdock(self):
        return os.path.join(self.root_dir, 'rifdock')

    @property
    def rifgen(self):
        return os.path.join(self.root_dir, 'rifgen')

    @property
    def all_job_info_paths(self):
        return glob.glob(os.path.join(self.focus_dir, '*.json'))

    @property
    def all_job_info(self):
        from . import big_jobs
        return [big_jobs.read_job_info(x) for x in self.all_job_info_paths]


class RIFWorkspace(Workspace):
    '''
    Class for handling paths needed to run RIFDock on a given target
    '''
    def __init__(self, root, target_path):
        Workspace.__init__(self, root)
        self._initial_target_path = target_path

    @staticmethod
    def from_directory(directory):
        root = os.path.join(directory, '..', '..')
        target_name = os.path.basename(directory)
        target_path = os.path.join(root, 'targets',
                '{}.pdb.gz'.format(target_name))
        return RIFWorkspace(root, target_path)

    @property
    def target_path(self):
        return os.path.join(self.focus_dir, 'target.pdb')

    @property
    def initial_target_path(self):
        return self._initial_target_path

    @property
    def focus_name(self):
        return self.basename(self.initial_target_path)

    @property
    def focus_dir(self):
        return os.path.join(self.rifdock_outdir, self.focus_name)

    @property
    def scaffold_prefix(self):
        return 'scaffold_'

    def scaffold_folder(self, scaffold_basename):
        return self.scaffold_prefix + scaffold_basename

    @property
    def patches(self):
        return sorted(glob.glob(os.path.join(
            self.focus_dir,
            'patch_*', '{}*'.format(self.scaffold_prefix))))

    def scaffold_dir(self, patch, scaffold):
        return os.path.join(self.focus_dir, 'patch_{}'.format(patch),
                self.scaffold_prefix + scaffold)

    def active_patch(self, job_info):
        return self.patches[job_info['task_id'] % len(self.patches)]

    def scaffold_name(self, scaffold_subfolder):
        '''Get the scaffold name from the scaffold subfolder.
        Useful when running RIFDock, where each task gets a folder and
        must figure out what scaffold to use.'''
        basename = os.path.basename(scaffold_subfolder)
        folder_fields = basename.split('_')
        scaffold_name = '_'.join(folder_fields[1:])
        return scaffold_name

    def scaffold_from_dir(self, scaffold_subfolder):
        '''
        Get the scaffold pdb file from the scaffold subfolder.
        '''
        pdb = os.path.abspath(os.path.join(self.helix_dir,
                self.scaffold_name(scaffold_subfolder) + '.pdb'))
        if os.path.exists(pdb):
            return pdb
        elif os.path.exists(pdb + '.gz'):
            return pdb + '.gz'
        else:
            e = "WARNING! Could not find scaffold with the name {} or {} "\
            "in the 'helices' folder.".format(
                    pdb, pdb + '.gz')
            raise PipelineError(e)

    def docking_directory(self, job_info):
        return os.path.join(self.active_patch(job_info), 'docked_full')

    @property
    def docking_directories(self):
        dirs = []
        for d in self.patches:
            dirs.append(os.path.join(d, 'docked_full'))
        return dirs

    def make_dirs(self):
        scripting.mkdir(self.focus_dir)
        scripting.mkdir(self.log_dir)
        pickle_path = os.path.join(self.focus_dir, 'workspace.pkl')
        with open(pickle_path, 'wb') as file:
            pickle.dump(self.__class__, file)

    def nstruct(self):
        return len(self.patches)

    @property
    def log_dir(self):
        return os.path.join(self.focus_dir, 'logs')

    def clear_outputs(self):
        for docking_dir in self.docking_directories:
            scripting.clear_directory(docking_dir)
        scripting.clear_directory(self.log_dir)
        for f in self.all_job_info_paths:
            os.remove(f)

    def clear_cluster_outputs(self):
        import shutil
        cluster_logs = glob.glob(
                os.path.join(
                    self.log_dir,
                    'cluster*'))
        for log in cluster_logs:
            os.remove(log)
        clustered = glob.glob(
                os.path.join(
                    self.target_clusters(self.initial_target_path), '*'
                    )
                )
        for clst in clustered:
            shutil.rmtree(clst)

    @property
    def unclaimed_inputs(self):
        inputs = set(self.patches)
        for params in self.all_job_info:
            inputs -= set(params['inputs'])
        return sorted(inputs)


class MatchWorkspace(Workspace):
    '''
    Class for handling paths needed to run the matcher on a given target
    '''
    def __init__(self, root, target_path):
        Workspace.__init__(self, root)
        self._initial_target_path = target_path

    @staticmethod
    def from_directory(directory):
        root = os.path.join(directory, '..', '..')
        target_name = os.path.basename(directory)
        target_path = os.path.join(root, 'targets',
                '{}.pdb.gz'.format(target_name))
        return MatchWorkspace(root, target_path)

    @property
    def target_path(self):
        return self.initial_target_path

    @property
    def initial_target_path(self):
        return self._initial_target_path

    @property
    def focus_name(self):
        return self.basename(self.initial_target_path)

    @property
    def focus_dir(self):
        return os.path.join(self.match_outdir, self.focus_name)

    @property
    def log_dir(self):
        return os.path.join(self.focus_dir, 'logs')

    @property
    def output_dir(self):
        return os.path.join(self.focus_dir, 'outputs')

    def make_dirs(self):
        scripting.mkdir(self.focus_dir)
        scripting.mkdir(self.log_dir)
        scripting.mkdir(self.output_dir)
        pickle_path = os.path.join(self.focus_dir, 'workspace.pkl')
        with open(pickle_path, 'wb') as file:
            pickle.dump(self.__class__, file)

    def clear_outputs(self):
        scripting.clear_directory(self.log_dir)
        scripting.clear_directory(self.output_dir)
        if os.path.exists(self.all_scaffold_dataframe):
            os.remove(self.all_scaffold_dataframe)
        for scaffold in self.all_scaffold_clusters:
            if os.path.exists(self.scaffold_dataframe(scaffold)):
                os.remove(self.scaffold_dataframe(scaffold))
        for f in self.all_job_info_paths:
            os.remove(f)

    @property
    def unclaimed_inputs(self):
        '''Unclaimed inputs = targets * dataframes'''
        return

    @property
    def rifdock_workspace(self):
        return os.path.join(self.rifdock_outdir, self.focus_name)

    @property
    def cluster_outputs(self):
        return os.path.join(self.rifdock_workspace,
                'cluster_representatives')

    @property
    def all_scaffold_clusters(self):
        return sorted(glob.glob(os.path.join(
            self.cluster_outputs, '*/'
            )))

    @property
    def all_scaffold_dataframe(self):
        return os.path.join(self.cluster_outputs, 'query_helices.pkl')

    def scaffold_clusters(self, scaffold):
        '''Returns the folder for clusters of a given scaffold. Should
        be able to take either a directory or pdb file as the scaffold
        argument.'''
        scaffold_dir = os.path.join(self.cluster_outputs,
                self.basename(scaffold))
        return scaffold_dir

    def scaffold_dataframe(self, scaffold):
        return os.path.join(self.scaffold_clusters(scaffold),
                'query_helices.pkl')



def big_job_dir():
    return os.path.join(os.path.dirname(__file__), '..', 'big_jobs')

def big_job_path(basename):
    return os.path.join(big_job_dir(), basename)

def predecessor_from_dir(workspace, step):
    """
    Get the workspace that precedes the current step.
    """
    searchstr = '{:02}*'.format(step - 1)
    results = []
    for folder in glob.glob(os.path.join(workspace.root_dir,
        searchstr)):
        if os.path.isdir(folder):
            results.append(folder)
    if len(results) > 1:
        e = 'WARNING! More than one predecessor workspace found. '\
                'Check filenames for conflicts.'
        raise PipelineError(e)
    elif len(results) == 0:
        e = 'WARNING! No predecessor workspaces found. '\
                'Check filenames for conflicts.'
        raise PipelineError(e)
    else:
        return results[0]

def workspace_from_dir(directory, recurse=True):
    """
    Construct a workspace object from a directory name.  If recurse=True, this
    function will search down the directory tree and return the first workspace
    it finds.  If recurse=False, an exception will be raised if the given
    directory is not a workspace.  Workspace identification requires a file
    called 'workspace.pkl' to be present in each workspace directory, which can
    unfortunately be a little fragile.
    """
    directory = os.path.abspath(directory)
    pickle_path = os.path.join(directory, 'workspace.pkl')

    # Make sure the given directory contains a 'workspace' file.  This file is
    # needed to instantiate the right kind of workspace.

    if not os.path.exists(pickle_path):
        if recurse:
            parent_dir = os.path.dirname(directory)
            print(parent_dir)

            # Keep looking for a workspace as long as we haven't hit the root
            # of the file system.  If an exception is raised, that means no
            # workspace was found.  Catch and re-raise the exception so that
            # the name of the directory reported in the exception is meaningful
            # to the user.

            try:
                return workspace_from_dir(parent_dir, parent_dir != '/')
            except WorkspaceNotFound:
                raise WorkspaceNotFound(directory)
        else:
            raise WorkspaceNotFound(directory)

    # Load the 'workspace' file and create a workspace.

    with open(pickle_path, 'rb') as file:
        workspace_class = pickle.load(file)

    return workspace_class.from_directory(directory)

def root_from_dir(directory, recurse=True):
    """
    Similar to workspace_from_dir, but this returns the root directory
    of a workspace rather than a workspace object. 
    """
    # This does not check whether this is a Workspace object and could
    # therefore be prone to errors if we aren't careful.

    directory = os.path.abspath(directory)
    pickle_path = os.path.join(directory, 'workspace.pkl')

    # Make sure the given directory contains a 'workspace' file.  This file is
    # needed to instantiate the right kind of workspace.

    if not os.path.exists(pickle_path):
        if recurse:
            parent_dir = os.path.dirname(directory)

            # Keep looking for a workspace as long as we haven't hit the root
            # of the file system.  If an exception is raised, that means no
            # workspace was found.  Catch and re-raise the exception so that
            # the name of the directory reported in the exception is meaningful
            # to the user.

            try:
                return root_from_dir(parent_dir, parent_dir != '/')
            except WorkspaceNotFound:
                raise WorkspaceNotFound(directory)
        else:
            raise WorkspaceNotFound(directory)

    # Return the directory in which the pkl file was found.

    return pickle_path[:-len('workspace.pkl')]

def load_loops(directory, loops_path=None):
    """
    Return a list of tuples indicating the start and end points of the loops
    that were sampled in the given directory.
    """

    if loops_path is None:
        workspace = workspace_from_dir(directory)
        loops_path = workspace.loops_path

    from klab.rosetta.input_files import LoopsFile
    loops_parser = LoopsFile.from_filepath(loops_path)

    # We have to account for some weird indexing behavior in the loops file
    # parser that I don't really understand.  It seems to shrink the loop by
    # one residue on each side.  At first I thought it might be trying to
    # convert the indices to python indexing, but on second thought I have no
    # idea what it's trying to do.

    return [(x-1, y+1) for x, y in loops_parser.get_distinct_segments()]

def load_resfile(directory, resfile_path=None):
    """
    Return a list of tuples indicating the start and end points of the loops
    that were sampled in the given directory.
    """

    if resfile_path is None:
        workspace = workspace_from_dir(directory)
        resfile_path = workspace.resfile_path

    from klab.rosetta.input_files import Resfile
    return Resfile(resfile_path)

def fetch_data(directory, remote_url=None, recursive=True, include_logs=False, dry_run=False):
    import os, subprocess

    workspace = workspace_from_dir(directory)

    # Try to figure out the remote URL from the given directory, if a
    # particular URL wasn't given.

    if remote_url is None:
        remote_url = workspace.rsync_url

    # Make sure the given directory is actually a directory.  (It's ok if it
    # doesn't exist; rsync will create it.)

    if os.path.exists(directory) and not os.path.isdir(directory):
        print("Skipping {}: not a directory.".format(directory))
        return

    # Compose an rsync command to copy the files in question.  Then either run
    # or print that command, depending on what the user asked for.

    rsync_command = [
            'rsync', '-av',
    ] +   (['--no-recursive'] if not recursive else []) + [
            '--exclude', 'rosetta',
            '--exclude', 'rsync_url',
            '--exclude', 'fragments',
            '--exclude', 'core.*',
            '--exclude', 'sequence_profile*',
            '--exclude', 'python',
            '--exclude', 'standard_params/database',
    ]
    if not include_logs:
        rsync_command += [
                '--exclude', 'logs',
                '--exclude', '*.sc',
        ]

    # This code is trying to combine the remote URL with a directory path.
    # Originally I was just using os.path.join() to do this, but that caused a
    # bug when the URL was something like "chef:".  This is supposed to specify
    # a path relative to the user's home directory, but os.path.join() adds a
    # slash and turns the path into an absolute path.

    sep = '' if remote_url.endswith(':') else '/'
    remote_dir = os.path.normpath(
            remote_url + sep + os.path.relpath(directory, workspace.parent_dir))
    rsync_command += [
            remote_dir + '/',
            directory,
    ]

    if dry_run:
        print(' '.join(rsync_command))
    else:
        subprocess.call(rsync_command)

def fetch_and_cache_data(directory, remote_url=None, recursive=True, include_logs=False):
    from . import structures
    fetch_data(directory, remote_url, recursive, include_logs)

    # Don't try to cache anything if nothing has been downloaded yet.
    if glob.glob(os.path.join(directory, '*.pdb*')):
        structures.load(directory)

def push_data(directory, remote_url=None, recursive=True, dry_run=False):
    import os, subprocess

    workspace = workspace_from_dir(directory)

    if remote_url is None:
        remote_url = workspace.rsync_url

    remote_dir = os.path.normpath(os.path.join(
            remote_url, os.path.relpath(directory, workspace.parent_dir)))

    rsync_command = [
            'rsync', '-av',
    ] +   (['--no-recursive'] if not recursive else []) + [
            '--exclude', 'rosetta',
            '--exclude', 'rsync_url',
            '--exclude', 'logs',
            '--exclude', 'python',
            directory + '/', remote_dir,
    ]

    if dry_run:
        print(' '.join(rsync_command))
    else:
        subprocess.call(rsync_command)


class PipelineError (IOError):

    def __init__(self, message):
        super(PipelineError, self).__init__(message)
        self.no_stack_trace = True


class PathNotFound (PipelineError):

    def __init__(self, path, *directories):
        if len(directories) == 0:
            message = "'{0}' not found.".format(path)

        elif len(directories) == 1:
            path = os.path.join(directories[0], path)
            message = "'{0}' not found.".format(path)

        else:
            message = "'{0}' not found.  Looked in:".format(path)
            for directory in directories:
                message += "\n    " + directory

        PipelineError.__init__(self, message)


class RosettaNotFound (PipelineError):

    def __init__(self, workspace):
        PipelineError.__init__(self, """\
No rosetta checkout found in '{0.root_dir}'.
Use the following command to manually create a symlink to a rosetta checkout:

$ ln -s /path/to/rosetta/checkout {0.rosetta_dir}""")


class WorkspaceNotFound (PipelineError):

    def __init__(self, root):
        message = "'{0}' is not a workspace.".format(root)
        PipelineError.__init__(self, message)


class UnspecifiedRemoteHost (PipelineError):

    def __init__(self):
        PipelineError.__init__(self, "No remote host specified.")
