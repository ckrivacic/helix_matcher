#!/usr/bin/env python2

"""\
HELIX (Helical ELements Interface eXplorer) is a new approach to
protein-protein interface design. For a given target, HELIX attempts to
(a) find positions on the target protein where helices may be likely to
bind using RIFDock, then (b) search the PDB or some other custom set of
proteins, such as a design library, for scaffolds that would enable
those helical positions.

Usage:
    helix <command> [<args>...]
    helix --version
    helix --help

Arguments:
    <command>
        The name of the command you want to run.  You only need to specify 
        enough of the name to be unique.

{command_table}

    <args>...
        The necessary arguments depend on the command being run.  For more 
        information, pass the '--help' flag to the command you want to run.

Options:
    -v, --version
        Display the version of PIP that's installed.

    -h, --help
        Display this help message.

HELIX's design pipeline has the following steps:

1. Setup a workspace to help deal with files created by HELIX.
    $ helix setup_workspace ...

2. Prepare and run RIFDock

    $ helix prep_rifgen ...

3. Match docked helices to potential scaffolds.

    $ helix match ...

4. Analyze...
(finish documentation later)

"""

import sys, re
from klab import scripting
import docopt
from . import __version__

def make_command_table(entry_points):
    """
    Return a nicely formatted table of all the PIP commands installed on this 
    system to incorporate into the help text.  The table will have two columns.  
    The first will list the commands that comprise the main pipeline and the 
    second will list all the other miscellaneous helper functions.
    """
    import itertools

    # Split every command installed on the system into two categories: those 
    # that are part of the main pipeline and those that are just utilities or 
    # helpers.  Pipeline scripts start with numbers, helper scripts don't.

    pipeline_commands = []
    helper_commands = []

    for command in sorted(entry_points):
        if re.match('\d+_', command):
            pipeline_commands.append(command)
        else:
            helper_commands.append(command)

    # Figure out what the longest pipeline command is, so we know how much to 
    # indent the helper commands.

    longest_pipeline_command = 0
    for command in pipeline_commands:
        longest_pipeline_command = max(len(command), longest_pipeline_command)

    # Make the table.

    rows = []
    columns = itertools.zip_longest(
            pipeline_commands, helper_commands, fillvalue='')

    for commands in columns:
        #row = '        {0[0]:{1}}   {0[1]}'.format(
        row = '        {0[0]}   {0[1]}'.format(
                commands, longest_pipeline_command)
        rows.append(row)

    return '\n'.join(rows)

def did_you_mean(unknown_command, entry_points):
    """
    Return the command with the name most similar to what the user typed.  This 
    is used to suggest a correct command when the user types an illegal 
    command.
    """
    from difflib import SequenceMatcher
    similarity = lambda x: SequenceMatcher(None, x, unknown_command).ratio()
    did_you_mean = sorted(entry_points, key=similarity, reverse=True)
    return did_you_mean[0]

@scripting.catch_and_print_errors()
def main():
    from pkg_resources import iter_entry_points, DistributionNotFound

    # Load every PIP command installed on this system.  This is cool because by 
    # using ``pkg_resources``, other packages can add commands to PIP!
    entry_points = {}
    for entry_point in iter_entry_points(group='helix.commands'):
        entry_points[entry_point.name] = entry_point

    # Read the command the user typed on the command line.
    command_table = make_command_table(entry_points)
    arguments = docopt.docopt(
            __doc__.format(**locals()),
            version=__version__,
            options_first=True,
    )
    command_name = arguments['<command>']

    # Find all the commands that match what the user typed.
    matching_entry_points = [
            name for name in entry_points
            if name.startswith(command_name)]

    # If no commands match, print out an error and suggest a command the user 
    # might have been trying to type.
    if len(matching_entry_points) == 0:
        scripting.print_error_and_die("""\
Unknown command '{0}'.  Did you mean:

    $ helix {1} {2}

""", command_name, did_you_mean(command_name, entry_points), ' '.join(arguments['<args>']))

    # If two or more commands match, print all the ambiguous commands and tell 
    # the user to be more specific.
    elif len(matching_entry_points) > 1:
        message = "Command '{0}' is ambiguous.  Did you mean:\n\n"
        for matching_entry_point in matching_entry_points:
            message += "    $ helix {0} {{1}}\n".format(matching_entry_point)
        message += '\n'
        scripting.print_error_and_die(message, command_name, ' '.join(arguments['<args>']))

    # If a unique command was given, make sure all of its dependencies are 
    # installed (because the dependencies for the analysis scripts are not by 
    # default).  If there is a problem, suggest how to fix it.  Otherwise, run 
    # the command.
    else:
        entry_point = entry_points[matching_entry_points[0]]

        try:
            entry_point.require()
        except DistributionNotFound as error:
            scripting.print_error_and_die("""\
The '{0}' command requires the '{1.req}' package.

Installation will be further streamlined in the future. In the meantime,
the helix.yml conda environment should get you most of the packages
(except PyRosetta and klab):

conda env create --file helix.yml

The latest version of klab can be installed via the following commands:

git clone https://github.com/kortemme-lab/klab
cd klab
python setup.py install

Any other dependencies should be available via conda or pip.
""".format(command_name, error))
    #'$ pip install 'helix [analysis]'

        sys.argv = sys.argv[:1] + matching_entry_points + arguments['<args>']
        entry_point.load()()
