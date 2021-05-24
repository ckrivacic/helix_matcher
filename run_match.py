'''
When submitting a parent directory for matching on the cluster, this
script figures out which sub-directory to use and runs the command.
'''
import sys, os, glob
import yaml
from utils import run_command

def main():
    folders = glob.glob(sys.argv[1] + '/*/cluster_representatives/')
    folders = sorted(folders)
    task = int(os.environ['SGE_TASK_ID'])
    folder = os.path.abspath(folders[task - 1])
    python = os.environ['ROSEASY_PYTHON']
    matcher_script = os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)
                ),
            'matcher.py'
            )
    settings = yaml.load(open('settings.yml', 'r'))

    cmd = python, matcher_script, 'match'
    cmd += 'match', folder
    cmd += '-a', settings['match']['-a']
    cmd += '-d',  settings['match']['-d']
    cmd += '--database', settings['match']['--database']
    cmd += '--out', os.path.join(settings['match']['--out'],
            os.path.basename(folder))

    print(cmd)

    run_command(cmd)


if __name__=='__main__':
    main()
