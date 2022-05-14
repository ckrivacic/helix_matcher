import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

def define_command(module, extras=None):
    entry_point = '{0} = helix.commands.{0}:main'.format(module)
    if extras is not None:
        entry_point += ' {0}'.format(extras)
    return entry_point

setuptools.setup(
    name="helix-codykrivacic", # Replace with your own username
    version="0.0.1",
    author="Cody Krivacic",
    author_email="krivacic@berkeley.edu",
    description="Helical elements interface explorer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ckrivacic/helix_matcher",
    packages=setuptools.find_packages(),
    package_data={
        'helix': [
            'standard_params/*.py',
            'standard_params/*.wts',
            'standard_params/*.yml',
            '*.png',
            '*.py',
            'rifdock/*.py',
            'rifdock/*.sh',
            'commands/*.py',
            'matching/*.py',
            'utils/*.py',
            'matching/*.py',
            'matching/*.sh',
            'analysis/*.py',
            'design/*.py',
            ]
        },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'docopt',
        ],
    extras_require={
        'analysis':[
            'numpy',
            'scipy',
            'pandas',
            'matplotlib==3.1.3',
            'nonstdlib',
            'seaborn',
            ]
        },
    entry_points={
        'console_scripts':[
            'helix=helix.main:main',
            ],
        'helix.commands': [
            define_command('fetch_data'),
            define_command('push_data'),
            define_command('combine_dataframes'),
            # define_command('submit'),
            # define_command('plot_funnels'),
            # define_command('violin_plot'),
            # define_command('generate_fragments'),
            define_command('setup_workspace'),
            define_command('analyze_designs'),
            define_command('build_database'),
            define_command('scan_pdb_folder'),
            define_command('bin_database'),
            define_command('view_matches'),
            define_command('prep_patchman'),
            define_command('patchman'),
            define_command('design_patchman'),
            define_command('01_prep_rifdock'),
            define_command('02_rifdock'),
            define_command('03_cluster'),
            define_command('04_match'),
            define_command('05_score_matches'),
            define_command('filter'),
            define_command("export_matches"),
            define_command('docking_benchmark'),
            define_command('plot_filter'),
            define_command('06_design_scaffolds')
            ],
        }
)

download = input("Install default database (22 MB)? (y/N) ")
yes = ['y', 'yes', 'oui', 'si', 'yeppers', 'yessir', 'heckyeah']
no = ['no', 'n', 'non', 'noway', 'fohgetaboudit']
if download.lower() in yes:
    import wget, tarfile
    url = 'https://guybrush.ucsf.edu/HELIX_default_database_2021-07-19.tar.gz'
    this_directory = os.path.dirname(os.path.realpath(__file__))
    output_directory = os.path.join(
            this_directory,
            'helix'
            )
    os.chdir(output_directory)
    fname = url.split('/')[-1]
    output_filename = os.path.join(output_directory,
            fname)
    filename = wget.download(url, out=output_filename)
    tar = tarfile.open(filename, 'r:gz')
    tar.extractall()
    tar.close()
    os.remove(filename)
    
    helix_dataframe_url = 'https://guybrush.ucsf.edu/helixdf_nrpdb_2021-07-22.pkl'
    output_directory = os.path.join(this_directory, 'helix', 'database')
    outfile = wget.download(helix_dataframe_url, output_directory)

elif download.lower() in no:
    pass

else:
    print("Input not recognized; if you wish to download the default "\
    "database, run setup again and enter 'y' when prompted.")
