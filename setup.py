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
            # define_command('submit'),
            # define_command('plot_funnels'),
            # define_command('violin_plot'),
            # define_command('generate_fragments'),
            define_command('setup_workspace'),
            define_command('review_matches'),
            define_command('01_prep_rifdock'),
            define_command('build_database'),
            define_command('02_rifdock'),
            define_command('match_helices'),
            define_command('03_align_rifdock_outputs'),
            define_command('03_cluster'),
            # define_command('pick_designs_to_validate'),
            # define_command('add_residues'),
            # define_command('make_table'),
            # define_command('web_logo'),
            ],
        }
)

download = input("Install default database (22 MB)? (y/N) ")
yes = ['y', 'yes', 'oui', 'si', 'yeppers', 'yessir']
no = ['no', 'n', 'non', 'noway', 'forgetaboutit']
if download.lower() in yes:
    import wget, tarfile
    url = 'https://guybrush.ucsf.edu/HELIX_default_database_2021-07-19.tar.gz'
    output_directory = os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)
                ),
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
elif download.lower() in no:
    pass
else:
    print("Input not recognized; if you wish to download the default "\
    "database, run setup again and enter 'y' when prompted.")
