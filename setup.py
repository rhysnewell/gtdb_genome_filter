import io
from os.path import dirname, join
from setuptools import setup


# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]


setup(
    name='aviary',
    version=get_version("gtdb_genome_filter/__init__.py"),
    url='https://github.com/rhysnewell/gtdb_genome_filter',
    license='BSD-3',
    author='Rhys Newell',
    author_email='rhys.newell@hdr.qut.edu.au',
    description='gtdb genome filter - Choosing out genomes from gtdb taxonomy.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['gtdb_genome_filter'],
    package_data={'': [
            "gtdb_genome_filter/*",
                       ]},
    data_files=[(".", ["README.md", "LICENSE"])],
    include_package_data=True,
    install_requires= [
        "pandas",
        "numpy"
    ],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={
          'console_scripts': [
              'gtdb_filter = gtdb_filter.gtdb_filter:main'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)