#!/usr/bin/env python
###############################################################################
# gtdb_filter.py - Info about gtdb_filter.py
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################

__author__ = "Rhys Newell"
__copyright__ = "Copyright 2020"
__credits__ = ["Rhys Newell"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Rhys Newell"
__email__ = "rhys.newell near hdr.qut.edu.au"
__status__ = "Development"

###############################################################################
# System imports
import sys
import argparse
import logging
import os
import shutil
from datetime import datetime
import subprocess

# Local imports
import pandas as pd
import numpy as np

# Debug
debug = {1: logging.CRITICAL,
         2: logging.ERROR,
         3: logging.WARNING,
         4: logging.INFO,
         5: logging.DEBUG}


###############################################################################
############################### - Exceptions - ################################

class BadTreeFileException(Exception):
    pass


###############################################################################                                                                                                                      [44/1010]
################################ - Functions - ################################


def phelp():
    print(
        """
    aviary
    
    SUBCOMMAND:
    recover
    """
    )


def main():
    ############################ ~ Main Parser ~ ##############################
    main_parser = argparse.ArgumentParser(prog='gtdb_filter',
                                          formatter_class=CustomHelpFormatter,
                                          add_help=False)
    main_parser.add_argument('--version',
                             action='version',
                             version=__version__,
                             help='Show version information.')
    main_parser.add_argument('--verbosity',
                             help='1 = critical, 2 = error, 3 = warning, 4 = info, 5 = debug. Default = 4 (logging)',
                             type=int,
                             default=4)
    main_parser.add_argument('--log',
                             help='Output logging information to file',
                             default=False)
    subparsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    ########################## ~ sub-parser ~ ###########################
    input_options = subparsers.add_parser('filter',
                                          description='Get a filtered list of GTDB genomes',
                                          formatter_class=CustomHelpFormatter,
                                          epilog='''
                                ~ FILTER ~
    How to use recover:

    gtdb_filter filter --input gtdb_db_metadata_file.tsv --completeness 90 --contamination 5 --taxonomy f_Bacillaceae

    ''')

    input_options.add_argument(
        '--input',
        help='GTDB metadata file',
        dest="input",
        required=True,
    )

    input_options.add_argument(
        '--completeness',
        help='Minimum genome completeness',
        dest='completeness',
        default=90
    )

    input_options.add_argument(
        '--contamination',
        help='Maximum genome contamination',
        dest='contamination',
        default=5
    )

    input_options.add_argument(
        '--taxonomy',
        help='If supplied, the output will only contain organisms from within this taxonomic group. e.g. f_Bacillaceae',
        dest='taxonomy',
        default="all"
    )

    input_options.add_argument(
        '--representatives_only',
        help='Return only GTDB representatives',
        dest='representative',
        default=True
    )

    input_options.add_argument(
        '--output',
        help='Output directory, outputs to current directory *DON"T CHANGE, relative paths currently broken for this*',
        dest='output',
        default='filtered_gtdb_output.tsv',
    )


    ###########################################################################
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parsing input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
    else:
        args = main_parser.parse_args()
        time = datetime.now().strftime('%H:%M:%S %d-%m-%Y')

        if args.log:
            if os.path.isfile(args.log):
                raise Exception("File %s exists" % args.log)
            logging.basicConfig(filename=args.log,
                                level=debug[args.verbosity],
                                format='%(asctime)s %(levelname)s: %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p')
        else:
            logging.basicConfig(level=debug[args.verbosity],
                                format='%(asctime)s %(levelname)s: %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p')
        logging.info("Time - %s" % (time))
        logging.info("Command - %s" % ' '.join(sys.argv))

        prefix = args.output
        if not os.path.exists(prefix):
            os.makedirs(prefix)

        filter(args.input,
               float(args.completeness),
               float(args.contamination),
               args.taxonomy,
               str2bool(args.representative),
               args.output)


def filter(input_file, completeness=90, contamination=5, taxonomy="all", representative=True, output="filtered_gtdb_metadata.tsv"):
    metadata = pd.read_csv(input_file, sep="\t")
    col_to_return = ["acccession", "gtdb_genome_representative", "checkm_completeness", "checkm_contamination", "ncbi_genbank_assembly_accession",
                     "ncbi_isolation_source", "coding_density", "gc_percentage", "genome_size", "gtdb_taxonomy",
                     "ncbi_taxonomy"]
    if representative:
        if taxonomy.lower() == "all":
            filtered = metadata[(metadata['checkm_completeness'] >= completeness)
                                & (metadata['checkm_contamination'] >= contamination)
                                & (metadata["gtdb_representative"] == "t")]

        else:
            filtered = metadata[(metadata['checkm_completeness'] >= completeness)
                                & (metadata['checkm_contamination'] >= contamination)
                                & (metadata["gtdb_representative"] == "t")
                                & (metadata.gtdb_taxonomy.str.contains(taxonomy))]
    else:
        if taxonomy.lower() == "all":
            filtered = metadata[(metadata['checkm_completeness'] >= completeness)
                                & (metadata['checkm_contamination'] >= contamination)
                                & (metadata["gtdb_representative"] == "f")]

        else:
            filtered = metadata[(metadata['checkm_completeness'] >= completeness)
                                & (metadata['checkm_contamination'] >= contamination)
                                & (metadata["gtdb_representative"] == "f")
                                & (metadata.gtdb_taxonomy.str.contains(taxonomy))]



    filtered = filtered[col_to_return]
    filtered.to_csv(output, sep="\t")

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


###############################################################################
################################ - Classes - ##################################

class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
                    action.default != [] and \
                    action.default != None \
                    and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL,
                                        argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])


if __name__ == '__main__':
    sys.exit(main())
