#!/usr/bin/env python 

## \file mesh_adaptation.py
#  \brief python script for mesh adaptation with SU2-MMG interface
#  \author Alberto Perlini, Riccardo Guglielmi
#  \version 8.0.1 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import sys
import SU2
from optparse import OptionParser

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main(): 

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=0,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-b", "--binary", dest="meditformat", default=1,
                      help="Medit mesh format (binary = 1, ascii = 0)", metavar="MEDITFORMAT")
    parser.add_option("-r", "--runcfd", dest="runcfd", default=1,
                      help="run CFD simulation after adaptation", metavar="RUNCFD")

    (options, args)=parser.parse_args()

    options.partitions = int( options.partitions )  
    options.meditformat = int( options.meditformat )
    options.runcfd = int( options.runcfd )

    sys.stdout.write(
        "\n-------------------------------------------------------------------------\n"
    )
    sys.stdout.write(
        "|    ___ _   _ ___                                                      |\n"
    )
    sys.stdout.write(
        '|   / __| | | |_  )   Release 8.0.1 "Harrier"                           |\n'
    )
    sys.stdout.write(
        "|   \\__ \\ |_| |/ /                                                      |\n"
    )
    sys.stdout.write(
        "|   |___/\\___//___|   Anisotropic Mesh Adaptation Script                |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "-------------------------------------------------------------------------\n"
    )
    sys.stdout.write(
        "| SU2 Project Website: https://su2code.github.io                        |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| The SU2 Project is maintained by the SU2 Foundation                   |\n"
    )
    sys.stdout.write(
        "| (http://su2foundation.org)                                            |\n"
    )
    sys.stdout.write(
        "-------------------------------------------------------------------------\n"
    )
    sys.stdout.write(
        "| Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)                |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| SU2 is free software; you can redistribute it and/or                  |\n"
    )
    sys.stdout.write(
        "| modify it under the terms of the GNU Lesser General Public            |\n"
    )
    sys.stdout.write(
        "| License as published by the Free Software Foundation; either          |\n"
    )
    sys.stdout.write(
        "| version 2.1 of the License, or (at your option) any later version.    |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| SU2 is distributed in the hope that it will be useful,                |\n"
    )
    sys.stdout.write(
        "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n"
    )
    sys.stdout.write(
        "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n"
    )
    sys.stdout.write(
        "| Lesser General Public License for more details.                       |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| You should have received a copy of the GNU Lesser General Public      |\n"
    )
    sys.stdout.write(
        "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n"
    )
    sys.stdout.write(
        "-------------------------------------------------------------------------\n"
    )
    
    # Run Mesh Adaptation
    mesh_adaptation ( options.filename,
                      options.partitions, 
                      options.meditformat, 
                      options.runcfd )

#: def main()


# -------------------------------------------------------------------
#  Mesh Adaptation Function
# -------------------------------------------------------------------

def mesh_adaptation( filename,
                     partitions = 0,
                     meditformat = 1,
                     runCFD = 1):
    
    if not filename:
        sys.stderr.write("  ## ERROR : a .cfg file must be provided.\n");
        sys.exit(1)
    
    # Set the name of the configuration file
    config_name = filename
    
    # Read the specified configuration file
    config = SU2.io.Config(config_name)
    
    # Set the number of partitions for parallel computations
    config.NUMBER_PART = partitions
    
    # Call CFD to generate a solution
    SU2.adap.mmg(config, meditformat, runCFD)
    
#: def mesh_adaptation()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
