#!/usr/bin/env python3

from sys import argv, stderr
# from genies import FeGenie, LithoGenie, RosGenie, MagnetoGenie, WspGenie, Lucifer, HmmGenie, GasGenie, MnGenie, CircGenie, PolGenie, YfGenie, RiboGenie
from genies_v2 import FeGenie, LithoGenie, OmniGenie, HmmGenie

"""
MagicLamp.py: A script for querying HMMs against provided datasets and processing output.
Installation requirements:
    *python3
    *libraries from the Python standard library: see FeGenie.py and HmmGenie
 """
__author__ = "Arkadiy Garber"
__version__ = "2"
__maintainer__ = "Arkadiy Garber"
__email__ = "ark@midauthorbio.com"

errorMessage = "Options: MagicLamp.py [ FeGenie | LithoGenie | OmniGenie | HmmGenie | help ]\n"

try:
    argv[1]
except IndexError:
    stderr.write(errorMessage)
    exit()

if argv[1] == "FeGenie":
    FeGenie.main()
elif argv[1] == "LithoGenie":
    LithoGenie.main()
elif argv[1] == "OmniGenie":
    OmniGenie.main()
elif argv[1] == "HmmGenie":
    HmmGenie.main()
elif argv[1] == "help":
    stderr.write("\tMagicLamp.py FeGenie: HMM-based identification and categorization of iron genes and iron gene operons in genomes and metagenomes.\n"
                 
                 "\tMagicLamp.py LithoGenie: HMM-based identification and categorization of genes and operons relevant to chemolithoautotrophic metabolisms.\n"
                 
                 "\tMagicLamp.py OmniGenie: HMM-based identification for a given genie.\n"
                 
                 "\tMagicLamp.py HmmGenie: Identification of a user-provided set of HMMs.\n")
    exit()
else:
    stderr.write(errorMessage)
    exit()

