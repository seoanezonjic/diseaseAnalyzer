#! /usr/bin/env python

import argparse
import os
import re

##############
## METHODS
##############

def  

##############
## OPTPARSE
##############

parser = argparse.ArgumentParser(description="")
parser.add_argument("-F", "--filter_value", type=float, help="Filter value", required=True, default=0.8)
parser.add_argument("-i", "--input_file", metavar="PATH", help="Input file", required=True)
parser.add_argument("-m", "--min_groups", type=int, help="Min diseases per gene", required=True, default=0)
parser.add_argument("-o", "--output_file", metavar="PATH", help="Output file", default="output.txt")
parser.add_argument("-t", "--orpha_genes", metavar="PATH", help="ORPHA-genes file (genes separated by commas)")
args = parser.parse_args()

##############
## MAIN
##############