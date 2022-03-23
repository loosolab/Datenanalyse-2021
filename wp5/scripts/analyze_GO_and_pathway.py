#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser(description='Evaluate differntial binding')
parser.add_argument('--in-dir', metavar="PATH/TO/DIR", help='Path to directory where the motif discovery runs are stored', required=True )
parser.add_argument('--motifs', metavar='MOTIF_CLUSER.yml', help='Output of the motif clustering with TOBIAS', required=True)
parser.add_argument('--out', metavar="FILENAME_prefix", help="Prefix of how the output files should be named.", required=True)
args = parser.parse_args()
