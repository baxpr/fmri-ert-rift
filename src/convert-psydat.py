#!/usr/bin/env python3
#
# We need a specific psychopy, depending on which version of psychopy
# was used to run the experiment in the first place:
#    ERT1: 2024.2.4
#
# Many of the psychopy function names/calls are very version dependent as well
#
# psychopy does not work on Mac Silicon

import argparse
import os
import psychopy.misc

parser = argparse.ArgumentParser()
parser.add_argument('--psydat_file', required=True)
parser.add_argument('--out_file', required=True)
args = parser.parse_args()


data = psychopy.misc.fromFile(args.psydat_file)

# Method here is to_csv in 2024.1.4 but saveAsWideText in 2024.2.4
data.saveAsWideText(args.out_file, index=False)
