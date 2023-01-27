#!/bin/bash
# 'RunMakeCondorJobs.sh'
# Derek Anderson
# 01.25.2023
#
# Short shell script to run
# the makeCondorJobs.py
# script

python3 makeCondorJobs.py --inputType PYTHIA8_PP_MB --run 50 --truth --calo --g4hit

# end -------------------------------------------------------------------------
