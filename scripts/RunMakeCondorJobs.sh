#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunMakeCondorJobs.sh'
# Derek Anderson
# 01.25.2023
#
# Short shell script to run
# the makeCondorJobs.py
# script
# -----------------------------------------------------------------------------

python3 makeCondorJobs.py --inputType JET_20GEV --run 6 --embed pau --truth --calo --g4hit

# end -------------------------------------------------------------------------
