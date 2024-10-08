#!/bin/bash
# -----------------------------------------------------------------------------
# @file   RunCorrelatorJetTreeMakerOnCondor.sh
# @author Derek Anderson
# @date   01.06.2022
#
# Script to run correlator jet tree module via condor.
# Derived from code by Cameron Dean and Antonio Silva.
# -----------------------------------------------------------------------------

# grab home directory and set up environment
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh

# initialize no. of events to run and verbosity
nEvents=0
verbose=0

# create array of file lists
inputFiles="{"
for fileList in $@; do
  inputFiles+="\"${fileList}\","
done
inputFiles=${inputFiles::-1}
inputFiles+="}"

# run script
echo running: RunCorrelatorJetTreeMakerOnCondor.sh $*
root -b -q "Fun4All_RunCorrelatorJetTreeMakerOnCondor.C(${inputFiles},$nEvents,$verbose)"
echo Script done

# end -------------------------------------------------------------------------
