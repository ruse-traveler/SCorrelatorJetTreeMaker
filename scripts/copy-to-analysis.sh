#!/bin/bash
# -----------------------------------------------------------------------------
# 'copy-to-analysis.sh'
# Derek Anderson
# 01.06.2023
#
# Script to automate copying files
# over to the sPHENIX analysis
# repository.
# -----------------------------------------------------------------------------

# declare filelist
declare -a files_to_copy

# top directory to copy from/to
copy_from="/sphenix/user/danderson/eec/SCorrelatorJetTreeMaker"
copy_to="/sphenix/user/danderson/analysis/AndersonAnalysisModules/ColdQcdCorrelatorAnalysis/SCorrelatorJetTreeMaker"

# what files to copy
files_to_copy[0]="README.md"
files_to_copy[1]="Fun4All_RunCorrelatorJetTreeMaker.C"
files_to_copy[2]="RunCorrelatorJetTreeMaker.sh"
files_to_copy[3]="macros/MergeFiles.C"
files_to_copy[4]="macros/Fun4All_RunCorrelatorJetTreeMakerOnCondor.C"
files_to_copy[5]="scripts/MergeFiles.rb"
files_to_copy[6]="scripts/MergeCondorFiles.rb"
files_to_copy[7]="scripts/SwitchToCorrectBuild.sh"
files_to_copy[8]="scripts/wipe-source.sh"
files_to_copy[9]="scripts/copy-to-analysis.sh"
files_to_copy[10]="scripts/makeCondorJobs.py"
files_to_copy[11]="scripts/RunCorrelatorJetTreeMakerOnCondor.sh"
files_to_copy[12]="src/SCorrelatorJetTreeMaker.cc"
files_to_copy[13]="src/SCorrelatorJetTreeMaker.h"
files_to_copy[14]="src/SCorrelatorJetTreeMaker.io.h"
files_to_copy[15]="src/SCorrelatorJetTreeMaker.jet.h"
files_to_copy[16]="src/SCorrelatorJetTreeMaker.sys.h"
files_to_copy[17]="src/SCorrelatorJetTreeMakerLinkDef.h"
files_to_copy[18]="src/SCorrelatorJetTreeMaker.cst.h"
files_to_copy[19]="src/SCorrelatorJetTreeMaker.evt.h"
files_to_copy[20]="src/autogen.sh"
files_to_copy[21]="src/configure.ac"
files_to_copy[22]="src/Makefile.am"
files_to_copy[23]="src/sphx-build"
files_to_copy[24]="scripts/HAddFilesFromList.rb"
files_to_copy[25]="scripts/HAddFilesFromPattern.rb"
files_to_copy[26]="scripts/MergeFilesFromList.rb"
files_to_copy[27]="scripts/SplitFileLists.rb"

# do copying
# TODO: automate detection/creation of sub-directories
(( nFile=0 ))
for file in ${files_to_copy[@]}; do
  source_file="$copy_from/$file"
  target_file="$copy_to/$file"
  rsync -azP $source_file $target_file
  (( nFile++ ))
done

# delete array
unset files_to_copy

# end -------------------------------------------------------------------------
