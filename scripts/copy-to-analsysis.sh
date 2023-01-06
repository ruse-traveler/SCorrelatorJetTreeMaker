#!/bin/bash
# 'copy-to-analysis.sh'
# Derek Anderson
# 01.06.2023
#
# Script to automate copying files
# over to the sPHENIX analysis
# repository.

# declare filelist
declare -a files_to_copy

# top directory to copy from/to
copy_from="/sphenix/user/danderson/eec"
copy_to="/sphenix/user/danderson/test"

# what files to copy
files_to_copy[0]="Fun4All_ForCorrelatorJetTree.C"
files_to_copy[1]="RunCorrelatorJetTree.sh"
files_to_copy[2]="macros/MergeFiles.C"
files_to_copy[3]="scripts/MergeFiles.sh"
files_to_copy[4]="scripts/SwitchToCorrectBuild.sh"
files_to_copy[5]="scripts/wipe-source.sh"
files_to_copy[6]="src/SCorrelatorJetTree.cc"
files_to_copy[7]="src/SCorrelatorJetTree.h"
files_to_copy[8]="src/SCorrelatorJetTree.io.h"
files_to_copy[9]="src/SCorrelatorJetTreeLinkDef.h"
files_to_copy[10]="src/autogen.sh"
files_to_copy[11]="src/configure.ac"
files_to_copy[12]="src/Makefile.am"

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
unset -a files_to_copy

# end -------------------------------------------------------------------------
