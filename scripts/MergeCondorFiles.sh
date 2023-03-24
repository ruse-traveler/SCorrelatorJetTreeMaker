#!/bin/bash
# -----------------------------------------------------------------------------
# 'MergeCondorFiles.sh'
# Derek Anderson
# 03.24.2023
#
# For merging files in chunks using 'hadd_files.C'.
# (Ideal for condor output!)
# -----------------------------------------------------------------------------

# declare arrays
declare -a chunks
declare -a labels

# input parameters
iPath="./condor/individual_files/pp200py8jet30run6_trksAndChrgPars_2023mar23"
iPref="outputData_CorrelatorJetTree_0"
iSuff=".root"

# output parameters
oPath="."
oPref="correlatorJetTree.pp200py8jet30run6_trksAndChrgPars_"
oDate="d24m2y2023"
oList=".list"
oRoot=".root"

# to select chunks
chunks[0]="00"
chunks[1]="01"
chunks[2]="02"
chunks[3]="03"
chunks[4]="04"
chunks[5]="05"
chunks[6]="06"
chunks[7]="07"
chunks[8]="08"
chunks[9]="09"
chunks[10]="10"
chunks[11]="11"
chunks[12]="12"
chunks[13]="13"
chunks[14]="14"
chunks[15]="15"
chunks[16]="16"
chunks[17]="17"
chunks[18]="18"
chunks[19]="19"

# for output files
labels[0]="0000to0099"
labels[1]="0100to0199"
labels[2]="0200to0299"
labels[3]="0300to0399"
labels[4]="0400to0499"
labels[5]="0500to0599"
labels[6]="0600to0699"
labels[7]="0700to0799"
labels[8]="0800to0899"
labels[9]="0900to0999"
labels[10]="1000to1099"
labels[11]="1100to1199"
labels[12]="1200to1299"
labels[13]="1300to1399"
labels[14]="1400to1499"
labels[15]="1500to1599"
labels[16]="1600to1699"
labels[17]="1700to1799"
labels[18]="1800to1899"
labels[19]="1900to1999"

# loop over chunks
(( iChunk=0 ))
for chunk in ${chunks[@]}; do

  # create list of files
  files="$iPath/$iPref$chunk*$iSuff"
  list="$oPath/$oPref${labels[$iChunk]}.$oDate$oList"
  nFiles=$(ls -1 $files | wc -l)
  ls $files > $list

  # merge files and increment counters
  output="$oPath/$oPref${labels[$iChunk]}.$oDate$oRoot"
  root -b -q "MergeFiles.C($nFiles, \"$list\", \"$output\")"
  (( iChunk++ ))

done

# delete arrays
unset chunk
unset label

# end -------------------------------------------------------------------------
