#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'QuickHAdd.rb'
# Derek Anderson
# 05.15.2023
#
# Makes hadding easier...
# -----------------------------------------------------------------------------

# i/o parameters
out_file = "correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA.d11m10y2023.root"
in_files = [
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk0.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk1.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk2.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk3.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk4.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk5.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk6.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk7.d11m10y2023.root",
  "./condor/intermediate_merge/pp200py8jet10run8_roughCutsWithTrkTupleQA_2023oct10/correlatorJetTree.pp200py8jet10run8_roughCutsWithTrkTupleQA_goodFiles_chunk8.d11m10y2023.root"
]

# create input argument
arg_input = ""
num_input = in_files.size
in_files.each_with_index do |file, iFile|
  arg_input += file
  arg_input += " " if iFile + 1 != num_input
end

# run hadd
exec("hadd #{out_file} #{arg_input}")

# end -------------------------------------------------------------------------
