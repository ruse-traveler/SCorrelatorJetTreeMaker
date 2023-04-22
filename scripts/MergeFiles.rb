#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'MergeFiles.rb'
# Derek Anderson
# 04.22.2023
#
# For merging files using 'hadd_files.C'
# -----------------------------------------------------------------------------

# modules to use
require 'fileutils'

# input parameters
in_path = "./intermediate_merge/"
in_pref = "sPhenixG4_forTrackCutStudy_embedOnly0t99_g4svtxeval.d"
in_suff = "m12y2022.root"

# output parameters
out_list = "testingRubyScript.list"
out_file = "testingRubyScript.root"

# create input matching pattern
in_pattern = in_path + "/" + in_pref + "*" + in_suff
in_pattern.gsub!("//", "/")

# create list of files to merge
File.open(out_list, "w") { |out|
  Dir[in_pattern].each do |file|
    puts file
    out.puts file
  end
}

# grab number of files to merge
num_files = Dir[in_pattern].size

# merge files
exec("root -b -q \'MergeFiles.C(#{num_files}, \"#{out_list}\", \"#{out_file}\")\'")

# end -------------------------------------------------------------------------
