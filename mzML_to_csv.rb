#!/usr/bin/env ruby 
require 'mspire/mzml' #install by "gem install mspire"

if ARGV.size == 0
  puts "usage: feed me an mzML file and I'll spit back the ms1 in csv\n\twith mz,intensity,retention_time for each scan, unseparated by headers"
  exit
end

def write_to_csv(io_obj, spectrum)
  rt = spectrum.retention_time
  spectrum.peaks do |mz, int|
    io_obj.puts [mz,int,rt].join(',')
  end
end

ARGV.each do |file|
  outfile_base = File.absolute_path(file).sub('.mzML','')
  Mspire::Mzml.open(ARGV.first) do |mzml|
    File.open(outfile_base + '_ms1.csv', 'w') do |ms1_out|
			mzml.each do |spectrum|
				if spectrum.ms_level == 1
					write_to_csv(ms1_out, spectrum)
        end
      end
    end
  end
end
