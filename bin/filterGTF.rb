#!/usr/bin/env ruby
##    Copyright 2010 Ernest Turro
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA


if ARGV.size != 2 then
  $stderr.puts "Usage: #{File.basename($0)} cdna_file gtf_file > new_gtf_file"
  $stderr.puts "Filter GTF file according to transcripts in cDNA file."
  $stderr.puts ""
  $stderr.puts "Mandatory arguments:"
  $stderr.puts "  cdna_file: Ensembl transcript FASTA file (e.g. Homo_sapiens.GRCh37.56.cdna.ref.fa)."
  $stderr.puts "  gtf_file: Ensembl GTF file (e.g. Homo_sapiens.GRCh37.56.gtf)."
  $stderr.puts ""
  exit
end

CDNA_FILENAME = ARGV[0]
GTF_FILENAME = ARGV[1]

rxp = Regexp.new('>(E\S+).*')

tids = {} # hash with transcript ID keys
File.open(CDNA_FILENAME) do |f|
  while(!f.eof?) do
    s = f.gets
    if s[0,1] == '>'
      tids[rxp.match(s)[1]]=""
    end
  end
end

rxp = Regexp.new('.*transcript_id \"(E\S+)\".*')

File.open(GTF_FILENAME) do |f|
  while(!f.eof?) do
    s = f.gets
    unless tids[rxp.match(s)[1]].nil?
      print s
    end
  end
end

