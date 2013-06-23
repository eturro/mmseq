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

tg_regexp = />(\S+).*gene:(\S+).*$/
t_rxp_ind = 1
g_rxp_ind = 2

if ARGV.size == 5 and ARGV[0] == "-m"
  ARGV.shift
  tg_regexp = Regexp.new(">" + ARGV.shift + "$")
  t_rxp_ind = ARGV.shift.to_i
  g_rxp_ind = ARGV.shift.to_i
else
  $stderr.puts "Usage: #{File.basename($0)} -m tg_regexp t_ind g_ind cdna_file"
  $stderr.puts "Test regular expression patterns for use with `bam2hits` and `haploref.rb`."
  $stderr.puts ""
  $stderr.puts "Arguments to flag -m:"
  $stderr.puts "  tg_regexp: regular expression matching FASTA entry names, where pairs of brackets"
  $stderr.puts "             are used to capture transcript and gene IDs. Default: \"(\\S+).*gene:(\\S+).*\""
  $stderr.puts "  t_ind:     index of bracket pair that captures the transcript ID. Default: 1."
  $stderr.puts "  g_ind:     index of bracket pair that captures the gene ID. Default: 2."
  $stderr.puts ""
  $stderr.puts "  cdna_file: reference FASTA file."
  exit
end

CDNA_FILENAME=ARGV[0]

File.open(CDNA_FILENAME) do |f|
  s = f.gets.chomp
  begin
    m = tg_regexp.match(s)
    tid = m[t_rxp_ind]
    gid = m[g_rxp_ind]
    puts "Transcript ID: #{tid}; Gene ID: #{gid}."
    seq = f.gets.chomp
    while(! f.eof? and (s = f.gets.chomp)[0,1] != '>') do seq = seq + s end
  end while not f.eof?
end


