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

def flip(allele)
  case allele
    when "A" then "T"
    when "T" then "A"
    when "G" then "C"
    when "C" then "G"
    else
      $stderr.puts "Error, flipping non-ATCG base."
      exit!
  end
end

def divide_string(str, len=60)
  res = str.scan(/.{#{len}}/)
  l = str.length
  res << str[l - (l % len) .. -1] if l % len > 0
  return res
end

tg_regexp = />(\S+).*gene:(\S+).*/
t_rxp_ind = 1
g_rxp_ind = 2

if ARGV.size == 8 and ARGV[0] == "-m"
  custom_ref = true
  ARGV.shift
  tg_regexp = Regexp.new(">" + ARGV.shift)
  t_rxp_ind = ARGV.shift.to_i
  g_rxp_ind = ARGV.shift.to_i
end

if ARGV.size != 4 then
  $stderr.puts "Usage: #{File.basename($0)} [-m tg_regexp t_ind g_ind] cdna_file gff_file pos_file hap_file > hapiso_fasta"
  $stderr.puts "Generate haplotype and isoform specific FASTA reference."
  $stderr.puts ""
  $stderr.puts "Mandatory arguments:"
  $stderr.puts "  cdna_file: default transcript FASTA file."
  $stderr.puts "  gff_file: GFF file containing structure annotation for transcripts in cdna_file."
  $stderr.puts "  pos_file: file containing, for each transcript, chromosome and positions of SNPs."
  $stderr.puts "  hap_file: file containing, for each transcript, two versions, e.g. suffixed _A and _B, on separate lines with the alleles for the two haplotypes at each position listed in the pos_file (A and B respectively)."
  $stderr.puts ""
  $stderr.puts "Arguments to option -m for use with non-Ensembl cDNA FASTA files:"
  $stderr.puts "  tg_regexp: regular expression matching FASTA entry names, where pairs of brackets"
  $stderr.puts "             are used to capture transcript and gene IDs. Default: \"(\\S+).*gene:(\\S+).*\""
  $stderr.puts "  t_ind:     index of bracket pair that captures the transcript ID. Default: 1."
  $stderr.puts "  g_ind:     index of bracket pair that captures the gene ID. Default: 2."
  $stderr.puts ""
  exit
end

CDNA_FILENAME = ARGV[0]
GFF_FILENAME = ARGV[1]
POS_FILENAME = ARGV[2]
HAP_FILENAME = ARGV[3]

$stderr.print "Reading FASTA file..."
fasta = {} # tid => [gid, seq]
File.open(CDNA_FILENAME) do |f|
  s = f.gets.chomp
  begin
    m = tg_regexp.match(s)
    tid = m[t_rxp_ind]
    gid = m[g_rxp_ind]
    seq = f.gets.chomp
    while(! f.eof? and (s = f.gets.chomp)[0,1] != '>') do seq = seq + s end
    fasta[tid] = [gid, seq]
  end while not f.eof?
end
$stderr.puts "done."


$stderr.print "Reading GFF file..."
gff = {}
strand = {}
File.open(GFF_FILENAME) do |f|
  s = []
  begin
    s = f.gets.split("\t") if s.empty? or s[2]=="gene"
    if s[2] == "gene" then next end
    if s[2] == "mRNA"
      tid = s[8].sub(/.*ID=([^;]*);.*/, '\1').chomp
      s = f.gets.split "\t"
      strand[tid] = s[6]
      positions = []
      begin
        positions << (s[3].to_i..s[4].to_i)
        s = f.gets.split "\t"
      end while !f.eof? and s[2] == "exon"
      gff[tid] = positions.sort do |a,b| a.first <=> b.first end # store ranges in order
    end
  end while !f.eof?
end
$stderr.puts "done."

f_pos = File.open(POS_FILENAME)
f_hap = File.open(HAP_FILENAME)


$stderr.print "Producing haplotype specific FASTA file..."

tids_without_snps = fasta.keys

begin
  l_pos = f_pos.gets.chomp.split(/[[:space:]]/)
  tid = l_pos[0]
  tids_without_snps.delete(tid)
  positions = l_pos[2..-1].map!{|x| x.to_i}
#  $stderr.puts "!!! #{tid} #{positions.join(" ")}"
  l_hap = f_hap.gets.chomp.split(/[[:space:]]/)
  if l_hap[0][0..-3] != tid then
    $stderr.puts "Hap and Pos files incompatible. Aborting."
    exit
  end
  hap1 = l_hap[1].split("")
  l_hap = f_hap.gets.chomp.split(/[[:space:]]/)
  if l_hap[0][0..-3] != tid then
    $stderr.puts "Hap and Pos files incompatible. Aborting."
    exit
  end
  hap2 = l_hap[1].split("")
#  $stderr.puts "HAPS: #{hap1} #{hap2}"

  pos_index=0
  cumulshift=0
  transcriptpositions=[]
  intronic_pos=[] # in case we encounter intronic positions (shouldn't happen in future)
  gff[tid].each do |p|
    while pos_index < positions.length and positions[pos_index] <= p.last 
      if positions[pos_index] < p.first
        $stderr.puts "#{positions[pos_index]} intronic in #{tid}"
        intronic_pos << pos_index
        pos_index+=1
        next
      end # in case this is an intronic position (shouldn't happen)
      transcriptpositions << positions[pos_index] - p.first + cumulshift
#      $stderr.puts "#{positions[pos_index]} - #{p.first} + #{cumulshift}"
      pos_index+=1
    end
    cumulshift += p.last-p.first+1
  end

  # in case we encounter intronic positions (shouldn't happen in future)
  intronic_pos.reverse_each{|i| hap1.delete_at(i); hap2.delete_at(i) }

#  $stderr.puts transcriptpositions.join(" ")
#  $stderr.puts hap1.join(" ")

  seq1 = String.new(fasta[tid][1])
  seq2 = String.new(fasta[tid][1])

  transcriptpositions.each_index do |i|
    if strand[tid] == "+"
 #     $stderr.puts "+Changing #{seq1[transcriptpositions[i],1]} to #{hap1[i]}"
 #     $stderr.puts "+Changing #{seq2[transcriptpositions[i],1]} to #{hap2[i]}"
      seq1[transcriptpositions[i],1] = hap1[i] 
      seq2[transcriptpositions[i],1] = hap2[i]
    else # reverse complement
  #    $stderr.puts "-Changing #{seq1[seq1.length-1-transcriptpositions[i],1]} to #{flip(hap1[i])}"
  #    $stderr.puts "-Changing #{seq2[seq2.length-1-transcriptpositions[i],1]} to #{flip(hap2[i])}"
      seq1[seq1.length-1-transcriptpositions[i],1] = flip(hap1[i])
      seq2[seq2.length-1-transcriptpositions[i],1] = flip(hap2[i])
    end
  end
#  $stderr.puts "DONE"
 # $stderr.puts seq1
#  $stderr.puts seq2

  if seq1==seq2
    puts ">#{tid} gene:#{fasta[tid][0]}"
    puts divide_string(seq1)
  else
    puts ">#{tid}_A gene:#{fasta[tid][0]}_A"
    puts divide_string(seq1)

    puts ">#{tid}_B gene:#{fasta[tid][0]}_B"
    puts divide_string(seq2)
  end
end while !f_pos.eof?

tids_without_snps.each do |tid|
  puts ">#{tid} gene:#{fasta[tid][0]}"
  seq = fasta[tid][1]
  puts divide_string(seq)
end
$stderr.puts "done."

