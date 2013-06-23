for c in bam2hits hitstools mmseq mmdiff mmcollapse extract_transcripts t2g_hits 'fastagrep.sh -h' testregexp.rb filterGTF.rb haploref.rb ensembl_gtf_to_gff.pl; do
  echo '`'$c'`' | sed -e 's/ -h//'
  echo ""
  if [[ $c != "ensembl_gtf_to_gff.pl" ]]; then 
    ./$c 2>&1 | grep -E -v "^Error|illegal" | sed -E 's/^/    /' 
  else
    echo "    convert gtf file from Ensembl to gff3 file"
    echo ""
    echo '    Usage: ensembl_gtf_to_gff.pl gtf_file > gff_file'
  fi
done
