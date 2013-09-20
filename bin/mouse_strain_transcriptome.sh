#!/bin/sh

## Script to create strain-specific mouse trancriptome and genome fastas
## Tweak variables below as required

ENSEMBLVERSION=70
SANGERFTP=ftp-mouse.sanger.ac.uk
SANGERDIR=REL-1303-SNPs_Indels-GRCm38
SANGERVERSION=v3
SNPFILE=mgp.$SANGERVERSION.snps.rsIDdbSNPv137.vcf.gz
INDELFILE=mgp.$SANGERVERSION.indels.rsIDdbSNPv137.vcf.gz
GENOMEFASTA=GRCm38_68.fa
GENOMEVERSION=GRCm38
STRAINS=(129P2 129S1 129S5 AJ AKRJ BALBcJ C3HHeJ C57BL6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ)

## Edit below at your own risk

echo "Welcome." >&2
echo "I will produce strain-specific genomes and transcriptomes using the Sanger mouse genomes project data" >&2
echo "Check the following software is installed:" >&2
echo -e "- MMSEQ (extract_transcripts and offsetGTF)" >&2
echo -e "- samtools" >&2
echo -e "- tabix" >&2
echo -e "- vcftools (vcf-isec)" >&2
echo -e "- GATK in path GATK_HOME" >&2
echo -e "- SnpEff in path SNPEFF_HOME" >&2
echo -e "- Picard in path PICARD_HOME\n" >&2

echo "I'm about to produce quite a lot of files; check you are in an appropriate working directory." >&2

while true; do
    read -p "Proceed? " yn
    case $yn in
        [Yy]* ) echo "" >&2; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

echo -e "ENSEMBLVERSION:\t$ENSEMBLVERSION" >&2
echo -e "SANGERFTP:\t$SANGERFTP" >&2
echo -e "SANGERDIR:\t$SANGERDIR" >&2
echo -e "SANGERVERSION:\t$SANGERVERSION" >&2
echo -e "SNPFILE:\t$SNPFILE" >&2
echo -e "INDELFILE:\t$INDELFILE" >&2
echo -e "GENOMEFASTA:\t$GENOMEFASTA" >&2
echo -e "GENOMEVERSION:\t$GENOMEVERSION" >&2
echo -e "STRAINS: ${STRAINS[@]}" >&2
echo "" >&2

for cmd in samtools tabix vcf-isec extract_transcripts offsetGTF; do
  if ! type "$cmd" 2> /dev/null; then
    if [[ "$cmd" == "extract_transcripts" || "$cmd" == "offsetGTF" ]]; then
      if [[ "`uname -s`" == "Linux" ]]; then
        MSSUF="-linux"
      elif [[ "`uname -s`" == "Darwin" ]]; then
        MSSUF="-mac"
      fi
      if ! type "$cmd${MSSUF}" 2> /dev/null; then
        echo "Error: please install $cmd" >&2
        exit 1
      else
        continue
      fi
    else
      echo "Error: please install $cmd" >&2
      exit 1
    fi
  fi
done
echo "" >&2

if [[ ! -d "$GATK_HOME" ]]; then
  echo "Error: please ensure the environment variable GATK_HOME contains the path to GATK" >&2
  exit 1;
fi

if [[ ! -d "$SNPEFF_HOME" ]]; then
  echo "Error: please ensure the environment variable SNPEFF_HOME contains the path to SnpEff" >&2
  exit 1;
fi

if [[ ! -d "$PICARD_HOME" ]]; then
  echo "Error: please ensure the environment variable PICARD_HOME contains the path to Picard" >&2
  exit 1;
fi

if [[ ! -e $SNPFILE ]]; then 
  wget ftp://$SANGERFTP/$SANGERDIR/$SNPFILE
  if [[ ! $? -eq 0 ]]; then
    echo "Error downloading file" >&2
    exit 1
  fi
else
  echo "WARNING: found $SNPFILE, not downloading. Delete and re-run to force download." >&2 
fi
if [[ ! -e $INDELFILE ]]; then
  wget ftp://$SANGERFTP/$SANGERDIR/$INDELFILE
  if [[ ! $? -eq 0 ]]; then
    echo "Error downloading file" >&2
    exit 1
  fi
else
  echo "WARNING: found $INDEFILE, not downloading. Delete and re-run to force download." >&2 
fi
if [[ ! -e $GENOMEFASTA ]]; then
  wget ftp://$SANGERFTP/ref/$GENOMEFASTA
  if [[ ! $? -eq 0 ]]; then
    echo "Error downloading file" >&2
    exit 1
  fi
else
  echo "WARNING: found $GENOMEFASTA, not downloading. Delete and re-run to force download." >&2 
fi
if [[ ! -e $Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf.gz ]]; then
  wget ftp://ftp.ensembl.org/pub/release-$ENSEMBLVERSION/gtf/mus_musculus/Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf.gz
  if [[ ! $? -eq 0 ]]; then
    echo "Error downloading file" >&2
    exit 1
  fi
else
  echo "WARNING: found $Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf.gz, not downloading. Delete and re-run to force download." >&2 
fi


echo -e "\n=== All required input files found; proceeding to generate transcriptome and genomes  ===\n" >&2
sleep 4

if [[ ! -e $GENOMEFASTA.fai ]]; then
  samtools faidx $GENOMEFASTA
fi

gunzip -c Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf.gz > Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf

java -Xmx4g -jar ${PICARD_HOME}/CreateSequenceDictionary.jar REFERENCE=$GENOMEFASTA OUTPUT=`basename $GENOMEFASTA .fa`.dict

extract_transcripts${MSSUF} $GENOMEFASTA Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf > C57BL6_transcriptome$ENSEMBLVERSION.fa

for STRAIN in ${STRAINS[@]}; do

echo -e "\n*** *** DOING STRAIN $STRAIN *** ***\n" >&2

vcf-subset -a -e -c ${STRAIN} $SNPFILE | bgzip -c > mgp.$SANGERVERSION.${STRAIN}.snps.annot.reformat.vcf.gz
zcat mgp.$SANGERVERSION.${STRAIN}.snps.annot.reformat.vcf.gz | java -jar ${SNPEFF_HOME}/SnpSift.jar filter "( GEN[0].FI = 1 )" | bgzip -c > mgp.$SANGERVERSION.${STRAIN}.snps.annot.reformat.filter.vcf.gz
tabix -p vcf mgp.$SANGERVERSION.${STRAIN}.snps.annot.reformat.filter.vcf.gz

vcf-subset -a -e -c ${STRAIN} $INDELFILE | bgzip -c > mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.vcf.gz
zcat mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.vcf.gz | java -jar ${SNPEFF_HOME}/SnpSift.jar filter "( GEN[0].FI = 1 )" | bgzip -c > mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.filter.vcf.gz
tabix -p vcf mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.filter.vcf.gz

# Merge SNPs and indels
# Remove redundant REF and ALT substrings and adjust POS
# Sort (this is needed because POS adjustment can destroy coordinate-sorting)
# Ensure only the first of a sequential set of overlapping indels is retained
# Two adjacent deletions must be separated by at least one base (otherwise there is an inconsistency, so I add a base to the left of each deletion range for checking for overlap even though this might remove a legit deletion that immediately follows an insertion)
# Ensure only first entry is kept if same POSition appears multiple times sequentially
# Set all FILTER values to PASS
# All this is necessary to make sure GATK's FastaAlternateReferenceMaker uses exactly what is in the VCF
# which is later used by offsetGTF
vcf-isec -n +1 mgp.$SANGERVERSION.${STRAIN}.snps.annot.reformat.filter.vcf.gz mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.filter.vcf.gz | 
awk 'BEGIN {OFS="\t"} /^#/ { print; next } {MIN=(length($5)<length($4))?length($5):length($4);
 if(MIN>1) {
    $4=substr($4,MIN);
    $5=substr($5,MIN);
    $2=$2+MIN-1
 }
 print
}' | vcf-sort -c | awk 'BEGIN {OFS="\t";CURPOS=0;CURSTART=0;CUREND=0;THISSTART=0;THISEND=0} /^#/ { print; next } 

{if(length($4) > length($5)) { 
  THISSTART=$2+length($5)-1; THISEND=$2+length($4)-1
}}
{if(length($5) > length($4)) {
  THISSTART=$2+length($4)-1; THISEND=$2+length($4)-1
}}
{if( length($4) != length($5) && ((THISSTART<=CUREND && THISEND>=CURSTART))) {
  CURSTART=THISSTART
  CUREND=THISEND
  next
}}
CURPOS==$2 {
  next;
}
{CURPOS=$2}
{if( length($4) != length($5)) {CURSTART=THISSTART; CUREND=THISEND }} $7!="PASS" {$7="PASS" } {print}' | bgzip -c > mgp.$SANGERVERSION.${STRAIN}.snps_indels.annot.reformat.filter.vcf.gz
tabix -p vcf mgp.$SANGERVERSION.${STRAIN}.snps_indels.annot.reformat.filter.vcf.gz

java -Xmx4g -jar ${GATK_HOME}/GenomeAnalysisTK.jar -R $GENOMEFASTA -T FastaAlternateReferenceMaker -o $GENOMEVERSION.${STRAIN}.gatk.fa --variant mgp.$SANGERVERSION.${STRAIN}.snps_indels.annot.reformat.filter.vcf.gz

# Fix chromosome names (which are now just numerical)
CHRS=(`grep '>' $GENOMEFASTA | awk '{print $1}'`)
unset AWKEXP
AWKEXP="/^[^>]/ { print; next } "
for I in `seq 1 ${#CHRS[@]}`; do
  AWKEXP+="/^>${I}$/ { print \"${CHRS[$(( $I - 1 ))]}\" } "
done

awk "$AWKEXP" $GENOMEVERSION.${STRAIN}.gatk.fa > $GENOMEVERSION.${STRAIN}.fa

bgzip -d -c mgp.$SANGERVERSION.${STRAIN}.snps_indels.annot.reformat.filter.vcf.gz | awk 'BEGIN { OFS="\t" } /^#/ {print;next } { if(length($4) != length($5)) { print }}' > mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.filter2.vcf

offsetGTF${MSSUF} Mus_musculus.$GENOMEVERSION.$ENSEMBLVERSION.gtf mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.filter2.vcf > Mus_musculus.${STRAIN}.$GENOMEVERSION.$ENSEMBLVERSION.gtf

rm mgp.$SANGERVERSION.${STRAIN}.indels.annot.reformat.filter2.vcf

extract_transcripts${MSSUF} $GENOMEVERSION.${STRAIN}.fa Mus_musculus.${STRAIN}.$GENOMEVERSION.$ENSEMBLVERSION.gtf > ${STRAIN}_transcriptome$ENSEMBLVERSION.fa
done
