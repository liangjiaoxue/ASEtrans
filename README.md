# ASEtrans
Codes for allele specific expression

Two Perl scripts were provided to generate the transcriptome and genes list for RSEM. 
get_allele_specific_transcripts.pl introduce variants from a VCF file into transcripts defined in one GFF file.
get_geneID.pl  generate the corresponding lists between allele IDs and transcripts IDs.

-----------------------------------------------------------------------------------------
# Usage:

Prepare vcf file

/usr/local/tabix/latest/bin/bgzip   xxx.vcf

/usr/local/tabix/latest/bin/tabix -p vcf xxx.vcf.gz

Get transcripts

perl get_allele_specific_transcripts.pl --gff  <gff>  --genomefile <genome>  --vcf <vcf>  --out <out>
vcf --> vcf data of hybrid sample
genome --> reference genome fasta file used to generate vcf file
this script requires Vcf.pm


Get gene list file
perl get_geneID.pl --gff  <gff>  --out <out>




# DEMO command line for RSEM
cd RSEM

/usr/local/rsem/latest/rsem-prepare-reference  \

            
   --transcript-to-gene-map P1979_gene_list.txt  \
            
                
       --bowtie2 --bowtie2-path /usr/local/bowtie2/latest/bin/ \
                       
                       
    P1979_allele_trans.fasta  ref/P1979_allele_ref &
    


/usr/local/rsem/latest/rsem-calculate-expression -p 5 --paired-end   \

 --bowtie2 --bowtie2-path /usr/local/bowtie2/latest/bin/   \
 
 --estimate-rspd     --bowtie2-mismatch-rate 0.2   \
 
  RNAseq_R1.fq RNAseq_R2.fq  \
  
  ref/P1979_allele_ref   exp/RNAseq_s1


