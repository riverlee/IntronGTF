# Make a intron GTF file from a bed file

For those who would like to calculate the 'intron expression' of a gene by cufflinks, this script helps convert a bed format file (12 columns) which describes the gene structure to a GTF format file where feature 'exon' actually represents the 'intron'.

# Note

* For genes with only one exon will be excluded since there is no intron regin
* For genes with multiple transcripts, we take the longest transcript with most exons to represent for its gene.


# Usage

perl  get_intron_gtf_from_bed.pl in.bed out.gtf

