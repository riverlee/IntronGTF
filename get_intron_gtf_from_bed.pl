#!/usr/bin/env perl
###################################
# Author: Jiang (River) Li
# Email:  riverlee2008@gmail.com
# Date:   Sun Mar 31 16:52:45 2013
###################################
use strict;
use warnings;

my ($inbed,$outgtf) = @ARGV;
my $usage="perl $0 gene.bed intron.gft\n";

die $usage if (! $inbed || ! $outgtf || ! -e $inbed);

# Bed file has 12 columns, the 4th columns is the gene name.
# Since gene has multiple isoforms, so for these genes, we take the longest transcript with most exons 
# to stand for the gene.

#1) load bed 
my %bed; #key are 4th column
open(IN,$inbed) or die $!;
my $totalline=0;
while(<IN>){
    s/\r|\n//g;
    $totalline++;
    my @array = split "\t";
    push @{$bed{$array[3]}},[@array];
}
close IN;
print "Total Line:   $totalline\n";
print "Unique Genes: ",scalar(keys %bed),"\n";

# Loop each gene and print out the gtf
open(OUT,">$outgtf") or die $!;
foreach my $gene(sort keys %bed){
    my $arrayref = $bed{$gene};
    my $strand = "+";
    my $start;
    my $end;
    
    my @array = get_longest_transcript_with_most_exons($arrayref);
    if($array[9] == 1){ #blockcount
        print "Gene '$array[3]' has no intron\n";
    }else{
        my @blocksizes = split ",", $array[10]; 
        my @blockstarts = split ",",$array[11];
        if($array[5] eq "+"){  #strand
            for(my $i=0;$i<$#blocksizes;$i++){
                my $intron_start = $array[1]+$blockstarts[$i]+$blocksizes[$i]+1;
                my $intron_end   = $array[1]+$blockstarts[$i+1];
                my $ii=$i+1;
                my $group = "gene_id \"$array[3]\"; transcript_id \"$array[3]\"; exon_number $ii\n";
                print OUT join "\t",($array[0],"intronGTF","exon",$intron_start,$intron_end,$array[4],$array[5],
                                         0,$group);
            }
        }else{
            for(my $i=0;$i<$#blocksizes;$i++){
                my $intron_start = $array[1]+$blockstarts[$i]+$blocksizes[$i]+1;
                my $intron_end   = $array[1]+$blockstarts[$i+1];
                my $ii=$#blocksizes-$i;
                my $group = "gene_id \"$array[3]\"; transcript_id \"$array[3]\"; exon_number $ii\n";
                print OUT join "\t",($array[0],"intronGTF","exon",$intron_start,$intron_end,$array[4],$array[5],
                                         0,$group);
            }
        }
    }
}


sub get_longest_transcript_with_most_exons{
    my ($ref) = @_;
    my $wanted_index=0;
    my $len=0;
    my $exons=0;
    for(my $i=0;$i<@{$ref};$i++){
        my $thislen=$ref->[$i]->[2]-$ref->[$i]->[1];
        if($thislen>$len){
            $wanted_index=$i;
            $len=$thislen;
            $exons=$ref->[$i]->[9];
        }elsif($thislen == $len){
            if($ref->[$i]->[9]>$exons){
                $exons=$ref->[$i]->[9];
                $wanted_index=$i;
            }
        }
    }
    return @{$ref->[$wanted_index]};
}




