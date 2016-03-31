#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;


my $version_num = '1.0';

my $command_line = join(" ", @ARGV);

###### Get  options ###############
my ($gff,$fas,$genome_vcf_file,$out) = '';

GetOptions(	     
	       'gff=s' => \$gff,      
	       'out=s' => \$out
	);
    
if( $gff  eq ''  or  $out eq ''){
    my $usage_message = usage_message($version_num);
	   die "$usage_message\n";
}


#####################################################
## get mRNA information
open (SRC,"$gff");
open (TGT,">$out") || die "No open";
my %gene2mRNA = ();
my %gene2site = ();

while(<SRC>) 
{
	if(m/^#/)
	{
		next;
	}	
    my @line = split(/\t/,$_);
    
    if($line[8] =~ m/([^\r\n]+)/){
    	  $line[8] = $1;
    }    
    my $id='';
    my $name='';
    my $parent = '';
       if($line[8] =~m/ID\=([^\;]+)/){
    		$id=$1;
    	}
    	if($line[8] =~m/Name\=([^\;]+)/){
    		$name =$1;
    	}
    	if($line[8] =~m/Parent\=([^\;]+)/){
    		$parent =$1;
    	}
    
    if (uc($line[2]) eq 'GENE'){
    	 $gene2site{$id}=$line[0].':'.$line[3].'-'.$line[4];    	
    }
    elsif(uc($line[2]) eq 'MRNA'){	
    	#$trans2gene{$id}   = $parent;
    	#$trans_id2name{$id} = $name;
    	$gene2mRNA{$parent}{$id} = 1;
    }
} ## while file

print "Read GFF done\n";

for my $gene(sort keys %gene2site){
	my $site = $gene2site{$gene} ;

	my $gene1 = $gene.'A';
	my $gene2 = $gene.'B';
	
	
	if(exists $gene2mRNA{$gene}){
		my @trans = sort keys %{$gene2mRNA{$gene}};
		for (@trans){
			print TGT $gene1 ,"\t",$_.'A',"\n";
		}
	    for (@trans){
			print TGT $gene2 ,"\t",$_.'B',"\n";
		}
	}
	
}

sub usage_message {
    my $version = shift;
    my $usage_message = "\n Trans extraction $version
Usage: perl get_geneID.pl --gff  <gff>  --out <out>

";
    return $usage_message;
}


 