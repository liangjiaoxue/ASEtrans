#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;
use Vcf;


#Usage: perl get_allele_specific_transcripts.pl --gff  <gff>  --genomefile <genome>  --vcf <vcf>  --out <out>

my $version_num = '1.0';
my $usage_message = usage_message($version_num);

my $command_line = join(" ", @ARGV);

###### Get  options ###############
my ($gff,$fas,$genome_vcf_file,$out) = '';

GetOptions(	     
	       'gff=s' => \$gff,
	       'genomefile=s' => \$fas,
	       'vcf=s' => \$genome_vcf_file,	      
	       'out=s' => \$out
);
    



if( $gff  eq '' or $fas  eq  '' or $genome_vcf_file eq '' or  $out eq ''){
    my $usage_message = usage_message($version_num);
	   die "$usage_message\n";
}

#####################################################
## get mRNA information
open (SRC,"$gff");
my %scafffold2mRNA = ();
my %trans2exon = ();
my %mRNA2gene = ();
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
    	next;	
    }
    elsif(uc($line[2]) eq 'MRNA'){	
    	$scafffold2mRNA{$line[0]}{$id} = 1;
    	$mRNA2gene{$id} = $parent;
    }
    else{
    	# exon
    	if(not exists $trans2exon{$parent}){
   	  	          my   @buff=();
   	  	          push @buff, [@line];
   	  	          $trans2exon{$parent} = \@buff;
   	    }
   	    else{
   	  	          push @{$trans2exon{$parent}}, [@line];
   	    }
    }
} ## while file

print "Read GFF done\n";


##############################################################
### merge exons from UTR and CDS
for(keys  %trans2exon){
	my @ttt=@{$trans2exon{$_}};
	#sort
	@ttt=sort{$a->[3]<=>$b->[3]}@ttt;
	#fusion
	my @sss=();
	for(@ttt){
		my @gff=@{$_};
		push @sss,$gff[3].'-'.$gff[4];
	}
	my $test=join('|',@sss);
	  @sss=split(/\-/,$test);
	  
	  for(my $i=0;$i<@sss;$i++){
	  	if($sss[$i]=~ m/\|/){
	  		my @ccc=split(/\|/,$sss[$i]);
	  		if($ccc[0] == $ccc[1]-1){
	  			$sss[$i]='';
	  		}
	  	}
	  }
	  
	  $test=join('-',@sss);
	  $test=~s/\-+/\-/g; ##  + impotant
	  @sss=split(/\|/,$test);
	  my @output=();
	  my @temp=@{$ttt[0]};
	  
	 for(@sss){
	    	my @ccc=split(/\-/,$_);
	  	    $temp[3]=$ccc[0];
	  	    $temp[4]=$ccc[1];
	  	    push @output,[@temp];
	  	    
	  	    ###### output exon
	  	    $temp[2] =   'exon';
	  	  #  print TGT join("\t", @temp),"\n"; 
	 }
	 $trans2exon{$_} = \@output;
	 
}

print "Merge exon done\n";
####################################################################

## read genome fasta 
my %hash_vcf =();
open (TGT,">$out");
my $hetero_num   = 0;
my $conflict_num = 0;
my $snp_site = 0;

my $seq_obj = new Bio::SeqIO(-format => 'largefasta',                                   
                          -file   => $fas); 
                          
my $pseq ;
my $paternal = 'M';

while( $pseq=$seq_obj->next_seq()){  
	
    my  $scaffold= $pseq->display_id;
    print $scaffold,"\n";
     
    if(not exists $scafffold2mRNA{$scaffold}){
     	  next;
    }
     ########################################################################################
     $paternal = 'M';        
     my $addr_out = &build_snp_hash($scaffold, $paternal);
     %hash_vcf = %{$addr_out};
     &write_transcripts ($scaffold,$paternal);
       
     
     #################################################################################
     %hash_vcf  = ();
     $paternal = 'P';   
     $addr_out = &build_snp_hash($scaffold, $paternal); 
     %hash_vcf = %{$addr_out};
     &write_transcripts ($scaffold,$paternal);
}# while

print "heterozygous total \t $hetero_num\n" ;
print "heterozygous conflict \t $conflict_num\n" ;


sub write_transcripts{
	  my ($scaffold, $paternal) = @_;
	
     my @mRNA_names = sort keys %{$scafffold2mRNA{$scaffold}};
     
     for my $mRNA_in (@mRNA_names){     	
     	  my $gene_out = $mRNA2gene{$mRNA_in };
     	  
     	  $snp_site = 0;
     	  if(not exists  $trans2exon{$mRNA_in}){
     	  	  next;
     	  }
     	  
     	  my @exon    = @{$trans2exon{$mRNA_in}};
     	  @exon       = sort {$a->[3] <=> $b->[3]} @exon;
     	  
     	  my $strand  = $exon[0][6];
     	  my $seq_out = '';
     	  if($strand eq '+'){
     	  	   $seq_out = &process_plus_mRNA(\@exon);
     	  }
     	  else{
     	  	   $seq_out  = &process_minus_mRNA(\@exon);
     	  }
     	  #print "to print \n";
     	  print TGT ">$mRNA_in".$paternal." $snp_site\n",$seq_out,"\n";
     } # mRNA
}

###################################################################

sub process_plus_mRNA{
	my @exon  = @{shift @_};
	
	my @out = ();
	
	for  (@exon) {
		my @line = @{$_};
		
		my $start  = $line[3];
		my $end    = $line[4];
		my $chr    = $line[0];
		my $seq_hit = $pseq->subseq($start,$end);  
		my $seq_order_addr = &get_all_site_from_hash($seq_hit,$chr,$start,$end); 	
		my @seq_out_order  = @{$seq_order_addr};
		
		## translte into sequence
		my $snp_out = &process_plus(\@seq_out_order);
		
		push @out,$snp_out;
	}
	
	my $string_out = join('', @out );
	return $string_out;
}

sub process_minus_mRNA{
    my @exon  = @{shift @_};
	
	my @out = ();
	
	for  (@exon) {
		my @line = @{$_};
		
		my $start  = $line[3];
		my $end    = $line[4];
		my $chr    = $line[0];
		my $seq_hit = $pseq->subseq($start,$end);  
		my $seq_order_addr = &get_all_site_from_hash($seq_hit,$chr,$start,$end); 	
		my @seq_out_order = @{$seq_order_addr};
		
		## translte into sequence
		my $snp_out = &process_minus(\@seq_out_order);
		
		push @out,$snp_out;
	}
	@out = reverse (@out);
	my $string_out = join('', @out );
	return $string_out;
}


		    
sub process_plus{
	my @all_in = @{shift @_};
	my @out = ();
	for my $each_site (@all_in){
		my $out_string = '';
		if(length($each_site) == 0){
				 $out_string = '';
		}
		elsif(length($each_site) == 1){
			   $out_string = $each_site;
		}
		else{ ## multiple or insetion
		        my @this_site  = split(/\,/,$each_site);
		        
		        if(@this_site>1){
		        	 die "Many hit this site\n";
		        }
		        $out_string = $this_site[0];
		}## multiple options
		
				        
	   if( $out_string =~ m/\-/){
		        $out_string = '';
	   }
		
		push @out, $out_string;
	} ## for
	my $out_final = join('',@out);
	
	return $out_final;
}


############################

sub process_minus{
	my @all_in = @{shift @_};
	my @out = ();
	for my $each_site (@all_in){
			my $out_string = '';
			if(length($each_site) == 0){
				 $out_string = '';
			}
			elsif(length($each_site) == 1){
				   $out_string = $each_site;
				   $out_string =~ tr/ACGT/TGCA/; # complementary
			}
			else{ ## multiple or insetion
			        my @this_site  = split(/\,/,$each_site);
			        
			        
			        if(@this_site>1){
			        	 die "Many hit this site\n";
			        }
			        
			        $out_string = $this_site[0];		        
			        $out_string = re_co($out_string);
			} #
			
		   if( $out_string =~ m/\-/){
			        $out_string = '';
		   }
		   push @out, $out_string;
	} # for loop
	
    @out = reverse (@out);
	my $out_final = join('',@out) ;
	return $out_final;
}
	
#########################################



sub get_all_site_from_hash{
	# check every nucleotide
    my ($ref_dna,$seqid,$start,$end) = @_;
    my @ref_all = split(//,$ref_dna);
    
    my @out_all = ();

    for(my $i=0; $i<@ref_all;$i++){
    		my $site_new = $start+$i;
    		my $ref_new  = $ref_all[$i];
    		my $index_new = $seqid.':'.$site_new.':'.$ref_new;
    		
    		if(not exists $hash_vcf{$index_new}){
    			 # no SNP
    			  my $ooo = $ref_new;
    			  push @out_all, $ooo;
    		}
    		else{
    			  # there is snp 
    			  $snp_site++;
    			 my  %hash_nt =   ();
    			 
    			 $hash_nt{$hash_vcf{$index_new}} = 1; # alt copy
    			 
    			 my @hits = sort keys %hash_nt;
    			 if(@hits >1){
    			 	 print "Really many hits :",join("\t",@hits),"\n";
    			 }
    			 
		    	 push @out_all, $hits[0];
    		} # there is snp
    }## for 
    return \@out_all;
}# sub 
    		
      


sub build_snp_hash{
	my ($seqid,$paternal) = @_;

	 my $vcf = Vcf->new(file=>$genome_vcf_file,region=>$seqid);
     $vcf->parse_header();
     
     my %vcf_hash = ();
     ### read SNP at region
     while (my $x=$vcf->next_data_array()) {
         my ($chr,$pos,$id,$ref,$alt,$QUAL,$FILTER,$INFO,$FORMAT,$sample) = 
            ($$x[0],$$x[1],$$x[2],$$x[3],$$x[4],$$x[5],$$x[6],$$x[7],$$x[8],$$x[9]);
            
        my  $index = $chr.':'.$pos.':'.$ref;
        
        
        my @sample = split(/\:/,$sample);
        #print $alt,"\n";
        #print $sample[0],"\n";
        my @geno  = split(/\/|\|/,$sample[0]);
        
        my $geno_site = $geno[0];
        
        if($paternal eq 'P'){
        	$geno_site = $geno[1];
        }

        my  @alt_original_old = split(/\,/,$alt);
        unshift @alt_original_old, $ref;  # put 0 values
         
        if (@alt_original_old > 1) {
         	    $hetero_num ++;
        }
        
         my  $first_alt  = $alt_original_old[$geno_site];## change as input
         
         my  $ref_len = length($ref);
         if($ref_len == 1){		    	
		      $vcf_hash{$index}{$first_alt}=1;
         } # ref eq 1
         else{ # ref length >1
	         	 my @ref_all = split(//,$ref);
	         	 if(@ref_all+0!= $ref_len ){
	    		      die "What's This Field";
	    	     }
			     my @alt_all = split(//,$first_alt);
			     for(my $i=0;$i<@ref_all;$i++){
			    			my $site_new = $pos+$i;
			    			my $ref_new  = $ref_all[$i];
			    			my $alt_out = '-';
			    			
			    			#### last record in ref
			    			if($i==(@ref_all-1)){
			    				if(@alt_all+0>0){
			    					$alt_out = join('',@alt_all); # all left
			    				}
			    			}
			    			else{
			    				if(@alt_all+0>0){
			    				   $alt_out = shift @alt_all;
			    			    }
			    			}
			    			my $index_new = $chr.':'.$site_new.':'.$ref_new;
			    		    $vcf_hash{$index_new}{$alt_out}=1;
			     }## for index
			   
        } ## ref lenght >1 
     }  ### read while  
     
     my %new_hash = %{&regenerate_hash(\%vcf_hash)};
     return \%new_hash ;
}# while




sub regenerate_hash {
	my $hash_addr = $_[0];
	
	my %vcf_hash = %{$hash_addr};
	my %vcf_site = ();
		for my $index (sort keys %vcf_hash){	
				my @hit = sort keys %{$vcf_hash{$index}};
				my @original = split(/\:/,$index);
				my $ref_nt = $original[-1];
				
				my @ooo = ();
				for( @hit ){
					 if($ref_nt ne $_){
					 	 push @ooo,$_;
					 }
				}
				@ooo  =  sort @ooo;
				
				if(@ooo > 1){
					 $conflict_num ++;
				}
				if(@ooo >0 ){
					my $out  =  $ooo[0];
					 $vcf_site{$index}=$out;
				}
			
	  }
	 return \%vcf_site;
}



sub re_co{
   my ($seq) = @_;
   $seq = reverse($seq);
   $seq =~ tr/ACGT/TGCA/;
   return($seq);
}


sub usage_message {
    my $version = shift;
    my $usage_message = "\n Trans extraction $version
Usage: perl get_allele_specific_transcripts.pl --gff  <gff>  --genomefile <genome>  --vcf <vcf>  --out <out>

";
    return $usage_message;
}

__END__

