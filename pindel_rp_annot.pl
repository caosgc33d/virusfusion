#!/usr/bin/perl
my $usage="pindel_rp_annot.pl input_file\n";

die $usage unless @ARGV == 1;

my ($f_rp)=@ARGV; 

my $all_sample_list="/gscuser/scao/gc2517/dinglab/virus_histogram_figures/pan23_all_set_virus_count.dat";
my $f_bed="/gscuser/scao/gc2524/dinglab/bed_maker/E75_gene_start_end_pos.tsv";

my $sn=(split(/\//,$f_rp))[-1];
$sn=~s/\_RP//g;

my %h=();
my %gbed=(); 
my %s_2_t=();
my %s_2_nt=();
my $f_out=$f_rp.".summary"; 

foreach my $a (`cat $f_bed`) 
{ 
	my @ss=split("\t",$a); 
    my $chr=$ss[0]; 
	my $l=$ss[1]; 
	my $r=$ss[2]; 
	my $g=$ss[3]; 
	chomp($g); 
	$gbed{$g}{$chr}{$l}{$r}=1; 
} 

foreach my $a (`cat $f_rp`) 
{ 
	my @ss=split("\t", $a); 
   	my $chr=$ss[5]; 
	my $pos_l=$ss[6]+1; 
	my $pos_r=$ss[7]+1;
	my $support=$ss[11];
	chomp($support);
	$support=~s/Support: //g;
	if($chr=~/gi/) { next; }
 	#print $chr,"\t",$pos_l,"\t",$pos_r,"\t",$support,"\n";
	#<STDIN>;
	$scount=0; 
	foreach my $g (keys %gbed) 
	{ 
		foreach my $c (keys %{$gbed{$g}}) 
		{
		#print $g,"\t",$c,"\n";
		 
		if($chr==$c) 
		{ 
		foreach $l (keys %{$gbed{$g}{$c}}) 
		{ foreach $r (keys %{$gbed{$g}{$c}{$l}}) 
		{  if(!($pos_r<$l || $pos_l>$r) && $scount==0) 
			{ $h{$g}+=$support; $scount=1; }}}}}}
} 

open(OUT,">$f_out");

foreach my $l (`cat $all_sample_list`) 
{ 
	my @ss=split("\t",$l); 
	$s_2_t{$ss[0]}=$ss[2]; 
	$s_2_nt{$ss[0]}=$ss[1]; 
} 

foreach my $g (sort { $h{$b} <=> $h{$a} } keys %h) 
{ print OUT $sn,"\t",$s_2_nt{$sn},"\t",$s_2_t{$sn},"\t",$g,"\t",$h{$g},"\n";}

close OUT;

