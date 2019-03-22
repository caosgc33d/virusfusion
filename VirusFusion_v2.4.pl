#########Song Cao###########

#VirusFusion_V2.1.pl
#updated date: Nov 18 2014
#V2.1: the version which can deal with UNC bams and normal bams
#V2.2. use new virus database
#V2.3: add hh6 to virus reference
#V2.4: add HPV30 to virus reference

#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
my $version = 2.4;

#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m"; 

#usage information
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the virus discovery pipeline on Sun Grid Engine.
Pipeline version: $version
$yellow		Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$green		[1]  Run bwa realignment using paired-end information againt human+virus genome
$red        [2]  Run annotation for fusion site

$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
	$run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 2) || ($step_number >= 12));


#####################################################################################
# values need to be modified to adapt to local environment
my $email = "scao\@wustl\.edu";

# software path
#my $cd_hit = "/gscuser/mboolcha/software/cdhit/cd-hit-est";
#my $repeat_masker = "RepeatMasker";
#my $blastn = "/gscuser/scao/tools/ncbi-blast+/bin/blastn";
#my $blastx = "/gscuser/scao/tools/software/ncbi-blast+/bin/blastx";

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# path and name of databases
#my $db_BN = "/gscuser/scao/gc3027/nt/nt";
#my $db_BX = "/gscuser/scao/gc3027/nr/nr";
#my $bwa_ref = "/gscuser/scao/gc3027/fasta/virus/virusdb_082414.fa";

#my $db_BN = "/gscmnt/gc3027/dinglab/medseq/nt/nt";
#my $db_BX = "/gscmnt/gc3027/dinglab/medseq/nr/nr";
#my $bwa_ref = "/gscmnt/gc3027/dinglab/medseq/fasta/";

# reference genome taxonomy classification and database location.
# It's better to change $refrence_genome_taxonomy and $reference_genome based on the data being analyzed.
#my $refrence_genome_taxonomy = "";
#my $reference_genome = "";

my $bwa_ref = "/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37PlusSelectedVirus2015.04/GRCh37Plus28Virus.fa";

#####################################################################################
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];

my $HOME1="/gscmnt/gc2524/dinglab";

#store job files here
if (! -d $HOME1."/tmp") {
    `mkdir $HOME1"/tmp"`;
}
my $job_files_dir = $HOME1."/tmp";

#store SGE output and error files here
if (! -d $HOME1."/SGE_DIR") {
    `mkdir $HOME1"/SGE_DIR"`;
}
my $lsf_file_dir = $HOME1."/SGE_DIR";

# obtain script path
my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

my $hold_RM_job = "norm"; 
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
# get sample list in the run, name should not contain "."
my $all_sample_list="/gscuser/scao/gc2517/dinglab/virus_histogram_figures/pan23_all_set_virus_count.dat"; 
my $anno_bed_b37="/gscuser/scao/gc2524/dinglab/bed_maker/E75_gene_start_end_pos.tsv";
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
&check_input_dir($run_dir);

# start data processsing
if ($step_number < 14 || $step_number>=22) {
	#begin to process each sample
	for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_dir_list[$i];
		if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
			$sample_full_path = $run_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
				$current_job_file="";
				if ($step_number == 0 || $step_number>=12) {#run the whole pipeline
					######################################################################
					#bwa
                    if($step_number==0) 
					{ &bsub_bwa();}
					if($step_number<=12) {
					&bsub_anno();}
				
				######################################################################
				#run the pipeline step by step
				}elsif ($step_number == 1) {
					&bsub_bwa();
				}elsif($step_number == 2) {
                    &bsub_anno(1);
                }
			}
		}
	}
}


#######################################################################
# send email to notify the finish of the analysis
if (($step_number == 0) || ($step_number == 2) || ($step_number>=12)) {
	print $yellow, "Submitting the job for sending an email when the run finishes ",$sample_name, "...",$normal, "\n";
	$hold_job_file = $current_job_file;
	$current_job_file = "Email_run_".$$.".sh";
	open(EMAIL, ">$job_files_dir/$current_job_file") or die $!;
	print EMAIL "#!/bin/bash\n";
    print EMAIL "#BSUB -n 1\n";
    print EMAIL "#BSUB -o $lsf_file_dir","\n";
    print EMAIL "#BSUB -e $lsf_file_dir","\n";
    print EMAIL "#BSUB -J $current_job_file\n";
	print EMAIL "#BSUB -w \"$hold_job_file\"","\n";	
	print EMAIL $run_script_path."send_email.pl ".$run_dir." ".$email."\n";
	close EMAIL;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
	#$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
	system ($bsub_com);
}
#######################################################################
if ($step_number == 0) {
	print $green, "All jobs are submitted! You will get email notification when this run is completed.\n",$normal;
}

exit;


########################################################################
# subroutines 

sub check_input_dir {
	my ($input_dir) = @_;
	my $have_input_sample = 0;
	
	# get sample list in the run, name should not contain "."
	opendir(DH, $input_dir) or die "Cannot open dir $input_dir: $!\n";
	my @sample_list = readdir DH;
	close DH;
	
	for (my $i=0;$i<@sample_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_list[$i];
		if (!($sample_name =~ /\./)&&!($sample_name =~/Analysis_/) && !($sample_name =~/worklog/)) {
			$have_input_sample = 1;
			$sample_full_path = $input_dir."/".$sample_name;
			print $sample_full_path,"\n";
			if (-d $sample_full_path) { # is a full path directory containing a sample
				my $input_file = $input_dir."/".$sample_name."/".$sample_name.".bam";
				if (!(-e $input_file)) { # input file does not exist
					print $input_file,"\n";
					print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
					die;
				}
			}
			else { # input sample directory does not exist
				print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
				die;
			}
		}
	}

	if (!($have_input_sample)) { # does not have any input sample directory
		print $red, "Do not have appropriate input directory structure. Please check your command line argument!", $normal, "\n\n";
		die;
	}

}

########################################################################
sub bsub_bwa{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_bwa_sampe_".$sample_name.".sh";
    my $IN_bam = $sample_full_path."/".$sample_name.".bam";
    if (! -e $IN_bam) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam is empty!", $normal, "\n\n";
    }
    open(BWA, ">$job_files_dir/$current_job_file") or die $!;
    my $f_fq_1=$sample_full_path."/".$sample_name.".1.fq";
	my $f_fq_2=$sample_full_path."/".$sample_name.".2.fq";
	
	print BWA "#!/bin/bash\n";
    print BWA "#BSUB -n 1\n";
    print BWA "#BSUB -R \"rusage[mem=30000]\"","\n";
    print BWA "#BSUB -M 30000000\n";
    print BWA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print BWA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print BWA "#BSUB -J $current_job_file\n";
    print BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
    #print BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
  	print BWA "BWA_sorted=".$sample_full_path."/".$sample_name.".sorted\n";
	print BWA "BWA_sorted_bam=".$sample_full_path."/".$sample_name.".sorted.bam\n";
	print BWA "BWA_sai_1=".$sample_full_path."/".$sample_name.".1.sai\n";
	print BWA "BWA_sai_2=".$sample_full_path."/".$sample_name.".2.sai\n";
    print BWA "BWA_fq_1=".$sample_full_path."/".$sample_name.".1.fq\n";
    print BWA "BWA_fq_2=".$sample_full_path."/".$sample_name.".2.fq\n";	
    print BWA "BWA_sam=".$sample_full_path."/".$sample_name.".realign.sam\n";
    print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".realign.bam\n";
	print BWA "VIRUSN=".$sample_full_path."/".$sample_name.".virus.sam\n";
	print BWA "DISCORDANT=".$sample_full_path."/".$sample_name.".discordant.sam\n";
    #print BWA "BWA_mapped_bam=".$sample_full_path."/".$sample_name.".mapped.bam\n";
    #print BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped.reads\n";
    #print BWA "BWA_fa=".$sample_full_path."/".$sample_name.".fa\n";
	#print BWA 
	print BWA 'if [ ! -s $BWA_bam ]',"\n";
    print BWA "    then\n";
	print BWA "rm \${BWA_sai_1}","\n";
	print BWA "rm \${BWA_sai_2}","\n";
    print BWA "rm \${BWA_sam}","\n";
    print BWA "rm \${BWA_fq_1}","\n";
    print BWA "rm \${BWA_fq_2}","\n";
    print BWA "rm \${BWA_sorted_bam}","\n";
	print BWA "samtools sort -n \${BWA_IN} \${BWA_sorted}","\n";
	print BWA "samtools view \${BWA_sorted_bam} | perl -ne \'\$l=\$_; \$f_q1=\"$f_fq_1\"; \$f_q2=\"$f_fq_2\"; if(\$first==0) { open(OUT1,\">\$f_q1\"); open(OUT2,\">\$f_q2\");  \$first=1;}  \@ss=split(\"\\t\",\$l); \$flag=\$ss[1]; \$cigar=\$ss[5]; if((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) { next; } \$id=\$ss[0]; \$seq=\$ss[9]; \$q=\$ss[10];  if(\$id=~/\\/1\$/ || (\$flag & 0x40) ) { \$r1=\$id; \$r1=~s/\\/1\$//g; \$seq1=\$seq; \$q1=\$q; if((\$flag & 0x4) || (\$cigar=~/S/)) { \$unmapped1=1; } else { \$unmapped1=0; }  } if(\$id=~/\\/2\$/ || (\$flag & 0x80)) { \$r2=\$id; \$r2=~s/\\/2\$//g; \$seq2=\$seq; \$q2=\$q; if((\$flag & 0x4) || (\$cigar=~/S/)) { \$unmapped2=1; } else { \$unmapped2=0;} } if((\$r1 eq \$r2) && (\$unmapped1==1 || \$unmapped2==1)) { print OUT1 \"\@\",\$r1,\"/1\",\"\\n\"; print OUT1 \$seq1,\"\\n\"; print OUT1 \"+\",\"\\n\"; print OUT1 \$q1,\"\\n\"; print OUT2 \"\@\",\$r1,\"/2\",\"\\n\"; print OUT2 \$seq2,\"\\n\"; print OUT2 \"+\",\"\\n\"; print OUT2 \$q2,\"\\n\";}\'","\n";
	print BWA "bwa aln $bwa_ref \${BWA_fq_1} > \${BWA_sai_1}","\n";
    print BWA "bwa aln $bwa_ref \${BWA_fq_2} > \${BWA_sai_2}","\n";
	print BWA "bwa sampe $bwa_ref \${BWA_sai_1} \${BWA_sai_2} \${BWA_fq_1} \${BWA_fq_2} > \${BWA_sam}","\n";
	print BWA "samtools view -b -S \${BWA_sam} > \${BWA_bam}","\n";
	print BWA "samtools view \${BWA_bam} | perl -ne \'\$l=\$_; \@ss=split(\"\\t\",\$l); if(\$ss[2]=~/\^gi/) { print \$l; }\' > \${VIRUSN}","\n";
	print BWA "samtools view \${BWA_bam} | perl -ne \'\$l=\$_; \@ss=split(\"\\t\",\$l); if(\$ss[2]=~/\^gi/  && (\!(\$ss[6]=~/\^gi/))  && (!(\$ss[6]=~/\\=/))) { print \$l; }\' > \${DISCORDANT}","\n";
	print BWA "rm \${BWA_fq_1}","\n";
	print BWA "rm \${BWA_fq_2}","\n";
	print BWA "rm \${BWA_sai_1}","\n";
    print BWA "rm \${BWA_sai_2}","\n";
	print BWA "rm \${BWA_sam}","\n";
	print BWA "rm \${BWA_sorted_bam}","\n";
	print BWA "fi\n";
    close BWA;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}

########################################################################
sub bsub_anno{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j2_".$sample_name."_anno.sh";
	my $ANNO_BED=$sample_full_path."/".$sample_name.".discordant.bed";
	my $ANNO_SUMMARY=$sample_full_path."/".$sample_name.".discordant.summary\n";
    open(ANNO, ">$job_files_dir/$current_job_file") or die $!;
    print ANNO "#!/bin/bash\n";
    print ANNO "#BSUB -n 1\n";
    #print ANNO "#BSUB -q ding-lab\n";
    print ANNO "#BSUB -R \"rusage[mem=10000]\"","\n";
    print ANNO "#BSUB -M 10000000\n";
    print ANNO "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print ANNO "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print ANNO "#BSUB -J $current_job_file\n";
    print ANNO "#BSUB -w \"$hold_job_file\"","\n";
    print ANNO "ANNO_IN=".$sample_full_path."/".$sample_name.".discordant.sam\n";
    print ANNO 'if [ -s $ANNO_IN ]',"\n";
    print ANNO "    then\n";
	print ANNO "awk \'{print \$7,\$8,\$4}\' \${ANNO_IN} | sort -V -k1,1 -k2,2 > $ANNO_BED","\n";
	print ANNO "perl -e \'\$f_g=\"$ANNO_BED\"; \$f_r=\"$anno_bed_b37\"; \$sn=\$f_g; \$sn=\"$sample_name\"; foreach \$a (\`cat \$f_r\`) { \@ss=split(\"\\t\",\$a); \$chr=\$ss[0]; \$l=\$ss[1]; \$r=\$ss[2]; \$g=\$ss[3]; chomp(\$g); \$gbed{\$g}{\$chr}{\$l}{\$r}=1; } foreach \$a (\`cat \$f_g\`) { \@ss=split(\" \",\$a); \$chr=\$ss[0]; \$pos=\$ss[1]; \$scount=0; foreach \$g (keys \%gbed) { foreach \$c (keys \%{\$gbed{\$g}}) { if(\$chr==\$c) { foreach \$l (keys \%{\$gbed{\$g}{\$c}}) { foreach \$r (keys \%{\$gbed{\$g}{\$c}{\$l}}) {  if(\$r>=\$pos && \$l<=\$pos && \$scount==0) { \$h{\$g}++; \$scount=1; }}}}}}} foreach \$l (\`cat $all_sample_list\`) { \@ss=split(\"\\t\",\$l); \$s_2_t{\$ss[0]}=\$ss[2]; \$s_2_nt{\$ss[0]}=\$ss[1]; } foreach \$g (sort { \$h{\$b} <=> \$h{\$a} } keys %h) { print \$sn,\"\\t\",\$s_2_nt{\$sn},\"\\t\",\$s_2_t{\$sn},\"\\t\",\$g,\"\\t\",\$h{\$g},\"\\n\";}\'> $ANNO_SUMMARY","\n"; 
    print ANNO "fi\n";
    close ANNO;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

}
