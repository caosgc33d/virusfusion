#!/usr/bin/perl
use strict;
use warnings;

my $version = 1.0;

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
This script will run the alignment by bwa.
Pipeline version: $version
$yellow		Usage: perl $0 <run_folder> $normal
<run_folder> = full path of the folder holding files for this sequence run
Run bwa
OUT

#my $bwa_ref="";
my $bwa_ref = "/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37PlusSelectedVirus2015.04/GRCh37Plus28Virus.fa";

die $usage unless @ARGV == 1;
my ($run_dir) = @ARGV;

if ($run_dir =~/(.+)\/$/) {
	$run_dir = $1;
}

#my $bwa_ref="/gscuser/scao/gc3027/fasta/bacteria/ref.fa";
#####################################################################################
# values need to be modified to adapt to local environment
# software path
#my $bwa = "bwa";

#####################################################################################
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];

# To run jobs faster, split large fasta files to small ones. Split to specific number of 
# files instead of specific sequences in each small file, because the number of job array 
# cannot be determined if spliting to specific number of sequences in each file. Job 
# number is required by qsub ${SGE_TASK_ID}. The minimum size of each file is 4kb. 
# The number of files should be determined accourding to CPUs available in the computer
# cluster.

# The number of small fasta files to split to from a large file for RepeatMasker
#my $file_number_of_RepeatMasker = 200; #default 
# the number of small fasta files to split to from a large file for Blast_Reference_Genome
#my $file_number_of_Blast_Ref_Genome = 200; #default
# the number of small fasta files to split to from a large file for Blast_N
#my $file_number_of_Blast_N = 200; #default
# the number of small fasta files to split to from a large file for Blast_X
#my $file_number_of_Blast_X = 200; #default

#store job files here
if (! -d $HOME."/tmp") {
	`mkdir $HOME"/tmp"`;
}
my $job_files_dir = $HOME."/tmp";

#store SGE output and error files here
if (! -d $HOME."/SGE_DIR") {
	`mkdir $HOME"/SGE_DIR"`;
}
my $lsf_file_dir = $HOME."/SGE_DIR";

# obtain script path
my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

my $qsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
my $current_job_file; 
#directory suffix constants

# get sample list in the run, name should not contain "."
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure

# start data processsing
	#begin to process each sample
for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
		$sample_name = $sample_dir_list[$i];
		if (!($sample_name =~ /\./)) {
			$sample_full_path = $run_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
			     &bsub_bwa(); 	
				print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
				}
			}
		}

########################################################################
sub bsub_bwa{

	#my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

	$current_job_file = "bwa_".$sample_name.$$.".sh";
		
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
	print BWA "#!/bin/bash\n";
        print BWA "#BSUB -n 1\n";
        print BWA "#BSUB -R \"rusage[mem=20000]\"","\n";
        print BWA "#BSUB -M 20000000\n";
        print BWA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
        print BWA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
        print BWA "#BSUB -J $current_job_file\n";
	
	print BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
	print BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
 	print BWA "BWA_sai=".$sample_full_path."/".$sample_name.".sai\n";
	print BWA "BWA_sam=".$sample_full_path."/".$sample_name.".sam\n";
    print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".realign.bam\n";
	print BWA "BWA_mapped_bam=".$sample_full_path."/".$sample_name.".mapped.bam\n";
	print BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped\n";
    print BWA "VIRUSN=".$sample_full_path."/".$sample_name.".virus.sam\n";
	print BWA 'if [ ! -f $BWA_mapped ]',"\n";
    print BWA "    then\n";
	print BWA 'bedtools bamtofastq -i ${BWA_IN} -fq ${BWA_fq}',"\n";
	print BWA "bwa aln $bwa_ref \${BWA_fq} > \${BWA_sai}","\n";
	print BWA "bwa samse $bwa_ref \${BWA_sai} \${BWA_fq} > \${BWA_sam}","\n";
	print BWA "samtools view -b -S \${BWA_sam} > \${BWA_bam}","\n";
	print BWA "samtools view -F 4 -b -o \${BWA_mapped_bam} \${BWA_bam}","\n";
	print BWA "samtools view \${BWA_mapped_bam} > \${BWA_mapped}","\n";	
    print BWA "samtools view \${BWA_bam} | perl -ne \'\$l=\$_; \@ss=split(\"\\t\",\$l); if(\$ss[2]=~/\^gi/) { print \$l; }\' > \${VIRUSN}","\n";
    print BWA "   fi\n";
	close BWA;
    $qsub_com = "bsub < $job_files_dir/$current_job_file\n";
	system ( $qsub_com );
}

