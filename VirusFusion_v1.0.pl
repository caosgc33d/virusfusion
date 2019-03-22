#########Song Cao###########

#VirusFusion_V1.0.pl.pl
#updated date: Nov 10 2014
 
#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
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
This script will run the virus discovery pipeline on Sun Grid Engine.
Pipeline version: $version
$yellow		Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$green		[1]  Run bwa realignment using paired-end information
$red		[2, or <=12]  Split files for RepeatMasker

$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
	$run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 17) || ($step_number >= 22));


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

my $bwa_ref = "/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37PlusSelectedVirus2014.11/GRCh37Plus19Virus.fa";

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
		if (!($sample_name =~ /\./)) {
			$sample_full_path = $run_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
				$current_job_file="";
				if ($step_number == 0 || $step_number>=12) {#run the whole pipeline
					######################################################################
					#bwa
                    if($step_number==0) 
					{ &bsub_bwa();}
					#if($step_number<=12) {
					#&summary_for_each_sample();}
				
				######################################################################
				#run the pipeline step by step
				}elsif ($step_number == 1) {
					&bsub_bwa();
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
		if (!($sample_name =~ /\./)&&!($sample_name =~/Analysis_/)) {
			$have_input_sample = 1;
			$sample_full_path = $input_dir."/".$sample_name;
			if (-d $sample_full_path) { # is a full path directory containing a sample
				my $input_file = $input_dir."/".$sample_name."/".$sample_name.".bam";
				if (!(-e $input_file)) { # input file does not exist
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
########################################################################
sub bsub_bwa{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_bwa_realign_".$sample_name.$$.".sh";
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
    print BWA "#BSUB -R \"rusage[mem=30000]\"","\n";
    print BWA "#BSUB -M 30000000\n";
    print BWA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print BWA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print BWA "#BSUB -J $current_job_file\n";
    print BWA "BWA_IN=".$sample_full_path."/".$sample_name.".bam\n";
    #print BWA "BWA_fq=".$sample_full_path."/".$sample_name.".fq\n";
    print BWA "BWA_sai_1=".$sample_full_path."/".$sample_name.".1.sai\n";
	print BWA "BWA_sai_2=".$sample_full_path."/".$sample_name.".2.sai\n";
    print BWA "BWA_sam=".$sample_full_path."/".$sample_name.".realign.sam\n";
    print BWA "BWA_bam=".$sample_full_path."/".$sample_name.".realign.bam\n";
    #print BWA "BWA_mapped_bam=".$sample_full_path."/".$sample_name.".mapped.bam\n";
    #print BWA "BWA_mapped=".$sample_full_path."/".$sample_name.".mapped.reads\n";
    #print BWA "BWA_fa=".$sample_full_path."/".$sample_name.".fa\n";
	#print BWA 
	print BWA 'if [ ! -s $BWA_bam ]',"\n";
    print BWA "    then\n";
	print BWA "rm \${BWA_sai_1}","\n";
	print BWA "rm \${BWA_sai_2}","\n";
    print BWA "rm \${BWA_sam}","\n";
	print BWA "bwa aln $bwa_ref -b1 \${BWA_bam} > \${BWA_sai_1}","\n";
    print BWA "bwa aln $bwa_ref -b2 \${BWA_bam} > \${BWA_sai_2}","\n";
    print BWA "bwa sampe $bwa_ref \${BWA_sai_1} \${BWA_sai_2} \${BWA_IN} \${BWA_IN} > \${BWA_sam}","\n";
	print BWA "samtools view -b -S \${BWA_sam} > \${BWA_bam}","\n";
	print BWA "rm \${BWA_sam}","\n";
	print BWA "fi\n";
    close BWA;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}

