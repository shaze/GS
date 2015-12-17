#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

###MODULE LOAD###
#module load perl/5.12.3
#module load ncbi-blast+/2.2.29
#module load BEDTools/2.17.0
#module load Python/2.7

sub checkOptions {
    my %opts;
    getopts('h1:2:r:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $sero_DB, $outDir, $outName);

    if($opts{h}) {
        $help = $opts{h};
        help();
    }

    if($opts{1}) {
        $fastq1 = $opts{1};
        if( -e $fastq1) {
            print "Paired-end Read 1 is: $fastq1\n";
        } else {
            print "The forward paired-end file name is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "No paired end 1 fastq file path argument given.\n";
        help();
    }

    if($opts{2}) {
        $fastq2 = $opts{2};
        if( -e $fastq2) {
            print "Paired-end Read 2 is: $fastq2\n";
        } else {
            print "The reverse paired-end file name is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "No paired end 2 fastq file path argument given.\n";
        help();
    }

    if($opts{r}) {
        $sero_DB = $opts{r};
        if (-e $sero_DB) {
            print "The serotype reference database sequence: $sero_DB\n";
        } else {
            print "The serotype reference sequence location is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The serotype reference sequence location (including full path) has not been given.\n";
        help();
    }

    $outDir = "./";
    if($opts{o}) {
        if (-d $opts{o}) {
            $outDir = $opts{o};
            print "The output directory is: $outDir\n";
        } else {
            $outDir = $opts{o};
            mkdir $outDir;
            print "The output directory has been created: $outDir\n";
        }
    } else {
        print "The files will be output into the current directory.\n";
    }

    if($opts{n}) {
        $outName = $opts{n};
        print "The output file name prefix: $outName\n";
    } else {
        $outName = `echo "$fastq1" | awk -F"/" '{print \$(NF)}' | sed 's/_S[0-9]\\+_L[0-9]\\+_R[0-9]\\+.*//g'`;
        print "The default output file name prefix is: $outName";
    }

    return ($help, $fastq1, $fastq2, $sero_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
emm_typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <reference databases directory: file path> -o <output directory name: string> -n <output name prefix: string>  [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   serotype reference sequence directory (including full path)
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $sero_DB, $outDir, $outName) = checkOptions( @ARGV );


##Start Doing Stuff##
chdir "$outDir";
my $serotype_output = "TEMP_SeroType_Results.txt";
open(my $fh,'>',$serotype_output) or die "Could not open file '$serotype_output' $!";
print $fh "Matched_Allele\tMatch_Type\tSerotype\tAvgDepth\n";

###Detect GAS serotype sequence###
#my $sero_outname = "SERO_$outName";
system("mod-srst2.py --input_pe $fastq1 $fastq2 --output $outName --log --save_scores --min_coverage 99 --max_divergence 4 --gene_db $sero_DB");

###mpileup the 'SERO_.*.sorted.bam and create the called variants file with freebayes.
opendir(DIR, ".") or die "Couldn't open directory for reading: $!";
my @seroTargets = grep (/SERO_.*__fullgenes__.*\.txt/,readdir(DIR));
closedir(DIR);
my @seroT_bam = glob("SERO_*\.sorted\.bam");

my $seroT_output = $seroTargets[0];
#print "\nseroT_output is: $seroT_output\n";
open(MYINPUTFILE, "$seroT_output");
my %srst2_seroT;
while(<MYINPUTFILE>) {
    my $line = $_;
    chomp($line);
    my @seroT_line = split('\t', $line);
    if ($seroT_line[2] eq "gene") {
	next;
    } else {
	if ($seroT_line[5] > 10) {
	    my $newLine = "$seroT_line[2]:$seroT_line[3]:$seroT_line[4]:$seroT_line[5]:$seroT_line[6]:$seroT_line[7]";
	    $srst2_seroT{$seroT_line[2]} = $newLine;
	}
    }
}
#print Dumper(\%srst2_seroT);

if (exists $srst2_seroT{III} && (! exists $srst2_seroT{II})) {
    #print "Serotype 3: $srst2_seroT{III}\n";
    my @seroT_outLine = split(':',$srst2_seroT{III});
    my $status = "identical";
    if (exists $seroT_outLine[4]) {
	$status = "imperfect";
    }
    print $fh "$seroT_outLine[1]\t$seroT_outLine[0]=$status\t$seroT_outLine[0]\t$seroT_outLine[3]\n";
    delete $srst2_seroT{III};
} elsif (exists $srst2_seroT{II}) {
    #print "Serotype 2: $srst2_seroT{II}\n";
    my @seroT_outLineII = split(':',$srst2_seroT{II});
    my @seroT_outLineIII = split(':',$srst2_seroT{III});
    my $status = "identical";
    if (exists $seroT_outLineII[4]) {
        $status = "imperfect";
    }
    print $fh "$seroT_outLineII[1]\t$seroT_outLineII[0]=$status\t$seroT_outLineII[0]\t$seroT_outLineII[3]\n";
    delete $srst2_seroT{II};
    delete $srst2_seroT{III};
}

foreach my $key (keys(%srst2_seroT)) {
    #print "None II or III serotype: ";
    #print "$srst2_seroT{$key}\n";
    my @seroT_outLine = split(':',$srst2_seroT{$key});
    my $status = "identical";
    if (exists $seroT_outLine[4]) {
        $status = "imperfect";
    }
    print $fh "$seroT_outLine[1]\t$seroT_outLine[0]=$status\t$seroT_outLine[0]\t$seroT_outLine[3]\n";
}

###Delete TEMP files### 
close $fh;

##Old Stuff
#(my $seroT_vcf = $seroT_bam) =~ s/\.bam/\.vcf/g;
#(my $seroT_bai = $seroT_bam) =~ s/\.bam/\.bai/g;
#print "@seroTargets || @seroT_bam || $seroT_bam || $seroT_vcf || $seroT_bai\n";

#my $seroT_output = "SERO_test_20153056__fullgenes__GBS_seroT_Gene-DB_Final__results.txt";
#my $seroT_output = "SERO_test_20153056__fullgenes__GBS_seroT_Gene-DB_Final__results_2.txt";
#my $seroT_output = "SERO_test_20153056__fullgenes__GBS_seroT_Gene-DB_Final__results_3.txt";
