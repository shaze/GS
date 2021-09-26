#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use Getopt::Long;
use Getopt::Std;

###MODULE LOAD###
#module load samtools/0.1.18
#module load bowtie2/2.1.0
#module load Python/2.7
#module load freebayes/0.9.21

sub checkOptions {
    my %opts;
    getopts('h1:2:r:p:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $ref_dir, $protein_DB, $outDir, $outName);

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
        $ref_dir = $opts{r};
        $ref_dir =~ s/\/$//g;
        if (-d $ref_dir) {
            print "Directory containing the surface and secretory protein reference sequences: $ref_dir\n";
        } else {
            print "The directory containing the surface and secretory protein references is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The surface protein reference sequence directory (including full path) has not been given.\n";
        help();
    }

    if($opts{p}) {
        $protein_DB = "$ref_dir/$opts{p}";
        if ($protein_DB) {
            print "The surface and secretory protein reference sequence file: $opts{p}\n";
        } else {
            print "The surface and secretory protein reference sequence file is not in the correct format or doesn't exist.\n";
            #print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The surface and secretory protein reference sequence file has not been given.\n";
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
        $outName = `echo "$fastq1" | awk -F"/" '{print $(NF)}' | sed 's/_S[0-9]\\+_L[0-9]\\+_R[0-9]\\+.*//g'`;
        print "The default output file name prefix is: $outName";
    }

    return ($help, $fastq1, $fastq2, $ref_dir, $protein_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
GBS_Surface_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <ref directory: dir> -m <misc. resistance seq: file> -v <vanc. resistance seq: file> -o <output directory name: string> -n <output name prefix: string> [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   reference sequence directory (including full path)
    -p   surface and secretory protein reference sequence file
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $ref_dir, $protein_DB, $outDir, $outName) = checkOptions( @ARGV );


###Subroutines###


###Start Doing Stuff###
chdir "$outDir";
my $surface_output = "TEMP_Surface_Results.txt";
open(my $fh,'>',$surface_output) or die "Could not open file '$surface_output' $!";
my $BIN_surf_out = "BIN_Surface_Results.txt";
open(my $bh,'>',$BIN_surf_out) or die "Could not open file '$BIN_surf_out' $!";
my @Bin_Feat_arr = (0) x 11;
#print $fh "Feature_Group\tTarget\n";
my $outNameFEAT = "SURF_".$outName;
system("srst2 --samtools_args '\\-A' --input_pe $fastq1 $fastq2 --output $outNameFEAT --log --save_scores --min_coverage 99.0 --max_divergence 8 --gene_db $protein_DB");
my @TEMP_SURF_fullgene = glob("SURF_*__fullgenes__*__results\.txt");
my $SURF_full_name = $TEMP_SURF_fullgene[0];

my %Feat_Col = (
    "ALPH" => "neg",
    "SRR" => "neg",
    "PILI" => "neg",
    "HVGA" => "neg",
    );

###Type the Presence/Absence Targets###
###############################################################################################
open(MYINPUTFILE, "$SURF_full_name");
while(<MYINPUTFILE>) {
    next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @feat_fullgene;
    @feat_fullgene = split('\t',$line);
    if ($feat_fullgene[5] >= 10) {
        if ($feat_fullgene[3] =~ m/(ALP|RIB)/) {
            if ($Feat_Col{"ALPH"} eq "neg") {
                $Feat_Col{"ALPH"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"ALPH"}.":".$feat_fullgene[2];
                $Feat_Col{"ALPH"} = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/SRR/) {
            if ($Feat_Col{"SRR"} eq "neg") {
                $Feat_Col{"SRR"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"SRR"}.":".$feat_fullgene[2];
                $Feat_Col{"SRR"} = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/PI/) {
            if ($Feat_Col{"PILI"} eq "neg") {
                $Feat_Col{"PILI"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"PILI"}.":".$feat_fullgene[2];
                $Feat_Col{"PILI"} = $new_val;
            }
        }
        if ($feat_fullgene[3] =~ m/HVGA/) {
            if ($Feat_Col{"HVGA"} eq "neg") {
                $Feat_Col{"HVGA"} = $feat_fullgene[2];
            } else {
                my $new_val = $Feat_Col{"HVGA"}.":".$feat_fullgene[2];
                $Feat_Col{"HVGA"} = $new_val;
            }
        }
    }
}
###############################################################################################

###############################################################################################
###Make Binary Output Table###
open(MYINPUTFILE, "$SURF_full_name");
my $dummy=<MYINPUTFILE>;
while(<MYINPUTFILE>) {
    #next if $. < 2;
    my $line = $_;
    chomp($line);
    #print "$line\n";
    my @feat_fullgene;
    @feat_fullgene = split('\t',$line);
    if ($feat_fullgene[5] >= 10 && $feat_fullgene[8] <= 8) {
        if ($feat_fullgene[3] =~ m/HVGA/) {
            $Bin_Feat_arr[0] = 1;
        }
        if ($feat_fullgene[3] =~ m/PI1/) {
            $Bin_Feat_arr[1] = 1;
        }
        if ($feat_fullgene[3] =~ m/PI2A1/) {
            $Bin_Feat_arr[2] = 1
        }
        if ($feat_fullgene[3] =~ m/PI2A2/) {
            $Bin_Feat_arr[3] = 1;
        }
        if ($feat_fullgene[3] =~ m/PI2B/) {
            $Bin_Feat_arr[4] = 1;
        }
        if ($feat_fullgene[3] =~ m/SRR1/) {
            $Bin_Feat_arr[5] = 1;
        }
        if ($feat_fullgene[3] =~ m/SRR2/) {
            $Bin_Feat_arr[6] = 1;
        }
        if ($feat_fullgene[3] =~ m/ALP1/) {
            $Bin_Feat_arr[7] = 1;
        }
        if ($feat_fullgene[3] =~ m/ALP23/) {
            $Bin_Feat_arr[8] = 1;
        }
        if ($feat_fullgene[3] =~ m/ALPHA/) {
            $Bin_Feat_arr[9] = 1;
        }
        if ($feat_fullgene[3] =~ m/RIB/) {
            $Bin_Feat_arr[10] = 1;
        }
    }
}
###Print GAS Binary Output###
print $bh join(',',@Bin_Feat_arr);

###############################################################################################
###Print GAS Features Output###
while (my ($key, $val) = each %Feat_Col) {
    my @val_arr = split(':',$val);
    #print "@val_arr\n";
    my @val_sort = sort { "\L$a" cmp "\L$b" } @val_arr;
    #print "@val_sort\n";
    my $val_out = join(':',@val_sort);
    print "$key\t$val_out\n";
    $Feat_Col{"$key"} = $val_out;
    #print $fh "$key\t$val_out\n";
}

print $fh "ALPH\t$Feat_Col{'ALPH'}\n";
print $fh "SRR\t$Feat_Col{'SRR'}\n";
print $fh "PILI\t$Feat_Col{'PILI'}\n";
print $fh "HVGA\t$Feat_Col{'HVGA'}\n";
###############################################################################################


