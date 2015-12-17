#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use File::Basename;
use File::Spec;

###MODULE LOAD###
#module load perl/5.12.3
#module load ncbi-blast+/2.2.29
#module load BEDTools/2.17.0
#module load Python/2.7

sub checkOptions {
    my %opts;
    getopts('h1:2:r:d:o:n:', \%opts);
    my ($help, $fastq1, $fastq2, $PBP_DB, $outDir, $outName);

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
        $PBP_DB = $opts{r};
        if (-e $PBP_DB) {
            print "The PBP reference database sequence: $PBP_DB\n";
        } else {
            print "The PBP reference sequence location is not in the correct format or doesn't exist.\n";
            print "Make sure you provide the full path (/root/path/fastq_file).\n";
            help();
        }
    } else {
        print "The PBP reference sequence location (including full path) has not been given.\n";
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

    return ($help, $fastq1, $fastq2, $PBP_DB, $outDir, $outName);
}

sub help
{

die <<EOF

USAGE
GBS_PBP-Gene_Typer.pl -1 <forward fastq file: fastq> -2 <reverse fastq file: fastq> -r <reference databases directory: file path> -o <output directory name: string> -n <output name prefix: string>  [OPTIONS]

    -h   print usage
    -1   forward fastq sequence filename (including full path)
    -2   reverse fastq sequence filename (including full path)
    -r   PBP reference sequence file path (including full path)
    -o   output directory
    -n   output name prefix

EOF
}

my ($help, $fastq1, $fastq2, $PBP_DB, $outDir, $outName) = checkOptions( @ARGV );





###SUBROUTINES###
sub PBP_blastTyper {
    my ($pbp_type,$pbp_name,$pbp_extract) = @_;
    chomp($pbp_type);
    chomp($pbp_name);
    chomp($pbp_extract);
    my $pbp_out;
    #print "type: $pbp_type || name: $pbp_name || LoTrac Extract: $pbp_extract\n";

    my $query_seq = extractFastaByID($pbp_name,$pbp_extract);
    #my $query_length = fasta_seq_length($query_seq);
    open ( my $qOUT, ">", 'TEMP_query_sequence.fna' ) or die "Could not open file TEMP_query_sequence.fna: $!";
    print $qOUT $query_seq;
    close $qOUT;
    ##Translate TEMP_query_sequence.fna to .faa file
    system("Translate_DNA-6Frame.pl -s TEMP_query_sequence.fna -f 1 -l 3 > TEMP_query_sequence.faa");
    my $query_length = fasta_seq_length("TEMP_query_sequence.faa");
    #print "PBP extracted seq length: $query_length\n";

    my $GBS_db_path = dirname($PBP_DB);
    my $GBS_blastDB_name = "Blast_bLactam_".$pbp_type."_prot_DB";
    my $GBS_blast_seq = "GBS_bLactam_".$pbp_type."-DB.faa";
    my $GBS_blast_out = "TEMP_".$outName."_blast-out_".$pbp_type."txt";
    print "Blast DB name: $GBS_db_path/$GBS_blastDB_name\n";
    if (!(glob("$GBS_db_path/$GBS_blastDB_name*"))) {
	print "Need to make a new Blast database\n";
	system("makeblastdb -in $GBS_db_path/$GBS_blast_seq -dbtype prot -out $GBS_db_path/$GBS_blastDB_name");
	system("blastp -db $GBS_db_path/$GBS_blastDB_name -query TEMP_query_sequence.faa -outfmt 6 -out $GBS_blast_out");
    } else {
	print "Blast database has already been created\n";
        system("blastp -db $GBS_db_path/$GBS_blastDB_name -query TEMP_query_sequence.faa -outfmt 6 -out $GBS_blast_out");
    }

    my $bestHit = `cat $GBS_blast_out | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
    print "best hit info: $bestHit";
    my @bestArray = split('\t',$bestHit);
    my $best_name = $bestArray[1];
    my $best_len = $bestArray[3];
    my $best_iden = $bestArray[2];

    print "name of best hit in the PBP database: $bestArray[1]\n";
    print "identity of best hit in the PBP database: $bestArray[2]\n";
    print "length of best hit in the PBP database: $bestArray[3]\n";

    if ($best_iden == 100 && $best_len == $query_length) {
	print "Found a match\n\n";
	($pbp_out = $best_name) =~ s/^([0-9]+)\|\|.*/$1/g;
    } else {
	if (-s "TEMP_query_sequence.fna") {
	    print "Didn't find match.  The sequence needs to be added to the database using 'bLactam-PBP_Updater.sh'\n\n";
	    $pbp_out = "NEW";
	    my $new_PBPseq = "PBP_".$pbp_type."_query_sequence.faa";
	    rename("TEMP_query_sequence.faa",$new_PBPseq);
	    my $fragPath = File::Spec->rel2abs("$new_PBPseq");
	    open(my $f_new,'>>',"TEMP_newPBP_allele_info.txt") or die "Could not open file 'TEMP_newPBP_allele_info.txt' $!";
	    print $f_new "$outName\t$fragPath\t$pbp_type\n";
	    close $f_new;
	} else {
	    $pbp_out = "ERROR";
	}
    }
return $pbp_out;
}

sub extractFastaByID {
    my ($lookup, $reference) = @_;
    open my $fh, "<", $reference or die $!;
    #print "lookup: $lookup\n";
    local $/ = "\n>";  # read by FASTA record

    my $output;
    while (my $seq = <$fh>) {
        chomp $seq;
        #print "while seq:\n$seq\n";
        my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        if ($id eq $lookup) {
            $seq =~ s/^>*.+\n//;  # remove FASTA header
            #$seq =~ s/\n//g;  # remove endlines
            #print ">$id\n";
            #print "$seq\n";
            $output = ">$id\n$seq\n";
            last;
        }
    }
    return $output;
}

sub fasta_seq_length {
    my ($seq) = @_;
    open ( my $q_seq, "<", $seq ) or die "Could not open file '$seq': $!";
    my @lines = split /\n/, $q_seq;
    my $final_line;
    #foreach my $line (@lines) {
    while (my $line = <$q_seq>) {
        chomp($line);
        if ($line =~ /^>/) {
            next;
        } else {
            #print "line: $line\n";
            $final_line = $final_line.$line;
        }
    }
    chomp($final_line);
    #print "final line: $final_line\n";
    return length($final_line);
}




###Start Doing Stuff###
print "\n";
chdir "$outDir";
##my $PBP_output = "PBP_".$outName."_Results.txt";
my $PBP_output = "TEMP_pbpID_Results.txt";
open(my $fh,'>',$PBP_output) or die "Could not open file '$PBP_output' $!";
##print $fh "Sample_Name\tPBP_1A\tPBP_2X\tPBP_Code\n";
print $fh "Sample_Name\tPBP_Code(1A:2X)\n";

system("perl ~/TEMP_GBS-Typing/LoTrac_target.pl -1 $fastq1 -2 $fastq2 -q $PBP_DB -S 2.2M -f -n $outName");

my $pbp_1A = glob("EXTRACT_*1A*.fasta");
my $pbp1A_fragName = `cat $pbp_1A | grep ">" | tail -n1 | sed 's/>//g'`;
my $pbp_2X = glob("EXTRACT_*2X*.fasta");
my $pbp2X_fragName = `cat $pbp_2X | grep ">" | tail -n1 | sed 's/>//g'`;

my $code_1A = PBP_blastTyper("1A",$pbp1A_fragName,$pbp_1A);
my $code_2X = PBP_blastTyper("2X",$pbp2X_fragName,$pbp_2X);
##print $fh "$outName\t$code_1A\t$code_2X\t$code_1A:$code_2X\n";
print $fh "$outName\t$code_1A:$code_2X\n";
close $fh;


###Delete Temp Files###
#unlink(TEMP*);
