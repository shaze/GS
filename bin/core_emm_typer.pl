#! /bin/env perl


my $emm_DB=$ARGV[0];
print "emm_DB=$emm_DB\n";
my $outName=$ARGV[1];
chomp($emm_DB);


my $emmType_output = "${outName}_emm_type_results.txt";
open(my $fh,'>',$emmType_output) or die "Could not open file '$emmType_output' $!";
print $fh "Sample_Name\temm_Type\temm_Seq\tPercent_Identity\tMatch_Length\n";


###Get the best blast hit by sorting the blast output by bit score, then % ID, then alignment length and select the first hit as the best match.###
my $frwd_bestHit = `cat contig-vs-frwd_nucl.txt | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
my @frwd_bestArray = split('\t',$frwd_bestHit);
my $best_frwd_name = $frwd_bestArray[0];
my $best_frwd_len = $frwd_bestArray[3];
my $best_frwd_iden = $frwd_bestArray[2];
my $query_strt = $frwd_bestArray[6];
my $query_end = $frwd_bestArray[7];
#my $target_strt = $frwd_bestArray[8];
#my $target_end = $frwd_bestArray[9];

print "\nname of best hit against the emm forward primers: $best_frwd_name\n";
print "% identity of best hit against the emm foward primers: $best_frwd_iden\n";
print "length of best hit against the emm forward primer: $best_frwd_len\n";

if ($best_frwd_len == 19 && $best_frwd_iden >= 94.5) {
    if ($frwd_bestArray[9] > $frwd_bestArray[8]) {
	##Extract the +500 sequence with Bedtools##
	my $query_extract = $query_strt + 500;
	open(my $fhe, '>', 'emm_region_extract.bed') or die "Could not open file 'emm_region_extract.bed' $!";
	print $fhe "$best_frwd_name\t$query_strt\t$query_extract\n";
	close $fhe;
	#print "done\n";
	system("bedtools getfasta -fi ./velvet_output/contigs.fa -bed emm_region_extract.bed -fo emm_region_extract.fasta");
	if ($? !=0) { die "Failed: bedtools getfasta -fi ./velvet_output/contigs.fa -bed emm_region_extract.bed -fo emm_region_extract.fasta" }
    } else {
	print "extract from reverse strand\n";
	##Extract the -500 sequence with Bedtools##
        my $query_extract = $query_strt - 500;
        open(my $fhe, '>', 'emm_region_extract.bed');
        print $fhe "$best_frwd_name\t$query_extract\t$query_strt\n";
        close $fhe;
	
	my $extract_emm1 = `bedtools getfasta -tab -fi ./velvet_output/contigs.fa -bed emm_region_extract.bed -fo stdout`;
        print "extract emm is:\n$extract_emm1\n";
	my @emm1_array = split('\t',$extract_emm1);
	my $rev_comp_emm1 = reverse($emm1_array[1]);
	$rev_comp_emm1 =~ tr/ATGCatgc/TACGtacg/;

        open(my $fh2, '>', 'emm_region_extract.fasta') or die "Could not open file 'emm_region_extract.fasta' $!";
        print $fh2 ">$emm1_array[0]";
	print $fh2 "$rev_comp_emm1\n";
        close $fh2;
    }
} else {
    print "The best blast hit ($best_frwd_name) obtained from querying the assembled contigs against the emm forward primers\ndidn't meet minimum criteria of length and identity to call a true match.\n";
    print $fh "$outName\tExtraction_Error\t--\t--\t--\n";
    my $old_name = "emm-Type_Results.txt";
    my $emm_out = $outName."__emm-Type__Results.txt";
    rename $old_name, $emm_out;
    exit
}

###Blast extracted emm region against database to find match###
if ( -s "emm_region_extract.fasta") {
    print "Doing blast blastn -db $emm_DB -query emm_region_extract.fasta -outfmt 6 -word_size 4 -out emm_vs_DB_nucl.txt\n\n";
  system("blastn -db $emm_DB -query emm_region_extract.fasta -outfmt 6 -word_size 4 -out emm_vs_DB_nucl.txt");
} else {
    print "Although the best blast hit found a true match against the 19bp primer, the matching contig didn't contain the\nfull 500bp region downstream of the primer that comprises the emm-typing region\n";
    print $fh "$outName\tExtraction_Error\t--\t--\t--\n";
    my $old_name = "emm-Type_Results.txt";
    my $emm_out = $outName."__emm-Type__Results.txt";
    rename $old_name, $emm_out;
    exit
}

my $emm_bestHit = `cat emm_vs_DB_nucl.txt | sort -k12,12 -nr -k3,3 -k4,4 | head -n 1`;
my @emm_bestArray = split('\t',$emm_bestHit);
my $best_emm_name = $emm_bestArray[1];
my $best_emm_len = $emm_bestArray[3];
my $best_emm_iden = $emm_bestArray[2];


$outName =~ /EMM_(.*)/;
my $finalName = $1;
if ($best_emm_iden == 100 && $best_emm_len == 180) {
    #$best_emm_name =~ /\d+__EMM(.*)__EMM.*__\d+/;
    $best_emm_name =~ /\d+__[A-Z]+(.*)__.*__\d+/;
    my $emmType = $1;
    print "\n\n$finalName\t$emmType\t$best_emm_name\t$best_emm_iden\t$best_emm_len\n";
    print $fh "$finalName\t$emmType\t$best_emm_name\t$best_emm_iden\t$best_emm_len\n";
} else {
    #$best_emm_name =~ /\d+__EMM(.*)__EMM.*__\d+/;
    $best_emm_name =~ /\d+__[A-Z]+(.*)__.*__\d+/;
    my $emmType = $1;
    print $fh "$finalName\t$emmType*\t$best_emm_name\t$best_emm_iden\t$best_emm_len\n";
    ###OPEN 'Check_Target_Sequence.txt' FOR APPENDING (CHECK FOR FAILURES)
    open ( my $exOUT, ">>", 'Check_Target_Sequence.txt' ) or die "Could not open file 'Check_Target_Sequence.txt': $!";
    ###OPEN 'emm_region_extract.fasta' for READING (CHECK FOR FAILURES)
    open ( my $newEx, "<", 'emm_region_extract.fasta' ) or die "Could not open file 'emm_region_extract.fasta': $!";
    ###READ EACH LINE OF FILE B.txt (BAR) and add it to FILE A.txt (FOO)
    #print $exOUT "##############################--NEW EMM TYPE SEQUENCE--##############################\n";
    print $exOUT '#' x 65;
    print $exOUT "--NEW EMM TYPE SEQUENCE--";
    print $exOUT '#' x 65;
    print $exOUT "\n\n";    
    while ( my $line = <$newEx> ) {
	print $exOUT $line;
    }
    #print $exOUT "#####################################################################################\n";    
    print $exOUT '#' x 150;
    print $exOUT "\n\n";
    close $exOUT;
}
close $fh;

my $old_name = "emm-Type_Results.txt";
my $emm_out = $outName."__emm-Type__Results.txt";
rename $old_name, $emm_out;
