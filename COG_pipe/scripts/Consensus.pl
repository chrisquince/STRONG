#!/usr/bin/perl

use strict;

my %hashContigAssign = ();

my @clusters = ();
my %contigs = ();

my $line = <STDIN>;
print "$line";

while(my $line = <STDIN>){
    chomp($line);

    my @tokens = split(/,/,$line);

    my $split = $tokens[0];

    my $assign = $tokens[1];

    if($split =~ /(.*)\.(.*)/){
        my $contig = $1;
        $contigs{$contig}++;
        $hashContigAssign{$contig}{$assign} += 1;

        $clusters[$assign] += 1;
    }
    else{
        $hashContigAssign{$split}{$assign} += 1;
        $contigs{$split}++;
        $clusters[$assign] += 1;
    }
}

my @newFreqs = ();
my %hashNewAssign = ();

foreach my $contig(keys %contigs){
    my %contigAssign = %{$hashContigAssign{$contig}};
    
    #print "$contig\n";
    my @sorted =  sort { $contigAssign{$b} <=> $contigAssign{$a} } keys %contigAssign;
    my $bestAssign = $sorted[0];

    $hashNewAssign{$contig} = $bestAssign;
#    print "$contig,$bestAssign\n";
    $newFreqs[$bestAssign] += 1;    

}

my %reMap = ();
my $newIdx = 0;
my $i = 0;
foreach my $freq(@newFreqs){
    if ($freq > 0){
        $reMap{$i} = $newIdx;
        $newIdx += 1;
    }
    else{
        #print "Remove $i\n";
    }
    $i += 1;
}

foreach my $contig(keys %hashNewAssign){
    print "$contig,$reMap{$hashNewAssign{$contig}}\n";
}



