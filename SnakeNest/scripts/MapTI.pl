#!/usr/bin/perl

my %mapDir = ();

open(FILE,$ARGV[0]) or die;

while($line = <FILE>){
    chomp($line);

    my @tokens = split(/,/,$line);

    $mapDir{$tokens[0]} = $tokens[1];
}

while($line = <STDIN>){
    chomp($line);

    if($line=~/>TI(\d+)/){
        $species = $mapDir{$1};
        print ">$species\n";
    }
    else{
        print "$line\n";
    }
}

