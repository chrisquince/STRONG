#!/usr/bin/perl

use strict;

my $nameFile = shift(@ARGV);

my @names = ();

open(FILE, $nameFile) or die "Can't open $nameFile\n";

while(my $line = <FILE>){
    chomp($line);
    my $name = $line;

    push(@names,$name);
}

close(FILE);

my @lengths = ();
my @cogs = ();
my %hashCogName = {};

foreach my $file(@ARGV){
    my @Seq = ();
    my @id       = ();

    my $count = 0;

    my $seq = "";

    open(FFILE, $file) or die "Can't open $file\n";
    while(my $line = <FFILE>){ 
        chomp($line);
    
        if($line =~ />(.*)_COG.*/){
#            print "$1\n";    
            $id[$count] = $1;
    
            if($seq ne ""){
                $Seq[$count - 1] = $seq;

                $seq = "";
            }

            $count++;
        }
        elsif($line =~ />(Cluster.*)/){
            $id[$count] = $1;
    
            if($seq ne ""){
                $Seq[$count - 1] = $seq;

                $seq = "";
            }

            $count++;
        }
        else{
            $seq .= $line;
        }
    }
    close(FFILE);
    $Seq[$count - 1] = $seq;
    my $total = $count;

    my $cog = "";
   if($file =~/.*\/(COG\d+)_all_al.gfa/){
        $cog = $1;
        push(@cogs,$cog);
        
        for(my $i = 0; $i < $total; $i++){
            my $seqlength = length($Seq[$i]);
            
            if($i == 0){
                push(@lengths,$seqlength);
            }
            my $name = $id[$i];
            
            #print "$name $cog $Seq[$i]\n";
            #print "$cog $id[$i] $seqlength\n";
            $hashCogName{$name}{$cog} = $Seq[$i];
            
        }
    }


}

my @combinedSeqs = ();

my $c = 0;

my @mapped = ();
foreach my $cog(@cogs){
    my $n = 0;
    foreach my $name(@names){
        if($hashCogName{$name}{$cog} ne undef){
            
            $combinedSeqs[$n] .= $hashCogName{$name}{$cog};
            $mapped[$n][$c] = 1;
            
        }
        else{
            print STDERR "$name $cog $n $c\n";
            $combinedSeqs[$n] .= "-" x $lengths[$c];
            $mapped[$n][$c] = 0;
        }
        $n++;
    }
    $c++;
}

my $n = 0;
foreach my $name(@names){
    print ">$name\n$combinedSeqs[$n]\n";
    my $cString = join(",",@{$mapped[$n]});
    print STDERR "$cString\n";
    $n++;
}
