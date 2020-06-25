#!/usr/bin/perl


$minCount = 10;
$minFrac = 0.5;
while($line = <STDIN>){
    chomp($line);
    #k141_373698,23041,1,480,23041,1.0
    my @tokens = split(/,/,$line);

    my $id = $tokens[0];

    my $count = $tokens[1];

    my $N = $tokens[2];

    if ($N >0){
#k141_508351,763,1,1308,763,1.0
        if ($count > $minCount){
            my $f = $tokens[4];
            if ($f > $minFrac){
                print "$id,$tokens[3]\n";
            }
        }
    }
}
