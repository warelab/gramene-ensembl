#! /usr/bin/env perl


my %count;

while(<>){
    chomp;
    $count{$_}++;
}

for (keys %count){
    print join("\t", $count{$_}, $_,), "\n";
}
