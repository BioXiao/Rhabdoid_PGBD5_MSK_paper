#! /usr/bin/perl

use strict;

open (input, "<somatic.SVs.counts") or die "Can't open somatic.SVs.counts since $!\n";
open (output, ">recurrent.regions.bed") or die "Can't open recurrent.regions.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    if ($a[3] >= 2) {
	print output "$line\n";
    }
}
close input;
close output;

my @files=<somatic_SVs_summary/*.SVs.summary>;
my %bk1=();
my %bk2=();
foreach my $file (@files) {
    open (input, "<$file") or die "Can't open $file since $!\n";
    open (output1, ">bk1.bed") or die "Can't open bk1.bed since $!\n";
    open (output2, ">bk2.bed") or die "Can't open bk2.bed since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	my @b=split(/\:/, $a[3]);
	my $lower=$b[2]-1;
	my $upper=$b[2]+1;
	print output1 "$b[1]\t$lower\t$upper\t$a[3]\n";

	$lower=$b[4]-1;
	$upper=$b[4]+1;
	print output2 "$b[3]\t$lower\t$upper\t$a[3]\n";
    }
    close input;
    close output1;
    close output2;

    system("bedtools intersect -a bk1.bed -b recurrent.regions.bed -wo > tmp.bed");
    open (input, "<tmp.bed") or die "Can't open tmp.bed since $!\n";
    while (my $line=<input>) {
        chomp($line);
        my @a=split(/\t/, $line);
	if (! defined $bk1{$a[3]}) {
	    $bk1{$a[3]}=$a[7];
	}
    }
    close input;

    system("bedtools intersect -a bk2.bed -b recurrent.regions.bed -wo > tmp.bed");
    open (input, "<tmp.bed") or die "Can't open tmp.bed since $!\n";
    while (my $line=<input>) {
        chomp($line);
        my @a=split(/\t/, $line);
	if (! defined $bk2{$a[3]}) {
	    $bk2{$a[3]}=$a[7];
	}
    }
    close input;

    system("rm tmp.bed bk1.bed bk2.bed");
}

open (output, ">recurrent.region.somatic.SVs.txt") or die "Can't open recurrent.region.somatic.SVs.txt since $!\n";
while ((my $key, my $value) = each (%bk1)) {
    if ($value >= 3 && $bk2{$key} >= 3) {
	print output "$key\n";
    }
}
close output;

system("rm recurrent.regions.bed");
