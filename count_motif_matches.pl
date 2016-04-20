#! /usr/bin/perl

use strict;

my $somatic_total=0;
my $germline_total=0;
my %somatic=();
open (input, "<RPE-GFP-PGBD5.bp.mechanisms") or die "Can't open RPE-GFP-PGBD5.bp.mechanisms since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    if ($a[3] eq "somatic") {
	$somatic{$a[0]}=1;
	$somatic_total++;
    }
    else {
	$germline_total++;
    }
}
close input;

my %soma=();
my %germ=();
my $somatic_hits=0;
my $germline_hits=0;
open (input, "<fimo.txt") or die "Can't open fimo.txt since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^#/) {next;}

    my @a=split(/\t/, $line);
    if ($a[7] >= 0.1) {next;}
    if (abs($a[2]-51) >= 15 && abs($a[3]-51) >= 15) {next;}

    my @b=split(/\:/, $a[1]);
    my $key="$b[0]\:$b[1]\:$b[2]\:$b[3]\:$b[4]";
    if (defined $somatic{$key}) {
	if (! defined $soma{$key}) {
	    $somatic_hits++;
	    $soma{$key}=1;
	}
    }
    else {
	if (! defined $germ{$key}) {
	    $germline_hits++;
	    $germ{$key}=1;
	}
    }
}
close input;

print "Somatic hits: $somatic_hits\t$somatic_total\nGermline hits: $germline_hits\t$germline_total\n";
