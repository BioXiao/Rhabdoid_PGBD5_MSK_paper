#! /usr/bin/perl

use strict;
use List::Util qw(sum);
use Bio::Seq;

my %somatic=();
my %mechanism=();
open (input, "<RPE-GFP-PGBD5.bp.mechanisms") or die "Can't open RPE-GFP-PGBD5.bp.mechanisms since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
#    if ($a[3] eq "somatic") {
	$somatic{$a[0]}=$a[2]; 
	$mechanism{$a[0]}=$a[1];
#    }
}
close input;

my %ref1=();
my %ref2=();

open (input, "<RPE-GFP-PGBD5.reference1.fa") or die "Can't open RPE-GFP-PGBD5.reference1.fa since $!\n";
open (output, ">left.side.seqs.fa") or die "Can't open left.side.seqs.fa since $!\n";
my $i=1;
while (my $line=<input>) {
    chomp($line);
    my $contig=$line;
    $contig =~ s/:Reference1//;
    $contig =~ s/>//;
    if (defined $somatic{$contig}) {
	print output "$line\:$i\n";
	$line=<input>;
	chomp($line);
	$line=substr($line, 450, 101);
	print output "$line\n";
	$ref1{$contig}=$line;
	$i++;
    }
}
close input;
close output;

open (input, "<RPE-GFP-PGBD5.reference2.fa") or die "Can't open RPE-GFP-PGBD5.reference2.fa since $!\n";
open (output, ">right.side.seqs.fa") or die "Can't open right.side.seqs.fa since $!\n";
my $i=1;
while (my $line=<input>) {
    chomp($line);
    my $contig=$line;
    $contig =~ s/:Reference2//;
    $contig =~ s/>//;
    if (defined $somatic{$contig}) {
        print output "$line\:$i\n";
        $line=<input>;
        chomp($line);
	$line=substr($line, 450, 101);
        print output "$line\n";
	$ref2{$contig}=$line;
        $i++;
    }
}
close input;
close output;

