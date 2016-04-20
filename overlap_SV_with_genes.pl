#! /usr/bin/perl

use strict;

my @files=<*.SVs.summary>;
foreach my $file (@files) {
    my $sample=$file;
    $sample =~ s/.SVs.summary//;

    my %somatic=();
    open (input, "<$sample.somatic") or die "Can't open $sample.somatic since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	$somatic{$a[0]}=1;
    }
    close input;

    my %processed=();
    open (input, "<$file") or die "Can't open $file since $!\n";
    open (output, ">$sample.SVs.summary.somatic") or die "Can't open $sample.SVs.summary.somatic since $!\n";
    open (output2, ">$sample.SVs.summary.somatic.inversion") or die "Can't open $sample.SVs.summary.somatic.inversion since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	my @b=split(/\:/, $a[3]);
	my $contig="$b[1]\:$b[2]\:$b[3]\:$b[4]";
	if (($somatic{$contig} == 1) && (! defined $processed{$a[4]})) {
	    if (($a[4] =~ /INVERSION/) || ($a[4] =~ /TRANSLOCATION/)) {
		chop($b[2]); chop($b[4]);
		my $lower1=$b[2]-1;
		my $lower2=$b[4]-1;
		print output "$b[1]\t$lower1\t$b[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\n";
		print output "$b[1]\t$lower2\t$b[4]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\n";

		if ($a[4] =~ /INVERSION/) {
		    print output2 "$line\n";
		}
	    }
	    else {
		print output "$line\n";
	    }
	    $processed{$a[4]}=1;
	}
    }
    close input;
    close output;
    close output2;
	
    system("intersectBed -a $sample.SVs.summary.somatic -b ~/reference_genomes/promoters.bed -wo | uniq > $sample.SVs.prom.overlap");
    system("intersectBed -a $sample.SVs.summary.somatic -b ~/reference_genomes/genes.bed -wo | uniq > $sample.SVs.gene.overlap");
    system("intersectBed -a $sample.SVs.summary.somatic -b ~/reference_genomes/exons.bed -wo | uniq > $sample.SVs.exon.overlap");
    system("intersectBed -a $sample.SVs.summary.somatic.inversion -b ~/reference_genomes/exons.bed -wo | uniq > tmp");

    my %eligible_inversion=();
    open (input, "<$sample.SVs.gene.overlap") or die "Can't open $sample.SVs.gene.overlap since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	if ($a[4] =~ /INVERSION/) {
	    $eligible_inversion{$a[4]} .= "$a[11],";
	}
    }
    close input;

    open (input, "<tmp") or die "Can't open tmp since $!\n";
    open (output, ">>$sample.SVs.exon.overlap") or die "Can't open $sample.SVs.exon.overlap since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	if ($eligible_inversion{$a[4]} =~ /$a[11]/) {
	    print output "$line\n";
	}
    }
    close input;
    close output;
	    
    system("rm tmp $sample.SVs.summary.somatic*");
}
