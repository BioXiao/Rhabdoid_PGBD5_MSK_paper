#! /usr/bin/perl

use strict;

my $total_len=0;
my %gene_len=();
my %gene_coor=();
open (input, "</home/zhuangj/reference_genomes/genes.bed") or die "Can't open genes.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    $gene_len{$a[3]}=$a[2]-$a[1];
    $gene_coor{$a[3]}="$a[0]\:$a[1]\:$a[2]";
    $total_len += $gene_len{$a[3]};
}
close input;

my %SV=();
my @files=<*.SVs.gene.overlap>;
foreach my $file (@files) {
    my $title=$file;
    $title =~ s/.SVs.gene.overlap//;
	
    my %num_affect_gene=();
    open (input, "<$file") or die "Can't open $file since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	if (! defined $num_affect_gene{$a[4]}) {$num_affect_gene{$a[4]}=$a[11];}
	elsif ($num_affect_gene{$a[4]} !~ /$a[11]/) {$num_affect_gene{$a[4]} .= ",$a[11]";}
    }
    close input;
    
    open (input, "<$file") or die "Can't open $file since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
#	my $id=$title;
	my $id="$title\:$a[4]";
	my @b=split(/\,/, $num_affect_gene{$a[4]});
	
#	if ($a[6] < 0.2) {next;}
	    
	if ($SV{$a[11]} !~ /$id,/) {
	    $SV{$a[11]} .= "$id,";
	}
    }
    close input;
}
    
my %num=();
while ((my $key, my $value) = each (%SV)) {
    chop($value);
    my @x=split(/\,/, $value);
    my $all=$#x+1;
    $num{$key}=$all;
#    if (defined $num{$key}) {$num{$key} .= ",$sample\:$all";}
#    else {$num{$key}="$sample\:$all";}
}

open (output, ">mutation.recurrent.table") or die "Can't open mutation.recurrent.table since $!\n";
print output "GENE\tSV number\tGene length\tGene coordinate\n";

while ((my $key, my $value) = each (%num)) {
    print output "$key\t$value\t$gene_len{$key}\t$gene_coor{$key}\n";
}
close output;
