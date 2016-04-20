#! /urs/bin/perl

use strict;

my %lines=();
my %freqs=();
my @samples=();
my @files=<PA*/*.SVs.vcf>;
foreach my $file (@files) {
    my @z=split(/\//, $file);
    push (@samples, $z[0]);
    my %tissue=();
    open (input, "<$z[0]/$z[0].bp.mechanisms") or die "Can't open $z[0].bp.mechanisms since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	$tissue{$a[0]}=$a[3];
    }
    close input;

    open (input, "<$file") or die "Can't open $file since $!\n";
    while (my $line=<input>) {
	if ($line =~ /^#/) {next;}
	chomp($line);
	my @a=split(/\t/, $line);
	my @b=split(/\:/, $a[2]);
	my $strand1=substr($b[2], length($b[2])-1, 1);
	my $strand2=substr($b[4], length($b[4])-1, 1);
	my $oppo1=$b[2];
	my $oppo2=$b[4];
	if ($strand1 eq "+") {$oppo1 =~ s/\+/\-/;}
	else {$oppo1 =~ s/\-/\+/;}
	if ($strand2 eq "+") {$oppo2 =~ s/\+/\-/;}
	else {$oppo2 =~ s/\-/\+/;}
	my $id1="$b[1]\:$b[2]\:$b[3]\:$b[4]";
	my $id2="$b[3]\:$b[4]\:$b[1]\:$b[2]";
	my $id3="$b[1]\:$oppo1\:$b[3]\:$oppo2";
	my $id4="$b[3]\:$oppo2\:$b[1]\:$oppo1";

	if (! defined $tissue{$a[2]}) {$tissue{$a[2]}="filtered";}

	if (defined $lines{$id1}) {
	    $freqs{"$id1\:$z[0]"}=$a[$#a].":$tissue{$a[2]}";
	}
	elsif (defined $lines{$id2}) {
	    $freqs{"$id2\:$z[0]"}=$a[$#a].":$tissue{$a[2]}";
	}
	elsif (defined $lines{$id3}) {
	    $freqs{"$id3\:$z[0]"}=$a[$#a].":$tissue{$a[2]}";
	}
	elsif (defined $lines{$id4}) {
	    $freqs{"$id4\:$z[0]"}=$a[$#a].":$tissue{$a[2]}";
	}
	else {
	    $freqs{"$id1\:$z[0]"}=$a[$#a].":$tissue{$a[2]}";
	    pop @a;
	    $lines{$id1}=join("\t",@a);
	    $lines{$id1} =~ s/$b[0]/$z[0]\:$b[0]/g;
	    $lines{$id1} =~ s/REFSUP/REFSUP:STATUS/;
	}
    }
    close input;
}

open (output, ">aggregated.SVs.vcf") or die "Can't open aggregated.SVs.vcf since $!\n";
print output "##fileformat=VCFv4.1\n##reference=hg19\n##assembly=$ARGV[0].contigs.fa\n";
print output "##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n";
print output "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n";
print output "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n";
print output "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
print output "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n";
print output "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n";
print output "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
print output "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print output "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n";
print output "##INFO=<ID=MATEORI,Number=1,Type=String,Description=\"Orientation of mate breakends\">\n";
print output "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">\n";
print output "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print output "##ALT=<ID=INS:ME,Description=\"Insertion of TE element\">\n";
print output "##FILTER=<ID=f0,Description=\"Frequency equals 0\">\n";
print output "##FILTER=<ID=r10,Description=\"Total supporting reads fewer than 10\">\n";
print output "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print output "##FORMAT=<ID=GF,Number=1,Type=Float,Description=\"Genotype frequency\">\n";
print output "##FORMAT=<ID=BKSUP,Number=1,Type=Integer,Description=\"Number of reads supporting the breakends\">\n";
print output "##FORMAT=<ID=REFSUP,Number=1,Type=Integer,Description=\"Number of reads supporting reference\">\n";
print output "##FORMAT=<ID=STATUS,Number=1,Type=String,Description=\"SV status\">\n";
print output "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach my $sample (@samples) {
    print output "\t$sample";
}
print output "\n";
while ((my $key, my $value) = each (%lines)) {
    print output "$value";
    foreach my $sample (@samples) {
	my $key1="$key\:$sample";
	print output "\t$freqs{$key1}";
    }
    print output "\n";
}
close output;

    
