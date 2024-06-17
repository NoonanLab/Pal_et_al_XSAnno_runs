#!/usr/bin/perl
use strict;
use warnings;

# Check for input BED file argument
if (@ARGV != 1) {
    die "Usage: $0 input.bed\n";
}

# Input BED file
my $inbed = $ARGV[0];

# Open and read the input BED file
open(IN, $inbed) or die "Could not open input BED file: $!";
while (<IN>) {
    chomp;
    my @t = split /\t/;

    # Extract information from the 4th column
    my @b = split /\|/, $t[3];
    my ($gene_id, $transcript_id, $gene_name, $transcript_name, $exon_id) = @b;

    # Adjust start coordinate to be 1-based
    my $st = $t[1] + 1;

    # Print in the desired GTF-like format
    print "$t[0]\tXSAnno\texon\t$st\t$t[2]\t$t[4]\t$t[5]\t.\t";
    print "gene_name \"$gene_name\" gene_id \"$gene_id\" transcript_name \"$transcript_name\" transcript_id \"$transcript_id\" exon_id \"$exon_id\"\n";
}
close(IN);
