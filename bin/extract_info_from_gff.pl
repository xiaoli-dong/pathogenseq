#!/usr/bin/env perl
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
my ($gff, $sample_name);

&GetOptions(
    "g=s" =>\$gff,
    "n=s" =>\$sample_name
    );

($gff) || 
    die "Name:\n".
    "by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  Get feature count from gff file\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -g <gff file>\n".
    "  -n <sample_name>\n";

    
my %stat = ();
open(IN, '<', $gff) or die "Could not opne $gff to read, $!\n";

while (<IN>) {
    
    last if m/^##FASTA/;
    next if /^#/;
    chomp;
    my @l = split(/\t/, $_);
    $stat{$l[2]}++;
}
print "Sample_id\t";
print join("\t", sort keys %stat);
print "\n";
my @values = ();
foreach my $k (sort keys %stat) {
    push(@values, $stat{$k});
}
print "$sample_name\t";
print join("\t", @values);
print "\n";



