#!/usr/bin/env perl
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

my ($fasta, $assembly_info);

&GetOptions(
    "f=s" =>\$fasta,
    "i=s" =>\$assembly_info
    );

($fasta && $assembly_info) || 
    die "Name:\n".
    "by Xiaoli Dong <xiaoli.dong\@albertaprecisionlabs.ca>\n".
    "Synopsis:\n".
    "  reset the start position for circular genome to the middle of the genome\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -f <fasta assembly >\n".
    "  -i <assembly info file from flye>\n";



open(INFO, $assembly_info) or die "Could not open $assembly_info file to read, $!\n";

=This is example assembly_info.txt file
#seq_name       length  cov.    circ.   repeat  mult.   alt_group       graph_path
contig_1        5567255 85      Y       N       1       *       1
contig_2        64213   500     Y       Y       5       *       2
contig_4        12254   672     Y       Y       8       *       4
contig_3        4160    926     Y       Y       11      *       3
=cut

my %circgenomes = ();

while (<INFO>) {
    next if !/\S/;
    next if /^#/;
    chomp;
    my @l = split(/\t/, $_);

    # reset circular genome
    if($l[3] eq "Y"){
        $circgenomes{$l[0]} = $l[1];
    }
}
close(INFO);

if ($fasta =~ /.gz$/) {
    open(FASTA, "gunzip -c $fasta |") or die "Could not open $fasta file to read, $!\n";
}
else{
    open(FASTA, $fasta) or die "Could not open $fasta file to read, $!\n";
}

$/ = "\n>";
while (<FASTA>) {
    chomp;
    if ( my ( $seqname, $other, $seq ) = /^>?(\S+)(.*?)\n(.*)/s ) {
        
        $seq =~ s/\s+//g;
        if ( exists $circgenomes{$seqname} ) {
            my $c_length =$circgenomes{$seqname};
            my $start = int ($c_length/2);
            my $first_part = substr($seq, 0, $start);
            my $second_part = substr($seq, $start);
            #print STDERR $seqname, "reset...\n";
            print ">$seqname$other\n$second_part$first_part\n";
        }
        else{
            print ">$seqname$other\n$seq\n";
        }
    }

}

$/ = "\n";
close(FASTA);