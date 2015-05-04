#!/usr/local/bin/perl -w 

=pod

=head1 NAME

rename_genes_inGFF3.pl - rename the genes and transcripts names in the gff3 files
                         according to a name mapping table

=head1 SYNOPSIS

  load_genes_from_jgi_gff3.pl rename_file gff_file

Maintained by Sharon Wei <weix@cshl.edu>

=cut


use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Readonly;
use List::Util qw( first );
use List::MoreUtils;

Readonly my @NAME_FILEDS => qw(NAME ALIAS ID);
Readonly my $UTR_REGEX   => qr{UTR}xms;

my $help=0;
my $man=0;

GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
       )
  or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

my $rename_file = $ARGV[0] || pod2usage("\nNeed the path to a name mapping file\n");
my $gff_file = $ARGV[1] || pod2usage("\nNeed the path to a gff file\n");

map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
} $gff_file, $rename_file;

open my $rfh, $rename_file or die "Cannot open $rename_file";
my %rename_table = map{ chomp; uc $_; split ' ';} (<$rfh>);

#map {print "$_ => $rename_table{$_}"} keys %rename_table ;

# Create a GFF stream
my  $GFF_HANDLE = IO::File->new("< $gff_file")
      or die( "Could not read $gff_file: $!" );

my $rename_regex = join '|', reverse sort keys %rename_table;
print "# This GFF3 file is extracted from $gff_file and rename the gene names according to mappings in $rename_file\n#\n";

while( my $line = $GFF_HANDLE->getline ){
  # Skip comment and empty lines
  next if ( $line =~ /\#/ ); 
  next if ( $line =~ /^\s+/ );

  next unless ( $line =~ /($rename_regex)\b/i);

  my $matched_name = $&;
  my $new_name = $rename_table{ uc $matched_name };

  $line =~ s/$matched_name\b/$new_name/ig;

  print $line;

}

