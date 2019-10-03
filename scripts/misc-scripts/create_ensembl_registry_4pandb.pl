#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use Pod::Usage;
use Readonly;

my ( $help, $man_page, $dbname_file, $user, $pass, $host, $port, $tax_species );
GetOptions(
    'help' => \$help,
    'man'  => \$man_page,
    'user=s' => \$user,
    'pass=s' => \$pass,
    'host=s' => \$host,
    'port=i' => \$port,
    'tax_species=s' => \$tax_species,
    '-dbname_file=s' => \$dbname_file,
) or pod2usage(2);

if ( $help || $man_page ) {
    pod2usage({
        -exitval => 0,
        -verbose => $man_page ? 2 : 1
    });
}; 


$user ||= 'gramene_web';
$pass ||= 'gram3n3';
$host ||=  'cabot';
$port ||= '3306';

Readonly my $header => q{
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;

};

my $db_params =  << "END_DBPAR";
my \$def_user = '$user';
my \$def_pass = '$pass';
my \$def_host = '$host';
my \$def_port = '$port';

END_DBPAR

print "tax_species=$tax_species\n";

my $dbfile_fh;
my %dbadaptor_hash;

open $dbfile_fh, $dbname_file or die "cannot open $dbname_file";

while (my $dbline = <$dbfile_fh>){

    chomp ($dbline);

    my @f = split ' ', $dbline;
    my $dbname = lc $f[0];
    $dbname =~ s/\s//g;

    if( $dbname =~ /([a-z0-9]+_[a-z0-9]+)_(core|variation|funcgen)_/ ){
	$dbadaptor_hash{$1}{group} = $2;
	$dbadaptor_hash{$1}{dbname} = $dbname;

    }else{
	warn ("not recognized dbname $dbname\n")
    }
    
}


my $dbconnection_string;
for my $species (keys %dbadaptor_hash){

    $dbconnection_string = qq{
Bio::EnsEMBL::DBSQL::DBAdaptor->new(
   -species => $tax_species ? $tax_species : $species,
   -group   => $dbadaptor_hash{$species}{group},
   -port    => \$def_port,
   -host    => \$def_host,
   -user    => \$def_user,
   -pass    => \$def_pass,
   -dbname  => $dbadaptor_hash{$species}{dbname},
	);
    };

   open my $fh, ">$species.reg" or die "Cannot open $species.reg for writing";
   print $fh "$header\n$db_params$dbconnection_string\n1;\n";
   close $fh;

};


__END__

# ----------------------------------------------------

=pod

=head1 NAME

create_ensembl_registry_4db.pl - a script

=head1 SYNOPSIS

  create_ensembl_registry_4db.pl 

Options:

  --help   Show brief help and exit
  --man    Show full documentation
  --dbname_file the file containing the database names you want to create dbAdaptor for.

=head1 DESCRIPTION

Describe what the script does, what input it expects, what output it
creates, etc.

=head1 SEE ALSO

perl.

=head1 AUTHOR

weix E<lt>weix@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2013 Cold Spring Harbor Laboratory

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
