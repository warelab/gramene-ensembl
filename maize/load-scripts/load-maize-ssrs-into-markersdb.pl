#!/usr/local/bin/perl

#marker_id       marker_name     map_name        marker_type     analysis        strand  start   end
#5892025 RM10001 chr_1   ssr     e-PCR   0       24509   24894
#5892026 RM10002 chr_1   ssr     e-PCR   0       58748   58919

# map_set_id 136 | GR TIGR Assm IRGSP Seq 2005 |


# vim: tw=78: sw=4: ts=4: et: 

# $Id: load-maize-ssrs-into-markersdb.pl,v 1.10 2007-09-11 19:48:56 pasternak Exp $

use strict;
use warnings;
use English qw( -no_match_vars );
use File::Basename;
use Getopt::Long;
use Gramene::CDBI::Markers;

#use lib "/home/weix/scripts/markers/lib/gramene/lib/perl";
use Gramene::Marker::DB;
use Pod::Usage;
use Readonly;
use Text::RecordParser::Tab;

Readonly my $VERSION => sprintf '%d.%02d', 
                        qq$Revision: 1.10 $ =~ /(\d+)\.(\d+)/;
Readonly my $MARKER_TYPE => 'ssr';
#Readonly my $MAP_SET_ID => 136;

my ( $help, $man_page, $show_version, $conffile, $analysisID, $MAP_SET_ID);
GetOptions(
    'help'    => \$help,
    'man'     => \$man_page,
    'version' => \$show_version,
   'config_file:s'  => \$conffile,
    'analysis_id:i'   => \$analysisID,
    'map_set_id:1'    => \$MAP_SET_ID,	   
) or pod2usage(2);

if ( $help || $man_page ) {
    pod2usage({
        -exitval => 0,
        -verbose => $man_page ? 2 : 1
    });
}; 

if ( $show_version ) {
    my $prog = basename( $PROGRAM_NAME );
    print "$prog v$VERSION\n";
    exit 0;
}

if( defined $conffile ){ $ENV{GrameneConfPath} = $conffile; }

my @files = @ARGV or die "No input files\n";
my $p     = Text::RecordParser::Tab->new;
my $mdb   = Gramene::Marker::DB->new;

my ( $num_files, $num_markers ) = ( 0, 0 );
for my $file ( @files ) {
    $p->filename( $file );
    
    my $line_num = 0;

    while ( my $rec = $p->fetchrow_hashref ) {
        $line_num++;

        my $marker_id   = $rec->{'marker_id'}      ||  0;
	if ( !$marker_id ) {
            complain("No marker id, $file line $line_num");
	    next;
        }

	my $analysis    = $rec->{'analysis'}       || '';

        my $start       = $rec->{'start'}   || '';
        my $end         = $rec->{'end'}     || '';
        my $strand      = $start > $end ? -1 : 1;

	my $map_name    = $rec->{map_name} || '';
	my $tmp         = $map_name;
	if($tmp =~ /\d+/){
	  $map_name = "Chr. $&";
	}
	print "map_set_id => $MAP_SET_ID, map_name   => $map_name\n";
#	next;
	my $map_id      = $mdb->find_or_create_Map(
					     map_set_id => $MAP_SET_ID,
					     map_name   => $map_name,
					    );

	if( !$map_id ) {
	  complain( "Cannot found map id for map_name $map_name, $file line $line_num");
	  next;
	};
	
	my $marker_obj;
	eval{
	  $marker_obj = $mdb->retrieve_Marker($marker_id);
	};

	if($@ || !$marker_obj){
	  complain("Cannot found marker for $marker_id, $file line $line_num\n$@");
	  next;
	}
	my $analysis_id = $analysisID ? $analysisID : $marker_obj->analysis_id();
	my $display_synonym_id = $marker_obj->display_synonym_id();
	my $mapping_id = $mdb->find_or_create_Mapping( {
					      marker_id  => $marker_id,
					      map_id     => $map_id,
					      start      => $start,
					      end        => $end,
					      strand     => $strand,
					      analysis_id => $analysis_id,
					      display_synonym_id => $display_synonym_id,
					      } );
	print "Added $mapping_id\n";


    }

    $num_files++;
    $num_markers += $line_num;
}

print "Done, processed $num_files files, $num_markers markers.\n";

sub complain {
    print STDERR @_, "\n";
}

__END__

# ----------------------------------------------------
=head1 NAME

load-mappings.pl - a script

=head1 VERSION

This documentation refers to load-mappings.pl version $Revision: 1.10 $

=head1 SYNOPSIS

  load-mappings.pl 

Options:

  --help        Show brief help and exit
  --man         Show full documentation
  --version     Show version and exit
  --config_file the config file
  --analysis_id user specified analysis id, othere wise the mapping record will use the analysis id associated with the marker
  --map_set_id  the map_set_id in the markers db

=head1 DESCRIPTION

Describe what the script does, what input it expects, what output it
creates, etc.

=head1 SEE ALSO

perl.

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut
