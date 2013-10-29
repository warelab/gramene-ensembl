#!/usr/local/bin/perl 

=head1 NAME

get_gi.pl - get list of accessions and gi number from MSU site

=cut



use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use HTML::Parser;
use LWP;

=head1 SYNOPSIS

get_gi.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --chroms		number of chromosomes [12]
    --url		url to fetch to get links (chromosome # will be appended)


=head1 OPTIONS

=over 4


=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS


=cut
my $nchroms=12;
my $script='http://rice.plantbiology.msu.edu/pseudomolecules/ordered_bac_';
    {  #Argument Processing
	my $help=0;
	my $man=0;
	GetOptions( "help|?"=>\$help,"man"=>\$man
		   ,"chroms=i" => \$nchroms,
		   ,"url=s" => \$script,
		  )
	  or pod2usage(2);
	pod2usage(-verbose => 2) if $man;
        pod2usage(1) if $help;
    }

my $user_agent=LWP::UserAgent->new;
my $parser=HTML::Parser->new(api_version=>3);
$parser->handler( start=>\&handle_link, 'self,tagname,attr');

for my $chrom (1 .. $nchroms) {
    my $request=HTTP::Request->new(GET=>$script.$chrom.".shtml");
    my $response =$user_agent->request($request);
    $response->is_success or print STDERR "Chromosome $chrom : ",$response->messge,"\n" and next ;

#print "===== $chrom =====\n",$response->content,"\n===========\n=========\n";
    $parser->parse($response->content);
    $parser->eof;
}

exit;
    
# end of main program

########################## subroutines ######################################

sub handle_link {

    my($parser,$tagname,$attr)=@_;

#    print "$tagname:",($attr->{href} || ''),"\n";


    return unless $tagname eq 'a'
              and $attr->{href};

   if ( 
       $attr->{href} =~ m!/cgi-bin/pseudoBAC_view.pl\?BAC=\S+! 
      ) {

     #/cgi-bin/pseudoBAC_view.pl?BAC=OSJNBb0077A02

       $parser->handler(text=> \&handle_text_clone_name, 'self,dtext' );

   } elsif ( $attr->{href} =~ /www.ncbi.nlm.nih.gov.*val=(\d*)/ ) {
     #http://www.ncbi.nlm.nih.gov:80/entrez/viewer.cgi?val=31414484&view=gb

       $parser->{ginum}=$1;
       $parser->handler(text=> \&handle_text, 'self,dtext' );

   } 
}


sub handle_text {
    my ($parser,$text)=@_;
    print "$text\t",$parser->{ginum},"\n";
    undef $parser->{ginum};
    $parser->handler(text => undef);
}

sub handle_text_clone_name {
    my ($parser,$text)=@_;
    print "$text\t";
    $parser->handler(text => undef);
}

__END__

=head1 OUTPUT

=item B<Standard Output>

=head1 NOTES
    

=head1 AUTHOR

   Sharon Wei

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

