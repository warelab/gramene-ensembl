# the following is actually due to transcript-less genes - i.e. just features!?!
# PROBLEM SAVING GENE PBANKA_010020: Can't use an undefined value as an ARRAY reference at
# /home/dsth/ensembl_source/ensembl_56/modules//Bio/EnsEMBL/DBSQL/GeneAdaptor.pm line 1176,

package VectorBase::Gff3; # {{{1
use strict;
use warnings;
use Carp;
use List::MoreUtils qw(any);

# fixing CDS
#r perl -pe 'my $line = $_; if (s/^(\S+\t\S+\t)exon(.+)/$1CDS$2/) {print $line}' plasmodium_chabaudi_1.gff3 
# getting grep codes
#r perl -ne 'chomp;print $1, qq{\t\t\t$_\n} if /(#\S+\.\S+#)/; ' Gff3_Parser_OO.pl

=head1 VectorBase::Gff3 

VectorBase::Gff3- Gff3 Stuff

=cut

=head1 VERSION

This document describes VectorBase::Gff3 version 0.0.1

=cut

use version; our $VERSION = qv('0.0.1');

=head1 AUTHOR

Daniel S. T. Hughes 
dsth@cantab.net
dsth@cpan.net
dsth@ebi.ac.uk

=cut

my $gene_save = 0;
my $gene_nosave = 0;
my $gene_skipped = 0;
my $trans_added = 0;
my $trans_notadded = 0;
my $trans_skipped = 0;
my $processed_genes = 0;
my $processed_trans = 0;
my $processed_exons = 0;
my $processed_cds = 0;
my $ignored_utrs = 0;
my $ignored_polys = 0;
my $ignored_others = 0;

#y///////////////////////// constructor ///////////////////////////////////////////////////////////

sub new { 
my ($class, $db, $host, $port, $user, $pw, $args) = @_;

    my $db_ad = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -dbname     => $db, 
        -host       => $host, 
        -port       => $port, 
        -user       => $user,
        -pass   => $pw,
    ) or die qq{\nCould not create database handle.};

print "dbconn parameter is        
        -dbname     => $db, 
        -host       => $host, 
        -port       => $port, 
        -user       => $user,
        -password   => $pw,\n";
    #y create DBConnect object - gives DBI handle functionality for direct sql usage
    my $db_handle = $db_ad->dbc;

    my $self = [undef, undef, [ $db_ad, $db_handle ], [undef, [$args]]];
    #$self = []
    bless $self, $class;

    return $self;
}

#y///////////////////////// parsing file into perl ////////////////////////////////////////////////

sub _parseGff3File {

    my ($self) = @_;

    my $args = $self->[3][1][0];

    my $log = $args->{logfile};
    my $infile = $args->{gff3file};
    my $report_limit = $args->{report_limit};
    my $noforceunique = $args->{noforceunique};

    die qq{\nYou must supply a Gff3 file.} if ( $infile eq q{} );

    my %ID_hashOhash;
    my %diagnostic;


#=fs counts and diagnostics initialisation for dev stages - retire most of this stuff
    $diagnostic{IDs_diagnostic} = [];
    $diagnostic{debuger1} = [];
    $diagnostic{count_uniqIDs} = 0;
    $diagnostic{_child_noID_} = 0;
    $diagnostic{_count_noID_noParent_} = 0;
    $diagnostic{_count_nonuniqueID_} = 0;
    $diagnostic{feature_count} = 0;
    $diagnostic{type_gene} = 0;
    $diagnostic{names_gene} = [];
#=fe

    open my $GFF_hndl, q{<}, qq{$infile} or die qq{\nCannot open Gff3 file for reading.};
    # open my $GFF_hndl, q{<}, qq{./$infile} or die qq{\nCannot open Gff3 file for reading.};
    open my $log_hndl, q{>}, 'gff3.log' or die qq{Cannot open log file gff3.log for writting.};
    open my $error_hndl, q{>}, 'gff3_error.log' or die qq{Cannot open log file gff3_error.log for writting.};

    $self->[3][0][0] = $log_hndl;
    $self->[3][0][1] = $error_hndl;

    my $line = q{====================================================================================================};
    #$self->_log_errors(qq{\n$line\nParsing Gff3 file $infile to native perl structures.\n$line\n}, 1); # turn setting none-error mode with 1
    $self->_log_errors(qq{\n$line\nParsing Gff3 file $infile.\n$line\n}, 1);

    LINE:
    while (my $line = <$GFF_hndl>) {

        chomp $line;

        #y ignore empty lines or those with just spaces/tabs
        if ($line =~ /^[\x20\x09]*$/x) { # not anchored - just if it has one at all
            print $error_hndl qq{* skipping line $. [it\x27s empty]\n};
            #$self->_log_errors(qq{* skipping line $. [it\x27s empty]\n});
            next LINE;
        }

        #y ignore comment lines
        if ($line =~ /^[\x20\x09]*\x23/x) { # not anchored - just if it has one at all
            print $error_hndl qq{* skipping line $. [it\x27s a comment]\n};
            next LINE;
        }

        #y superficial check for Gff3 legal format - check that lines end with ';' allow following spaces
        if ($line !~ /\x3B\x20*$/) {
            die qq{\nMalformed feature line $. (lines must terminate with '\x3B')}
            . qq{\n\nPerhaps try: \"perl -i -pe \'s/([\\w\\d]+)[^\\x3B]\\x20*\$/\$1;\\n/x\' $infile\"};
            next LINE;
        }

        #g parse line into HASH REF - error messages are printed by _lineFeatureSyntaxChecknParse so just next
        my $feature;
        next LINE if !($feature = $self->_lineFeatureSyntaxChecknParse($line));

        #print qq{\nhave a look at the seq id now: }, $feature->{landmark};
        #b put after line ignorers so only counts processed lines - relic?
        $diagnostic{feature_count}++;

        my $uniqID;

        #y unpack only those vars we will actually use
        my $type = $feature->{type};
        my $comments = $feature->{comments};

        # ID search MUST be non-greedy!?!

        #y check for parents using comments column (couldn't be bothered to nest these in attributes
        my @parents = ($feature->{comments} =~ (/Parent=(\S+?)\;/g));
            my $children_store;

        #y check if this feature has an ID
        if (exists $feature->{attributes}{ID}) {
            
            #y unpack ID
            $uniqID = $feature->{attributes}{ID};

#f/ vital that you check on the uniqID as well as children otherwise you
#f/ autovivify the value and the program will stop at the first gene
#y/ let's save the children till later - the single key above should be them
#/ alternatively we can assign individual pieces to the hash instead of feature as a whole that wipes it
            $children_store = $ID_hashOhash{$uniqID}->{children} 
              if ( (exists $ID_hashOhash{$uniqID} ) && (exists $ID_hashOhash{$uniqID}->{children}) );
            #y check if this ID has already been encountered

#f/ fundamental floor!?! you have allowed for different orders - but then
#f/ when you store the data later you wipe OUT the children again
#/ major problem is if the mRNA comes after the things that declare it as a parent. 
#/ i.e. it already existand we carp - but we don't want to wipe them OR change their IDs
#w we know that it exists - but we have to check that its not just cos it was already declared as a parent before itself being declared
#f thus this can go here - or above - i.e. above and it applied to everything so best there
#/ should make it test that the keys is children
            if ( (exists $ID_hashOhash{$uniqID}) && (scalar (keys %{$ID_hashOhash{$uniqID}}) != 1) ) {
            #if  (exists $ID_hashOhash{$uniqID}) {

                #b we absolutely do not allow gene or mRNA types to have non-unique
                #IDs but this of course allows other gen-level or transcript level
                #features to have non-unique IDs - we don't allow that atm - need and
                #option for that with other elsif between elsif and else 

                #y will won't tolerate genes or mRNA with non-unique IDs
                if ( $type eq q{gene} || $type eq q{mRNA} ) { 

                    #b use dies statements - we want to know where the code dies not what routine calls dying code
                    die qq{\nLine $. with ID \'$uniqID\'. This $type cannot have non-unique ID: ${uniqID}:\n$line};
                    #next LINE;
                }
                #y to retain exons/CDSs WITH IDs but NOT unique give them unique ID if they don't have one
                #non-unique we give it a unique one just to retain the data in a hash
                elsif ($type eq q{exon} || $type eq q{CDS}) {

                    #y if its the first time we are seeing this ID start the counter
                    $ID_hashOhash{$uniqID}->{count} = 0 if (!exists $ID_hashOhash{$uniqID}->{count});

                    #/ this is an issue - most of the replicate exons are cos they are in multiple transcripts
                    #/ thus they are identical - perhaps check that all details are exactly the same?!? 
                    #/ except parent?!?

                    #y increment counter (first) and assign 'now' unique ID
                    $uniqID .= q{_nonunique_}.++$ID_hashOhash{$uniqID}->{count}; 

                }
                #b/ turns off forcing of unique ids for stuff - prolly ought to bin if
                #b/ exon or CDS part - i.e. either force it or don't?!?
                elsif ($noforceunique) {
                    $ID_hashOhash{$uniqID}->{count} = 0 if (!exists $ID_hashOhash{$uniqID}->{count});
                    $uniqID .= q{_nonunique_}.++$ID_hashOhash{$uniqID}->{count}; 
                }
                #b we don't want to over-write features - we should probably die
                else { 
                    #/ need to change this to accomodate anything at all?!? - i.e. give unique to anything 
                    #/ that hasn't got a uniq - i.e. not just cds and exons?!? and let it get filtered later?!?
                    #/ HOWEVER mRNAs should always be unique...
                    $self->_log_errors(qq{#ID.NotUnique# $uniqID\n * Skipping feature with non-unique ID: $line}); 
                    $diagnostic{_count_nonuniqueID_}++;
                    next LINE;
                } 

            } # End ID already exists (is unique check)

            #b we know the thing with an id now has a unique one so increment counter
            $diagnostic{count_uniqIDs}++;   # retire
            $diagnostic{type_gene}++ if $type eq q{gene};   # definitely retire
            #b this is a pointless-ish list used to retain all unique ids for genes later cross-ref if there's a prob
            push @{$diagnostic{names_gene}}, $uniqID if $type eq q{gene};   # retire - don't use genes now   

        } 
        #y has no ID, so we check if it has any parents at least - i.e. id-less and parent-less will be ignored
        elsif ($comments =~ (/Parent=(\S+)\;/)) { # we are not actually using the value atm so no need for non-greedy

            $diagnostic{_child_noID_}++;    # retire

            #r create id for idless children
            $uniqID = q{_noID_childOf_}.$1.q{_}.$diagnostic{_child_noID_};
        }
        #r/ has no ID or parent so we ignore it
        else { 
            $self->_log_errors(qq{\n#ID.NOID_NOPARENT#\n* Skipping line $. [We will not allow ID and Parentless features]:\n$line});
            $diagnostic{_count_noID_noParent_}++;
            next LINE;
        } # End of has ID check

        #y definitely have a unique ID even if its just for retention we put into
        # if it doesn't exist we make it exist - condition totally superfluous as is the undefined HASH ref...
        $ID_hashOhash{$uniqID}= {} if (!exists $ID_hashOhash{$uniqID}); 
        #$ID_hashOhash{$uniqID}= {};

        #f/ this is the reason WE LOOSE CHILDREN
        #y/ probably ought to do this in a verbose fashion as you loose everything every time
        #r/ MUST feed the $feature FIRST - i.e. before adding new details or you wipe the hash of previous details
        $ID_hashOhash{$uniqID} = $feature;
        #print qq{\nhave a look at the seq id now: }, $feature->{landmark};

        #b diagnostic function - i.e. create array with all names (can use it if there are too many unique IDs)
        push @{$diagnostic{IDs_diagnostic}}, $uniqID; # retire

        #y we it parents so make it a child of them - we can of course just make list of parents - tree structure...
        #b potential danger - i.e. we are assuming that only bottom level features don't have unique ids...
        push @{$ID_hashOhash{$_}->{children}}, [$uniqID,$ID_hashOhash{$uniqID}] for (@parents);
        #push @{$ID_hashOhash{$_}->{children}}, {$uniqID => $ID_hashOhash{$uniqID}} for(@parents);
        
        #y only put parent array in if it exists
        $ID_hashOhash{$uniqID}->{parents} = [@parents] if (scalar @parents > 0);
        #re-store children
        $ID_hashOhash{$uniqID}->{children} = $children_store if ($children_store);

        #y don't think we've modified these but let's overwrite it anyway?!?
        $ID_hashOhash{$uniqID}->{type} = $type;
        $ID_hashOhash{$uniqID}->{comments} = $comments;

    } # end of line processing using diamond operator

    close $GFF_hndl;

    #y do some reporting and preliminary checks
    $self->_parsing_summary(\%ID_hashOhash, \%diagnostic, $log_hndl, $report_limit);
    
    #f/ feed the object - not copying just passing it - perhaps copy
    #f/ but if copy need to deep copy as its nested...
    $self->[4] = \%ID_hashOhash;

    #return ($diagnostic{names_gene}, $log_hndl);
    return;
    # return (\%ID_hashOhash, $diagnostic{names_gene}, $log_hndl);
}

sub _lineFeatureSyntaxChecknParse {

    my ($self, $line) = @_;  

    my $log_hndl = $self->[3][0][0];
    my $error_hndl = $self->[3][0][1];

    my $args = $self->[3][1][0];

    #y for Gff3 syntax requirements
    if (
      my ($landmark, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = 
      ($line =~/^\x20*([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\x20*$/x)) {

        #r series of tedious health checks

        #y check that coordinates are numeric
        #for ($start, $end) { die qq{start, end must be numeric and defined } if ($_ !~ /^\d+$/); }
        for ($start, $end) { 
            if ($_ !~ /^\d+$/) {
                # print error and return 0 - i.e. that skips
                $self->_log_errors(qq{* skipping line $. as \'start\', \'end\' } 
                  . q{must be numeric and defined for line.}); 
                return 0;
            }
        }
        
        #y score must be numeric or not defined
        $self->_log_errors(qq{* score must be numeric or not defined '.'}) if $score !~ /^(\d+|\.)$/;

        #y check strand orientation
        if ($strand !~ /^[+-?\.]$/) {
            $self->_log_errors(qq{strand must be must '+', '-', '?' or undefined '.' (currently :\'$strand\')});
            return 0;
        }

        #y check phases are permissible - stricly '.' with phase is not allowed for type eq cds
        if ($phase !~/^[012\x2e]$/) { # ignore \x2e '.'
            $self->_log_errors(qq{* skipping line $. as phase must be 0, 1, 2 or \x27\x2e\x27 }
              . qq{(currently: \'$phase\')});
            #next if ($phase !~/^[123\.]$/);
            return 0;
        }

        #y removing trailing space/tab in comments to avoid an attribute ' '
        $attributes =~ s/[\x20\x09]*$//;

        #r/ atm its being read in in interpolation context so you loose the escapes anyway
        my @attrs = split (q{;}, $attributes);  
        
        #y split up attributes into key=>value pairs
        my %attr_hash;
        for my $attribute (@attrs) {
           
            #y clean front and back - why are they separate?!?
            $attribute =~ s/^[\x20\x09]*//; 
            $attribute =~ s/[\x20\x09]*$//;

            if ($attribute =~ /\x20*([^\x09]+)\x20*\x3D\x20*([^\x09]+)\x20*/) { 

                #y store these before you loose them with more matching
                # must allocate all matches before cleaning them - i.e. re-allocate match otherwise
                my $att = $1;
                my $val = $2;

                #y clean up front and back of each - not sure why they're seperate?!?
                $att =~ s/^[\x20\x09]*//; 
                $att =~ s/[\x20\x09]*$//;
                $val =~ s/^[\x20\x09]*//; 
                $val =~ s/[\x20\x09]*$//;
                $attr_hash{$att} = $val; 


            }
            else { 
                $self->_log_errors(qq{* Line $. has problem with attribute [$attribute]}
                . qq{\n\tComment columns should have the form of a series of }
                . qq{\x27attribute_name=atrribute_value\x27 pairs separated by \x27\x3b\x27.}
                . qq{\n\tAny of \x27\x3b \x3d \x25 \x26 \x2c\x27 that are not part of the Gff3 syntax must be escaped.}
                . qq{\n\tIgnoring error - hope the attribute wasn\x27t important.});
            }
        }

        #y attributes won't handle multiple parents and can't be arsed to move the parent code here so undef
        delete $attr_hash{Parent} if (exists $attr_hash{Parent}); # no need for conditional...

        my $feature = {
            landmark => $landmark,
            source => $source,
            type => $type,
            start => $start,
            end => $end,
            score => $score,
            strand => $strand,
            phase => $phase,
            comments => $attributes,
        };

        #print qq{\n this is the seqid before recovering }, $feature->{landmark};
        #y put this one on if it exists
        $feature->{attributes} = {%attr_hash} if %attr_hash;

        return $feature;
        } # end of gff3 match 
        #y they gave us a shitty line
        else { 
            $self->_log_errors(qq{\n* skipping line $. - Malformed feature - 9 cols sep by tab line has a problem - }
            . qq{printing to log file and skipping: $line\n}); 
            return 0;
    }
}

sub _separate_tripartite_and_singleFeatures {
    my $self= shift;

    my %isolatedFeatures;

    #y to fetch all elements forming part of gene models can either delete singleFeatures from [4][0]
    #y or can use has children OR has parents - those - thus [0] head-nodes for those
    #y structures. [1] the single features and [4] will now be the individual features of [0]
    my %tripartiteFeaturesIndividually;
    my %tripartite;

    while ( my ($k, $v) = each %{$self->[4]}) {  
        #b we are making a big change - we want all top level nodes that also have
        #b children - not kids -> basic feature and not gene/biotype p'way
        #b thus compound conditional - using same thing as hode_node_list
        #if ($v->{type} eq q{gene}) {
        #    $just_genes{$k} = $v;
        #b abandon assumption about biotype
        #    $v->{EnsEMBL_biotype} = q{protein_coding};
        
        #o thus now it can best post-fix
        $tripartite{$k} = $v if ( (!exists $v->{parents}) && (exists $v->{children}) );
        $isolatedFeatures{$k} = $k if ( (!exists $v->{parents}) && (!exists $v->{children}) );
        $tripartiteFeaturesIndividually{$k} = $v if ( (exists $v->{parents}) || (exists $v->{children}) );
    }

    #b diagnostic 
    my $n = scalar (keys %tripartite);
    $self->_log_errors(qq{* There are $n gene-level features (features with derives-from relationships }
      . q{to other features but are not derivatives of any others).});

    $self->_log_errors(qq{* Separating $n gene-level features.}, 1);

    #f/ feed the object - that is move the all the data tripartites
    #f/ to one place and single features to another
    #w again we don't copy we send the actual data by ref
    $self->[0] = \%tripartite;
    $self->[1] = \%isolatedFeatures;

    #y re-pack it without the isolated features
    $self->[4] = \%tripartiteFeaturesIndividually;

    #b/ temporarily maintain pointers to ALL the features individually
    #/ i.e. [0] is JUST top-level structures and 
    #/ [1] is JUST bottom structure
    #/ all others are accessed downstream of top-levels in [0]
    #delete $self->[4];

    return;
    # return \%tripartite;
}

#y///////////////////////// misc. parsing methods/subs/////////////////////////////////////////////

sub _gff3ID_inconsistency_check {

    my ($self, $diagnostic) = @_;
        
    my %diagnostic = %{$diagnostic};

    my $error_hndl = $self->[3][0][1];

    my $args = $self->[3][1][0];
    my $report_limit = $args->{report_limit};

    my $line = q{====================================================================================================};
    $self->_log_errors(qq{\n$line\nPROBLEM\n$line\n\nThere are more features stored in the Perl structure than there }
      . qq{we features in the Gff3 file.\nThis is likely to be caused by inconsistencies between feature ID/Parent tags.\nPerhaps some have }
      . qq{an identifier e.g. 'vectorbase|' prepended and others not?\n});

    my @diff = @{&_find_array_differences( \@{$diagnostic{IDs_diagnostic}}, \@{$diagnostic{hashOhash_keys}})};

    my $limit = 0;
    for my $example (@diff) {
        my @things =  grep { /$example/ && $_ !~ /_noID_childOf_/ } @{$diagnostic{hashOhash_keys}};
        @things = map  { $_->[0] } sort { $a->[1] <=> $b->[1] } map  { [$_, length($_)] } @things;
        
        #y print up until report_limit for error_log
        print $error_hndl q{There appears to be discrepencies in ID: };
        print $error_hndl qq{\'$_\', } for (@things);

        #y print up to 10 times for STDOUT
        if ($limit < 10) {
            print STDOUT q{There appears to be discrepencies in ID: };
            print STDOUT qq{\'$_\', } for (@things);
        }

        if ($things[1] =~ /(.*)$things[0](.*)/) { 

            my $prefix = $1;
            my $suffix = $2;

            if ($prefix) { 
                #y print to (STDOUT and log file) inconsistent prefix at end of line for each
                print STDOUT qq{[Appears to have inconsistent prefix \'$prefix\'].\n} if ($limit < 10); 
                print $error_hndl qq{[Appears to have inconsistent prefix \'$prefix\'].\n}; 
            }
            elsif ($suffix) {
                #y print to (STDOUT and log file) inconsistent suffix at end of line for each
                print STDOUT qq{[Appears to have inconsistent suffix \'$suffix\'].\n} if ($limit < 10); 
                print $error_hndl qq{[Appears to have inconsistent suffix \'$suffix\'].\n}; 
            }
            else {
                #y print to (STDOUT and log file) inconsistent prefix/suffix at end of line for each
                print STDOUT qq{[Appears to have inconsistent prefix \'$prefix\' and suffix \'$suffix\'].\n} if ($limit < 10); 
                print $error_hndl qq{[Appears to have inconsistent prefix \'$prefix\' and suffix \'$suffix\'].\n}; 
            }
        }
        else { print STDOUT qq{\n} }

        #y print limit to STDOUT
        print STDOUT qq{\n-----\nReached STDOUT error reporting limit\n-----\n} if ($limit == 10);
        $limit++;
        if ($limit == $report_limit) {
            #y print limit to file
            print $error_hndl qq{\n-----\nReached error log error reporting limit\n------\n};
            last;
        }
    }
    close $error_hndl;
    #die qq{Check ".error_log" and try to fix these errors before re-parsing};
    return;
}

#y/////////////////////////////////////////////////////////////////////////////////////////////////

sub _gene_prefix_detection {

    my $self = shift;

    my $args = $self->[3][1][0];

    my %list = %{$self->[0]};
    #y safety check - we get the gene names as an ARRAY (i.e. not just distinct values)
    my @other_gene_names = (keys %list);
    #b duh - _log_errors only takes a single arg
    #$self->_log_errors(q{* Separating }, scalar @other_gene_names, qq{ genes});
    my $n = scalar @other_gene_names;

    my $prefixes = &_check_for_prefixes(\@other_gene_names,5, 11);

    my $line = q{====================================================================================================};
    print STDOUT qq{\n$line\nChecking sub-set of gene-level IDs for prefixes between 5 and 11 chars.\n$line\n};

    if (scalar (keys %{$prefixes}) > 0) {
        print STDOUT qq{\nPossible gene prefixes detected:\n};
    
        #r prefix hash:
        # HASH REFERENCE (0)
        #   |  
        #   |__'8'=>HASH REFERENCE (1) [ '->{8}' ]
        #   |    |  
        #   |    |__'vectorba'=>SCALAR = '30' (2)  [ '->{8}{vectorba}' ]
        #r thus $prefix->{length}{prefixfound} = timefound
        #r or $prefix->{length} = { prefixfound => timesfound }
        #y we know exactly what the length params are - duh don't need key sorting
        # for (sort { $b <=> $a } keys %{$prefixes}) {

        print STDOUT qq{\n* Char length }, 11-$_, qq{; }, (keys %{$prefixes->{11-$_}}) for (0..6);
        print STDOUT qq{\n};

        #y just pass the first one - i.e. pass one scalar in LIST
        my ($p) = (keys %{$prefixes->{11}});

        #/ don't want nasty heredoc in middle so factor into sub - i.e. make it a theredoc
        print &_theredoc_gene_prefix_removal($p);
        die;
        #draw($prefixes->{$_});
        #y for now we just want the keys not the values - i.e. the count of times the prefix was encountered
        # for (0..6) { my @Ps = (keys %{$prefixes->{11-$_}}); print qq{\n* Char length }, 11-$_, qq{; @Ps};; }
        #draw($prefixes->{$_});
        #y for now we just want the keys not the values - i.e. the count of times the prefix was encountered
    }
    else {
        print STDOUT qq{\nNo prefixes auto-detected.\n};
    }
    return;
}

sub _theredoc_gene_prefix_removal {

my $prefix = shift;
my $line = q{----------------------------------------------------------------------------------------------------};

#y want control over tedious layout here so using where doc

my $message = <<"PREFIX";

$line

Remove the prefixes from gene names in Gff3 file and restart.

If the genes have prefixes such as 'vectorbase|E000020' try using something like:

    > perl -pe 's/(ID=)vectorbase\\|(.+?;)/\$1\$2/' filename.gff3

Thus in this case something like (make sure you escape any special chars):

    > perl -pe 's/(ID=)$prefix(.+?;)/\$1\$2/' filename.gff3

Or for the uber-conservative:

    > perl -pe 's/^\\x20*([^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\tID=)vectorbase\\|(.+?;.*)\\x20*\$/\$1\$2/' filename.gff3

$line
     *** If you received this message in error re-run with the '-ignoregeneprefix' option. ***
$line

PREFIX

    return $message;

}

#y/////////////////////////////////////////////////////////////////////////////////////////////////

#b the dullest sub in the world. just don't care...
sub _analyse_comments_column { # {{{1

    my $self = shift;

    my $args = $self->[3][1][0];
    
    #y print tedious message about reclassifying biotypes as Gff3 type
    print &_theredoc_type_reclassification_text;

    chomp (my $exit = <STDIN>);
    #exit if ($exit eq q{n} || $exit eq q{N});
    return if ($exit eq q{n} || $exit eq q{N});

    print qq{* Isolating features that form part of gene-structures and checking attribute types.\n};
 
    #y get comment headers/types
    my @headers = @{$self->_get_geneLevel_comment_types()};

    my $header_number = scalar @headers;

    #/ why is it copying the data?!?
    my %genes = %{$self->[4]};
    # my %genes = %{$self->[0]}; # duh, this is acting only on gene-level features

    COMMENTS:
    while (1) {

        #b this is here just to reprint the list...
        if ($header_number == 0) { 
            print STDOUT qq{\n* The comments section of this Gff3 file contains only ID and Parent headers - there is no need to check for gene types in comments.\n};
            return;
        }
        else {
            print STDOUT qq{\nThe attributes column of features forming gene-structures includes $header_number attribute type(s) [excluding ID and Parent types]:\n\n};
            print STDOUT qq{* \x27$headers[$_]\x27\n} for (0..$#headers);
        }
        
        print STDOUT qq{\nEnter the name of one of these header types to check its values for possible mis-classification of biotype (enter \'c\' to continue)? };

        chomp (my $choice = <STDIN>);

        return if ($choice eq q{c} || $choice eq q{C});
        # last if $choice eq q{q};

        if (!any {$_ eq $choice} @headers) {
            if ($choice eq q{}) { print STDOUT qq{\n* You must enter something!\n}; }
            else { print STDOUT qq{\n\* '$choice\' is not a valid option!\n}; }
            next COMMENTS; # redo;
        }

        #if ($choice eq q{r} || $choice eq q{R}) {
        print STDOUT qq{\nIsolating all distinct values of $choice.\n};

        my %descriptions;
        for my $i (values %genes) { 
            while ( my ($attrs, $val) = each %{$i->{attributes}} ) {
                if ($attrs eq $choice) { 
                    #push @warnings, $val if ($val =~ /ps[eu][ue]do.{0,3}gene/i); # often spelt wrong or has '-'...
                    $descriptions{$val} = 1;
                }
            } 
        } 

        print STDOUT qq{\n* There are }, scalar (keys %descriptions), qq{ distinct values for attribute \'$choice\'.\n};
        my @values = sort (keys %descriptions);

        CHOICE_MAIN:
        while (1) { # not sure why put it in a loop rather than just remaining with goto

            print STDOUT qq{\nDo you wish to (s)earch these values for particular \'biotypes\' or (l)ist }
              . q{all the distinct values? [S/l] };

            chomp (my $what = <STDIN>);
            # local $, = qq{\n};
            # return if ($what eq q{q} || $what eq q{Q});

            if ($what =~ /^[Ll]$/) {
                CHOICE_LIST:
                print STDOUT qq{\nValues of $choice:\n\n};
                print STDOUT qq{[$_] \x27$values[$_]\x27\n} for (0..$#values);

                &_type_check_choice;
            
            }
            else {
                CHOICE_REGEXP:
                print qq{\nEnter a regular expression to search for [default search is for pseudogenes]: };
                chomp (my $regexp = <STDIN>);
                $regexp = qr{ps[eu][ue]do.?gene} if !$regexp;

                my @warnings;

                for my $val (@values) {
                    push @warnings, $val if ($val =~ /ps[eu][ue]do.{0,3}gene/i) ; # often spelt wrong or has '-'...
                    #push @warnings, $val if ($val =~ /ps[eu][ue]do.{0,3}gene/i); # often spelt wrong or has '-'...
                }

                if (scalar @warnings > 0) {
                    my $line = q{----------------------------------------------------------------------------------------------------};
                    print qq{\n$line\n                         *** Regular expression found in }, scalar @warnings, 
                      qq{ entries! ***\n$line\n\nDo you wish to list these entries? [Y/n] };
                    chomp ($what = <STDIN>);
                    if ($what !~ /^[nN]$/) {
                        print STDOUT qq{\nListing matching entries:\n};
                        print STDOUT qq{\n* }, $_ for (@warnings);
                    }

                    #y print one-liner message
                    print &_theredoc_type_reclassification_oneliner($choice);

                    chomp ($what = <STDIN>);
                    if ($what =~ /^[Qq]/) { exit; }
                    next COMMENTS;

                    }
                else {
                    #print qq{\n* Regular expression \'$regexp\' was not found in in any entries!\n\n}
                    print qq{\n* Regular expression was not found in in any entries!\n};

                    &_type_check_choice;

                    }
            }
            # else { next CHOICE_MAIN; }
            # else { goto CHOICE }
        }
    }
    return;
}        

sub _get_geneLevel_comment_types {
    my $self = shift;

    #y scan ALL levels for now for attribute types
    my $data = $self->[4]; # this operates on just the genes - but that's fine
    #my $data = $self->[0]; # this operates on just the genes - but that's fine

    my %descriptions;

    my $args = $self->[3][1][0];

    while ( my ($k, $v) = each %{$data}) {  

            for my $attrs (keys %{$v->{attributes}} ) {
                $descriptions{$attrs} = 1;
        }
    }

    delete $descriptions{ID} if (exists $descriptions{ID});

    my @headers = (keys %descriptions);

    return \@headers;
}

# Also, at this time the only time genes with multiple transcript-level features with differing 
# biotypes are tolerated is when the biotypes are a combination of 'pseudogene' and 
# 'mRNA/protein_coding' - in which case the parent gene is assigned the 'protein_coding' biotype.

sub _theredoc_type_reclassification_text { # {{{1

my $line1 = q{====================================================================================================};
my $line2 = q{----------------------------------------------------------------------------------------------------};

#y must use interpolation-context here! - i.e. "TXT"/qq{TXT}
my $message = <<"TXT";

$line1
Attributes column check.
$line1

Inconsistencies in Gff3 format usage means that the 'biotype' of gene/transcript-level features may 
be encoded as either:

* In the Gff3 'type' column (column 3) e.g.

  'source   source   BIOTYPE   start   end   score   strand   phase   ID=SomeName;Parent=SomeGene;'

* As an attribute value in the Attribute column (column 9) e.g. 
    
  'source   source   mRNA   start   end   score   strand   phase   ID=SomeName;Parent=SomeGene;Note=BIOTYPE;'

Furthermore, 'biotype' may only be encoded at the transcript-level and not gene-level, such that 
gene-level features have type 'gene' in column 3 irrespective of their biotype.

Due to such inconsistencies, this script ignores type for gene-level features and simply identifies 
entities corresponding to EnsEMBL Genes based solely upon a two-tier (gene- & transcript-level 
features) or tripartite structure (gene-, transcript- & exon-level [and cds-level where gene has 
biotype protein_coding] features).

The script then assigns a biotype to newly created EnsEMBL Gene features based upon the bioptype 
given at the transcript level. In the case of genes with transcript-level features with Gff3 'type' 
mRNA they are assigned the EnsEMBL protein_coding biotype - as long as they possess BOTH exons and 
CDS (in which case they VERBOSELY are ignored). 

$line2
CONSEQUENTLY IN ORDER TO AVOID EITHER THE MISCLASSIFICATION OF BIOTYPES OR PSEUDOGENES BEING IGNORED
DUE TO THEIR POSSESSING THE 'Mrna' BIOTYPE WITHOUT HAVING BOTH EXON AND CDS FEATURES IT IS ESSENTIAL 
THAT IN CASES WHERE BIOTYPE IS ASSIGNED AS A FEATURE ATTRIBUTE IT IS CONVERTED TO GFF3 TYPE.
$line2

- Additionally any features at the exon-/cds-level that do not have type 'exon' or 'cds' e.g. 
'five_prime_utr', 'three_prime_utr', 'start_codon'... are ignored.

Do you want to search for mis-classified biotypes in Gff3 features ('n' to exit)? [Y/n] 
TXT

return $message;
}

sub _theredoc_type_reclassification_oneliner {

my $attr = shift;
my $line = q{----------------------------------------------------------------------------------------------------};
my $message = <<"TYPE";


Finding the match may mean that there are genes with 'biotype' imporperly entered as an attribute 
field and not as  the transcript-level feature type.

If this is the case such categorisation will cause these genes to be mis-categorised within the 
EnsEMBL database.

Please fix such problems before commiting this data. Perhaps try a regular expression substitution. 
In the case of pseudogenes categorised by attribute 'Note' you can use something like one of the 
following:

    > perl -pe 's/(.+)(gene|mRNA)(\\t.+Note=pseudogene.+)/\$1pseudogene\$3/' filename.gff3

    > perl -pe 's/(.+)(gene|mRNA)
    > (\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\tNote=pseudogene.+)/
    > \$1pseudogene\$3/x' filename.gff3

    > perl -pe 's/^\\x20*([^\\t]+\\t[^\\t]+\\t)(gene|mRNA)
    > (\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\tID=vectorbase\\|.+?;.*)\\x20*\$/
    > \$1pseudogene\$3/x' 
    > filename.gff3

Just replace 'Note' with this attribute '$attr' and pseudogene with appropriate biotype!

Do you wish check further (a)ttribute mis-classifications or (q)uit now and correct this problem? [A/q]
TYPE

return $message;

}

sub _type_check_choice {

    print STDOUT qq{\nDo you want to (s)earch within these values for a pattern, (l)ist }
      . q{the values again, check a different (a)ttribute or (c)ontinue? [S/l/a/c] };

    chomp (my $what = <STDIN>);

    if ($what =~ /^[Aa]$/) { next COMMENTS; } 
    elsif ($what eq q{l} || $what eq q{L}) { goto CHOICE_LIST; }
    # not in same sub anymore so can't just use return
    elsif ($what eq q{c} || $what eq q{C}) { last COMMENTS; }
    #elsif ($what eq q{c} || $what eq q{C}) { return; }
    #elsif ($what eq q{q} || $what eq q{Q}) { return; }
    else { goto CHOICE_REGEXP }
    # next CHOICE_MAIN;

}
           
#y/////////////////////////////////////////////////////////////////////////////////////////////////

sub _fetch_slices_quick {
    
    #/ bunch of code re-used - factorise out

    my ($self, $cs_name) = @_;

    my %uniq_landmarks;

    # for my $i (values %{$self->[0]}) { $uniq_landmarks{$i->{landmark}} = 1; }
    ################################################################################################
    #/ is this due to passing the tripartitite hierachy and not individually accessible - why now?!?
    # for my $i (values %{$self->[0]}) { 
    for my $i (values %{$self->[4]}) { 
    ################################################################################################
        
#        print qq{\na seqied }, %{$i};
        $uniq_landmarks{$i->{landmark}} = 1; 
    
    }

    my @uniq_landmarks = (keys %uniq_landmarks);
    print qq{\n\nthe unique marks: }, @uniq_landmarks;

    #y generate api slice adaptor
    my $db_slice_adaptor = $self->[2][0]->get_adaptor("Slice") or carp qq{\nCoudln't make a slice adaptor};

    my $line = q{====================================================================================================};
    $self->_log_errors(qq{\n$line\nAttempting to generate EnsEMBL slice objects for all seqids:\n$line\n}, 1);

    my %slices;

    my $CS_flag = 1;

    for my $landmark (@uniq_landmarks) {

    #### why now?!?    
    next if ($landmark eq q{} || $landmark eq q{ });
        #my $slice = $db_slice_adaptor->fetch_by_region($cs, $landmark) or carp qq{\ncould not fetch slice object for $landmark};
        #eval {
        if (my $slice = $db_slice_adaptor->fetch_by_region($cs_name, $landmark)) {
            $self->_log_errors(qq{* Created slice adaptor for $landmark in $cs_name coordinate system.}, 1); 
            $slices{$landmark}  = $slice;
        }
        else {
            #y send this just to STDOUT as it iterative!?!
            #carp qq{could not fetch slice object for $landmark in $cs coordinate system\n}; 
            print STDOUT qq{* Could not fetch slice object for $landmark in $cs_name coordinate system.\n}; 
            #die;
            #r switch flag if anything fails
            $CS_flag = 0;

        }

    }
#    };

    if ($CS_flag == 0) {
        my $line = q{----------------------------------------------------------------------------------------------------};
        print STDOUT qq{\nCould not fetch all slices for these SEQIDs in using Coordinate System $cs_name.\n}
          . qq{\n$line\n                   *** Try re-running without the \'-coordsystem\' option. ***\n$line\n};
        die;
    }

    $self->[2][2][0] = \%slices;

    return;
}

sub _fetch_slices_with_checks {
    my ($self, $db) = @_;

    my $args = $self->[3][1][0];

    #w unpack - don't copy!?!
    #my $db_ad = $vb_Obj->[2][0]; # comment out and qs - neeed 1x
    my $db_handle = $self->[2][1]; 
    #my $genes = $vb_Obj->[0]; # comment out and qs - neeed 1x 
    #my $genes = \%{$vb_Obj->[0]};

    # don't copy the data use %{$genes} within the code e.g. my %genes = %{$genes};

    my $line = q{====================================================================================================};
    print STDOUT qq{\n$line\nCoordinate system check.\n$line\n};

#print STDOUT qq{\n$db has available gene biotypes:\n};
#my @biotypes =  @{(&VectorBase::Gff3::biotypes_list($db_handle))};
#print STDOUT qq{[$_] \x27$biotypes[$_]\x27\n} for (0..$#biotypes);

    my @CSs = @{(&_coordsystems_list($db_handle))};
    my $cs_num = scalar @CSs;

    die qq{\nThis database has no valid Coordinate Sytems} if scalar @CSs == 0;

    print STDOUT qq{\nDatabase \'$db\' has $cs_num coordinate systems:\n\n};

    print STDOUT qq{* name: \x27}, $CSs[$_]->[0], qq{' / rank: '}, $CSs[$_]->[1], qq{' / coord_system_id: '}, 
      $CSs[$_]->[2], qq{\x27\n} for (0..$#CSs); # -> is unessary here but let's be anal

    print STDOUT qq{\nComputing all unique scaffold SEQIDs\n};

    
    #y find all unique landmarks
    my %uniq_landmarks;
    ################################################################################################
    #/ is this due to passing the tripartitite hierachy and not individually accessible - why now?!?
    # for my $i (values %{$self->[0]}) { 
    for my $i (values %{$self->[4]}) { 
    ################################################################################################

        #print map { @{$_},qq{\n} } (values %{$i});
        #use Data::TreeDraw; draw($i); 
        $uniq_landmarks{$i->{landmark}} = 1; 
    
    }
    #delete $uniq_landmarks{''};
    # for my $i (values %{$self->[0]}) { $uniq_landmarks{$i->{landmark}} = 1; }

    # let's make slices for those landmarks - one of the sensitive stages so we do it early. 
    my @uniq_landmarks_orig = (keys %uniq_landmarks);
    my @uniq_landmarks_work = @uniq_landmarks_orig;
    my $uniq_num = scalar @uniq_landmarks_work;

    my %slices;
    my $prefix_removed = 0;

    # next or re-do that's about loop counting only!!! need to buffer the data in other array
    PREFIXES:
    while (1) {

        my @uniq_landmarks_buffer;

        #y print out the unique landmarks found
        print STDOUT qq{\nThere are $uniq_num unique seqids in the Gff3 file:\n\n};
        print STDOUT qq{* \x27$uniq_landmarks_work[$_]\x27\n} for (0..$#uniq_landmarks_work);

        print STDOUT qq{\nChecking for prefixes on SEQIDs\n};

        my $prefixes = &_check_for_prefixes(\@uniq_landmarks_work,5, 11);

        my $prefixes_detected = 1;

        #/ this code is re-used - factor into sub
        if (scalar (keys %{$prefixes}) > 0) {
            print STDOUT qq{\nPossible SEQID prefixes detected:\n};

            print qq{\n* Char length }, 11-$_, qq{; }, (keys %{$prefixes->{11-$_}}) for (0..6);
            print qq{\n};

            #draw($prefixes->{$_});
            #y for now we just want the keys not the values - i.e. the count of times the prefix was encountered
            # for (0..6) { my @Ps = (keys %{$prefixes->{11-$_}}); print qq{\n* Char length }, 11-$_, qq{; @Ps};; }
            #draw($prefixes->{$_});
            #y for now we just want the keys not the values - i.e. the count of times the prefix was encountered
        }
        else {
            print STDOUT qq{\nNo prefixes auto-detected.\n};
            $prefixes_detected = 0;
        }

        print STDOUT qq{\nWARNING: aberant prefixes/suffixes on seqid will prevent }
          . q{you from inserting data into EnsEMBL.};

        #print STDOUT qq{\n\nWould you like to remove a prefix from these seqid names? [Y/n] };

        #chomp (my $ans = <STDIN>);
        #if ($ans =~ /^[Yy]$/) { 

        #if ($ans !~ /^[Nn]$/) { 
        #if ($ans ne q{[Yy]}) {

        #y in place transformation with for not map
        print STDOUT qq{\n\n}, q{Enter the regular-expression to remove (remember to escape appropriate }
        . q{chars [default = 'vectorbase\|'] };

        chomp (my $ans = <STDIN>);  # duh - my $ans = chomp (my $choice = <STDIN>);

        #b don't use in place transformation as we want the original data maintained until a descision is made!
        # if ($ans eq q{}) { s/vectorbase\|// for (@uniq_landmarks); }
        # else { s/$ans// for (@uniq_landmarks); }
        #b we're getting an error as its returning the T/F result and not the - string don't modify $_ with map
        #if ($ans eq q{}) { @uniq_landmarks_buffer = map { s/vectorbase\|// } @uniq_landmarks; }
        #else { @uniq_landmarks_buffer = map { s/$ans// } @uniq_landmarks; }
        #b with map the last statement is the return value
        
        if ($ans eq q{}) { 
            @uniq_landmarks_buffer = 
                map { my $temp = $_; $temp =~ s/vectorbase\|//; $temp } @uniq_landmarks_work; 

            $prefix_removed = q{vectorbase\|}; # keep it out of int context pra agora
        }
        else { 
            @uniq_landmarks_buffer = 
                map { my $temp = $_; $temp =~ s/$ans//; $temp } @uniq_landmarks_work; 
            
            $prefix_removed = $ans;
        }

        print qq{\n\nSubstitution with supplied string resulted in the following SEQIDs:\n\n};
        print STDOUT qq{* \x27$uniq_landmarks_buffer[$_]\x27\n} for (0..$#uniq_landmarks_buffer);
        
        #}

        print STDOUT qq{\nAre these landmark names acceptable? [Y/n] };
        #print STDOUT qq{\nAre these landmark names acceptable? [y/N] };

        chomp ($ans = <STDIN>);

        if ($ans !~ /^[Nn]$/) { 
        #if ($ans =~ /^[yy]$/) { 
        #if ($ans eq q{[yy]}) { last prefixes; }

            @uniq_landmarks_work = @uniq_landmarks_buffer;    
            #b last prefixes; there's no need to have multiple exits - just one at the end - but have multiple restarts
        }
        #  in this context redo vs next is irrelevant - that's about loop counting?!?
        else { redo PREFIXES; }
        #else { next prefixes; }

        if ($prefixes_detected == 1) {

            print &_theredoc_seqid_prefix_removal ;
            #print q{----------------------------------------------------------------------------------------------------};
            chomp ($ans = <STDIN>);

            if ($ans eq q{y} || $ans eq q{Y}) { exit }

        }

        #y generate api slice adaptor
        my $db_slice_adaptor = $self->[2][0]->get_adaptor("Slice") or carp qq{\nCoudln't make a slice adaptor};

        #y see if we can find the first seqid in the list in any of the CSs available
        print STDOUT qq{\nAttempting to autodetect the CS type:\n};

        my @the_CS;
        for my $i (0..$#CSs) {
            my $turn1 = $CSs[$i]->[2];
            my $turn2 = $uniq_landmarks_work[0];
            my $sql = qq{select * from seq_region where coord_system_id=$turn1 and name='$turn2';}; # just try the first of them
            my $db_query = $db_handle->prepare($sql);
            $db_query->execute;
            my $blah = $db_query->fetchall_arrayref;
            #r returns an array reference - so can use existence of ref or its actual value!!!  
            #$CSs[$i]->[0], qq{ coordinate system\n} if ($blah);

            if (scalar @{$blah} > 0) {
                #print STDOUT q{* }, $cs->[0], qq{ coordinate system: Possible\n};
                print STDOUT q{* }, $CSs[$i]->[0], qq{ coordinate system: Possible\n};
                push @the_CS, $CSs[$i]->[0];
            }
            else {
                print STDOUT q{* }, $CSs[$i]->[0], qq{ coordinate system: No\n};
            }
        }    

        #b/ in case of if and else - list example seqids for each CS 
        
        my $cs = q{chromosome};

        if (scalar @the_CS == 0) {
            print STDOUT qq{\nCould not auto-detect CS. This means you probably have some prefix/suffice on your names\n}; 
            
            print STDOUT qq{\nFetching examples of SEQIDs for all coordinate systems for database.\n\n};

            for my $i (@CSs) {
                #use Data::TreeDraw;
                my $cs_id =  $i->[2]; # id
                my $sql = qq{select name from seq_region where coord_system_id = $cs_id limit 5;};
                my $db_query = $db_handle->prepare($sql);
                $db_query->execute;
                my $blah = $db_query->fetchall_arrayref;
                #r returns an array reference - so can use existence of ref or its actual value!!!  
                if (scalar @{$blah} > 0) {
                    print STDOUT qq{Example names from }, $i->[0], qq{ coord_system include: \n\n};
                    print q{* }, $blah->[$_][0], qq{\n} for (0..$#{$blah});
                    print qq{\n};
                }
                else {
                    print $i->[0], qq{ coord_system seems to be empty\n};
                }
            }

            print qq{Would you like to re-check for prefixes? [Y/n] };
            chomp ($ans = <STDIN>);
            exit if ($ans =~ /[Nn]/);
            #y go back to the start of prefix check
            my $line = q{====================================================================================================};
            print qq{\n\n$line\nRestarting prefix checking step.\n$line\n};

            #y recover original landmark data
            @uniq_landmarks_work = @uniq_landmarks_orig;
            redo PREFIXES; 
        }
        elsif (scalar @the_CS == 1) {
            print STDOUT qq{\n* The seqids seems to be from $the_CS[0] coordinate system.\n}; 
            $cs = $the_CS[0];
        }
        else {
            print STDOUT qq{\n* There seem to be several matching CSs please select the appropriate one;};
            print STDOUT qq{[$_] $the_CS[$_]} for (0..$#the_CS);
            CHECK:
            while (1) {

                chomp (my $ans = <STDIN>);
                
                if (!any {$_ eq $ans} @the_CS) {

                    print STDOUT qq{\n\* '$ans\' is not a valid option.\n\n};
                    next CHECK;
                }
                # make it the answer given
                $cs = $ans;
            }
        }

        # swtich if there are any slices that cannot be created and redo while loop based on this
        my $CS_flag = 1;

        my $line = q{====================================================================================================};
        $self->_log_errors(qq{\n$line\nAttempting to generate EnsEMBL slice objects for all seqids:\n$line\n}, 1);
        for my $landmark (@uniq_landmarks_work) {
            #my $slice = $db_slice_adaptor->fetch_by_region($cs, $landmark) or carp qq{\ncould not fetch slice object for $landmark};
            if (my $slice = $db_slice_adaptor->fetch_by_region($cs, $landmark)) {
                $self->_log_errors(qq{* Created slice adaptor for $landmark in $cs coordinate system.}, 1); 
                $slices{$landmark}  = $slice;
            }
            else {
                #y send this just to STDOUT as it iterative!?!
                #carp qq{could not fetch slice object for $landmark in $cs coordinate system\n}; 
                print STDOUT qq{* Could not fetch slice object for $landmark in $cs coordinate system.\n}; 

                #r switch flag if anything fails
                $CS_flag = 0;
            }
        }

        # the default action is redo so we need to force escape...
        last PREFIXES if ($CS_flag == 1);

    } # end PREFIX while loop 

    #/ feed object
    $self->[2][2][0] = \%slices;
    $self->[2][2][1] = $prefix_removed;



    # hand over slices for later usage
    return;
    # return \%slices, $prefix_removed;
}

sub _theredoc_seqid_prefix_removal {

my $line = q{----------------------------------------------------------------------------------------------------};

#y want control over tedious layout here so using where doc

my $message = <<"PREFIX";

$line

Aberrant prefixes may continue to cause commitment problems. Please remove any prefixes from the Gff3 file!

If SEQIDs are variable and of the form: 

    supercont3.1000:1-152674.gff3:vectorbase|supercont3.1000
    supercont3.1000:1-152674.gff3:vectorbase|supercont3.1000 
    supercont4.1000:1-152674.gff3:vectorbase|supercont3.1000

Perhaps use a conservative and verbose regular-expression substitution to remove the prefix:

    > perl -pe 's/^supercont\\d\\.\\d+:\\d+-\\d+.gff3:vectorbase\\|//' filename.gff3

Else if they are very simple such as:

    vectorbase|X
    vectorbase|X 
    vectorbase|X  

Try using (remember to escape special chars):

    > perl -pe 's/^vectorbase\\|//' filename.gff3

$line
Would you like to quit and fix any prefix problems in the Gff3 file now? [N/y]
PREFIX

    return $message;

}

sub _coordsystems_list {
    my $db_handle = shift;
    my $sql = 'select name, rank, coord_system_id from coord_system';
    my $db_query = $db_handle->prepare($sql);
    $db_query->execute;
    # return (&_dbi_unwrap($db_query->fetchall_arrayref));
    return $db_query->fetchall_arrayref;
}

#y///////////////////////// EnsEMBL commitment ////////////////////////////////////////////////////

sub _geneLevel_processing_loop {

    my ($self) = @_;

    my $args = $self->[3][1][0];
    
    my $log_hndl = $self->[3][0][0];

    #f/ this is required to keep the output identical - but we don't change the
    #f/ numbers processed or it'd appear in diff
    my %just_genes = %{$self->[0]};
    my $just_genes = \%just_genes;

    my $slices = $self->[2][2][0];

    #o but that should be identical too (except perhaps by virtue of never assigning
    #o to another var it is never copied its still the same hashref?!? why should it matter?!?
    #my $just_genes = \%{$vb_Obj->[0]};

    my $line = q{====================================================================================================};
    $self->_log_errors(qq{\n$line\nStarting commitment.\n$line\n}, 1);

    #y create analysis object
    my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => 'ensembl',) 
      or die qq{\nCould not instantiate analysis Adaptor}; 

    #y create the parameters to use each time - feed it to the object
    my $ens_params = {
        -ANALYSIS  =>  $analysis,
        -VERSION   =>  $args->{version},
    };

    $self->[3][3][1] = $ens_params;

    my $db_ad = $self->[2][0]; # handle is [2][1] # slices are [2][2][0]
    my $prefix_removed = $self->[2][2][1];
    
    ENSGENE:
    while (my ($gene_id, $gene_hashref) = each %{$just_genes}) {

        #y increment processed genes counter
        $processed_genes++;

        #use Data::TreeDraw; draw($gene_hashref);
        #y fix seqid name in-place removal of prefix
        print qq{\nthe seq_id of the gene }, $gene_hashref->{landmark};
        $gene_hashref->{landmark} =~ s/$prefix_removed// if ($prefix_removed);
        if (!exists $slices->{$gene_hashref->{landmark}}) {
            $self->_log_errors(qq{#Slice.NOSLICE# $gene_id\n}
              . qq{* The slice for $gene_id doesn\x27t exist\nSkipping gene.});
            $gene_skipped++;
            next ENSGENE;
        }

        #y check the gene's sense and make sure its decent
        if ($gene_hashref->{strand} ne  '-' && $gene_hashref->{strand} ne '+') {
            $self->_log_errors(qq{#Strand.ILEGAL# $gene_id\n}
              . qq{* Gene $gene_id has an illegal value for its strand orientation }
              . qq{\x27\x27. Skipping Gene.});
            $gene_skipped++;
            next ENSGENE;
        }

        #f as we may have allowed for gene-level features (that weren't of type gene we need
        #f to remove the names that preface them if they exist e.g. non-unique pseudos have em
        # we don't allow ID less at gene level so just grab its id directlt - i.e. can be for loop
        my $stable_id = $gene_hashref->{attributes}{ID};
        #my $stable_id = $gene_id;

        #y remove possible prefix - i.e. if it didn't have type eq gene and had non-unuqe id
        $stable_id =~ s/(_nonunique_\d*|_noID_childOf_.*)$//;

        my $gene_EnsObj;
        if (
            $gene_EnsObj = Bio::EnsEMBL::Gene->new(
                # bored of passing **** around
                %{$ens_params},
                -START      =>  $gene_hashref->{start}, # really no need the API bins it and calculates it from transcripts
                -END        =>  $gene_hashref->{end},
                -STRAND     => ($gene_hashref->{strand} eq '+' ? 1 : $gene_hashref->{strand} eq '-' ? -1 : 0),
                -SLICE      =>  $slices->{$gene_hashref->{landmark}},
            )
        ) { print $log_hndl qq{* Created Bio::EnsEMBL::Gene object for $stable_id.\n}; }
        else {
            $self->_log_errors(qq{* Cannot create Bio::EnsEMBL::Gene object for gene $stable_id});
            $gene_skipped++;
            next ENSGENE;
        }

        #y set stable id
        $gene_EnsObj->stable_id($stable_id);

        #y copy top level addresses for transcript-level features
        my @trans_hash_refs = @{$gene_hashref->{children}};

        #b junk
        print STDOUT qq{\nGENE: $gene_id/\n} if ($args->{diagnostic}); 

        #w gene-level features were selected on basis of children and no parent
        #g call _transcriptLevel_processing_loop
        my @transcript_biotypes = $self->_transcriptLevel_processing_loop($gene_hashref, \@trans_hash_refs, $gene_EnsObj, $gene_id);

        #b junk
        print qq{\ndiagnostic: @transcript_biotypes}  if ($args->{diagnostic});

        #g perform checks for allowed biotypes at gene level
        my $gene_biotype = $self->_geneLevel_biotype_checks(\@transcript_biotypes, $stable_id);

        #b junk
        print qq{\ngenebiotype: $gene_biotype} if ($args->{diagnostic});

        #y set gene biotype from returned value from transcript-level analysis
        $gene_EnsObj->biotype($gene_biotype);

        $gene_EnsObj->start($gene_hashref->{start});
        $gene_EnsObj->seq_region_start($gene_hashref->{start});
        $gene_EnsObj->end($gene_hashref->{start});
        $gene_EnsObj->seq_region_end($gene_hashref->{start});

        #b junk
        if ($args->{diagnostic}) {
            print qq{\n\nCODE: hashref start: }, $gene_hashref->{start};
            print qq{\nCODE: ensobj start: }, $gene_EnsObj->start, qq{\n\n};
            print qq{\n\nCODE: hashref end: }, $gene_hashref->{end};
            print qq{\nCODE: ensobj end: }, $gene_EnsObj->end, qq{\n\n};
        }

        #b put in message to avoid that horrendous one that happens when you save genes w/o transcripts!?!if
        #b could just check trans_added?!?
#        my $t_num;
#        eval { $t_num = scalar @{$gene_EnsObj->get_all_Transcripts} };
#        if ($t_num == 0) {
#        if ($@) {

        eval {
                if (!$gene_EnsObj->get_all_Transcripts) {
                $self->_log_errors(qq{#Gene.NOTRANSCRIPTS# $gene_id\n}
                . qq{* Gene $gene_id has no transcript-level features to store. Skipping gene.});
                $gene_skipped++;
                next ENSGENE;
            }
        };

        #r/ temp debugging
        #if ($gene_id eq q{PBANKA_010030}) { draw($gene_hashref); }
        #if ($gene_id eq q{PBANKA_010030}) { draw($gene_EnsObj); }
        #if ($gene_id eq q{PBANKA_010030}) { draw($gene_EnsObj->{_transcript_array}[0]); }

        #y save gene if validation option wasn't passed
        if (!$args->{validate}) { #f top level validation - i.e. just don't store the gene after the checks

            eval { $db_ad->get_GeneAdaptor->store($gene_EnsObj); }; #r eval needs ';'

            if ($@) {
                my $line = q{====================================================================================================};
                $self->_log_errors(qq{#Gene.NOSAVE# $gene_id\n$line\nPROBLEM SAVING GENE $gene_id: } . $@ . qq{\n$line});
                $gene_nosave++;
                next ENSGENE; # unecessary its leaving the gene loop now...
            }
            else { 
                print $log_hndl qq{* Saved gene $gene_id.\n};
                $gene_save++;
            }
        }
    } # end of gene loop for integrating code
    return;
} # end of _geneLevel_processing_loop

sub _transcriptLevel_processing_loop {

    my ($self, $gene_hashref, $trans_hashrefs, $gene_EnsObj, $gene_id) = @_;

    my $slices = $self->[2][2][0];
    my $prefix_removed = $self->[2][2][1];
    my $log_hndl = $self->[3][0][0];
    my $error_hndl = $self->[3][0][1];

    my $args = $self->[3][1][0];

    #y declare array to store all transcript biotypes for later gene biotype assignment
    my @transcript_biotypes;

    ENSTRANS:
    for my $tr (@{$trans_hashrefs}) {

    #y increment processed trans counter
    $processed_trans++;

        #y grab the child/transcript's name and HASH ref
        my $trans_id = $tr->[0];
        my $trans_hashref = $tr->[1];
        
        #y assign Gff3 type as transcript biotype - but if mRNA needs to be protein_coding
        my $trans_biotype = $trans_hashref->{type} eq 'mRNA' ? 'protein_coding' : $trans_hashref->{type};

        #b junk
        print STDOUT qq{\nTRANSCRIPT: $trans_id/$trans_biotype (gene }, $gene_EnsObj->stable_id, qq{)\n}
          if ($args->{diagnostic});

        #y push biotype to array for later gene biotype assignment
        push @transcript_biotypes, $trans_biotype;

        #r should this also be for gene_biotype?!? make this an arguement in the config file
        if ($trans_biotype !~ /^lincRNA|miRNA|miRNA_pseudogene|misc_RNA|misc_RNA_pseudogene|
                          Mt_rRNA|Mt_tRNA|Mt_tRNA_pseudogene|ncRNA|
                          processed_transcript|protein_coding|pseudogene|rRNA|
                          rRNA_pseudogene|scRNA_pseudogene|snlRNA|snoRNA|snoRNA_pseudogene|
                          snRNA|snRNA_pseudogene|tRNA|tRNA_pseudogene/x) { 

            $self->_log_errors(qq{#Biotype.UKNOWN# $trans_id\n* unknown biotype $trans_biotype, but allowing.\n}); 
        
        }

        #y removal of prefix from seqid in place
        $trans_hashref->{landmark} =~ s/$prefix_removed// if ($prefix_removed);

        if (!exists $slices->{$trans_hashref->{landmark}}) {
            $self->_log_errors(qq{\n#Slice.NOSLICE# $trans_id\n* The slice doesn\x27t exist\nSkipping gene: $trans_id.});
            $trans_skipped++;
            next ENSTRANS;
        }

        #y unpack id
        my $stable_id = $trans_hashref->{attributes}{ID};
        #my $stable_id = $trans_id;

        #y we DO not allow for ID less gene-level or transcript-level things - thus  we can actuall take the ID from the attributes?!?
        #f as we may have allowed for transcript-level features (that weren't of type mRNA may need to remove the names that preface them 

        #y in case it was something with non-unique id we preserved using a suffix
        $stable_id =~ s/(_nonunique_\d*|_noID_childOf_.*)$//;

        #y create an ensembl transcript object
        my $trans_EnsObj;
        if (
            $trans_EnsObj = Bio::EnsEMBL::Transcript->new(
                %{$self->[3][3][1]},
                -SLICE  =>  $slices->{$trans_hashref->{landmark}},
                -START  => $trans_hashref->{start},
                -END    => $trans_hashref->{end},
            )
        ) { print $log_hndl qq{* Created Bio::EnsEMBL::Transcription object for $stable_id\n}; }
        else {
            $self->_log_errors(qq{* Cannot create Bio::EnsEMBL::Transcription object for transcript $stable_id});
            $trans_skipped++;
            next ENSTRANS;
        }

        #b junk
        print qq{\ntranscript biotype: }, $trans_biotype if ($args->{diagnostic});

        #y store stable id
        $trans_EnsObj->stable_id($stable_id);

        #y store transcript-level type as its biotype 
        #r we should have already checked and forced correction of errors here
        $trans_EnsObj->biotype($trans_biotype);
        
        #y grab all the kiddies - make case insensitive
        my @all_children = map { $_->[1] } @{$trans_hashref->{children}};
        my @CDS_refs = grep { uc($_->{type}) eq q{CDS} } @all_children;
        my @exon_refs = grep { lc($_->{type}) eq q{exon} } @all_children;

        #y unpack gene's strand orientation - this shouldn't bve here
        my $gene_strand = $gene_hashref->{strand}; 

        #y process tripartite structure - i.e. only allows for exons/cds
        # if it has children it can still be any old thing - its only when it has CDS if must be protein_coding 
        if (scalar @all_children > 0) {

            #b junk
            print STDOUT qq{* Transcript-level feature $trans_id has }, scalar @exon_refs, 
            q{ exons and }, scalar @CDS_refs, q{ CDS. (Total number of derived features = }, 
            scalar @all_children, qq{)\n} if ($args->{diagnostic});

            #y quick message about presence of silly CDS/EXON level features
            #$self->_log_errors(qq{* Transcript-level feature $trans_id has excess derivative features (i.e. not exons or CDS).})
            #  if ( $trans_hashref->{type} eq q{mRNA} && ( (scalar @exon_refs + scalar @CDS_refs) != scalar @all_children) );
            print $error_hndl qq{#Transcript.EXCESSFEATURES# $trans_id\n}
              . qq{* Transcript-level feature $trans_id has excess derivative features (i.e. not exons or CDS).\n}
              if ( $trans_hashref->{type} eq q{mRNA} && ( (scalar @exon_refs + scalar @CDS_refs) != scalar @all_children) );

            #y grab any possible UTR exon/CDS level features
            my @utr;
            push @utr, grep { $_->{type} =~ /[Uu][Tt][Rr]/ } @all_children;
            #y grab possible polypeptide exon/CDS level feature
            my @polys;
            push @polys, grep { $_->{type} =~ /polypeptide/i } @all_children;

            #y count ignored utr features - i.e. each loop we add them up
            $ignored_utrs += scalar @utr;
            $ignored_polys += scalar @polys;

            #y are the excess cds/exon-level features pointless UTRs?
            if ((scalar @exon_refs + scalar @CDS_refs + scalar @utr + scalar @polys) != scalar @all_children) {
            #if ((scalar @exon_refs + scalar @CDS_refs + scalar @utr + scalar @polys) != scalar @all_children);
                $self->_log_errors(qq{\n* Excess features are not UTRs or polypeptides.});
                my $others = scalar @exon_refs + scalar @CDS_refs + scalar @utr + scalar @polys;
                my $all = scalar @all_children;
                $ignored_others += ($others - $all);
            }

            #y order exons and CDS - using starting position but can equally be end
            @exon_refs = sort { $a->{start} <=> $b->{start} } @exon_refs;
            @CDS_refs  = sort { $a->{start} <=> $b->{start} } @CDS_refs;

            #y depending on genes strand sense reverse or not reverse sorted cds/exons
            if ($gene_strand eq q{-}) { 
                @exon_refs = reverse @exon_refs;
                @CDS_refs  = reverse @CDS_refs;
            }

            #/ chado etc use minimalistic form of gff3 where exons are re-used thus
            #/ while CDS have a single parent the exons have multiple - one for each
            #/ transcript it appears in
            #b strictly they should be separated by ',' but they are using multiple Parents declarations - better?!?
            
            #/ consequently for shared exons we have already processed the strand so
            #/ when we process the next transcript we get an error

            #y check transcript, exons, CDS have same sense
            STRAND:
            for my $ref ($trans_hashref, @exon_refs, @CDS_refs) { 

                my $i = $ref->{strand};
                #my $i = \$ref->{strand};
                #print qq{here: $i and }, $ref->{attributes}{ID}, qq{\n};
                
                #r make sure we don't attempt to convert the values again
                next STRAND if (exists $ref->{strand_converted});

                # in place checks and modification - these are + or - so string context is better
                if ($i eq $gene_strand) {
                # if (${$i} eq $gene_strand) {

                    #y convert phase to ensembl style - or make false if its shite
                    my $i = $i eq '-' ? -1 : $i eq '+' ? 1 : 0;

                    if (!$i) {
                        $self->_log_errors(qq{#Strand.ILEGAL# $trans_id\n}
                          . qq{* Do not recognise this value of strand for exon. Skipping transcript $trans_id.});
                        $trans_skipped++;
                        next ENSTRANS;
                    }

                    #y its fine so re-pack the scalar value into the data structure - could be else
                    $ref->{strand} = $i;

                    #r flag for already converted
                    $ref->{strand_converted} = 1;
                }
                else {
                    $self->_log_errors(qq{#Strand.DISCORDANT# $trans_id #\n}
                    . qq{* The strand sense of the gene, transcript and exons must } 
                    . qq{be in concordance! skipping transcript $trans_id});
                    $trans_skipped++;
                    next ENSTRANS;
                }
            }

            #y already converted the strand ofr transcript etc rest and checked against this so now change this
            $gene_strand = $gene_strand eq '-' ? -1 : $gene_strand eq '+' ? 1 : 0;

            #g check/correct CDS phase data in files
            $self->_check_cds_phase(\@CDS_refs) if (scalar @CDS_refs > 0);

            #g build transcripts 
            $self->_cds2exon_phases_translation_processing(\@CDS_refs, \@exon_refs, $trans_EnsObj, $gene_strand); 

        }

        my $line = q{====================================================================================================};

        #b just a little namespace protection while i find out who/when internal var is recycled
        {

        #y grab coding sequence
        eval {
            my $seq = $trans_EnsObj->translateable_seq;
            #my $seq = $trans_EnsObj->spliced_seq;

            #y check for ATG start
            if (&VectorBase::Gff3::_check_CodingSeq_starts_with_ATG(\$seq)) { 
                #print STDOUT qq{\nTranscript sequence starts with ATG}; - uber verbose printing?
            } else { 
                $self->_log_errors(qq{#CodingSeq.ATG# $stable_id\n} 
                . qq{* Transcript \x27$stable_id\x27 sequence does not start with ATG. }
                . q{Will allow but address this issue.});
                #$trans_EnsObj->stable_id, qq{\x27\n\tWill allow but address this issue!};
            }

            #y check for in frame stops
            my $trans_stable_id = $trans_EnsObj->stable_id;
            my @list = @{&VectorBase::Gff3::_check_CodingSeq_for_stop_codons(\$seq)};
            my $length = length $seq;

            if ( (scalar @list == 1) && ($list[0][1] != $length) ) { 
            #if ( (scalar @list == 1) &&  (($list[0][1] + 2) != $length) ) { 
            #if ( (scalar @list == 1) &&  (($_->[1] + 2) != $length) ) { 
                #print STDOUT qq{\nLast codon is stop codon in translated sequence for: }, $trans_stable_id, qq{\n};
            }
            elsif (scalar @list > 0)  {
                $self->_log_errors(qq{#CodingSeq.STOP# $stable_id\n}
                  . qq{* In frame stop codons found within coding seqeunce for: \x27$trans_stable_id\x27.});
                #print STDOUT q{Codon: }, $_->[0], q{ at base: }, $_->[1], qq{\n} for @list;
                for (@list) {
                    my $codon = $_->[0];
                    my $base = $_->[1];
                    $self->_log_errors(qq{Codon: \x27$codon\x27 at base $base});
                }
                $self->_log_errors(qq{ Will allow but address this issue!});
            } 
            else { 
                #f here we will pull the next 3 bases and check for stop
                $self->_log_errors(qq{#CodingSeq.NOSTOP# $stable_id\n}
                  . q{* Coding sequence does not terminate with a stop codon. Translation seems }
                  . qq{to be incomplete for \x27$trans_stable_id\x27. Will allow but address this issue.});
            }

        #r/ duh - eval requires ';' 
        }; 

        if ($@) { $self->_log_errors(qq{#CodingSeq.NOTRANSLATABLE# $stable_id\n}
                    . qq{$line\nTranscript $stable_id does not seem to form translatable DNA sequence.\n$line}); }

        #b namespace protection ending
        }

        #y seems good so map it to gene
        eval { $gene_EnsObj->add_Transcript($trans_EnsObj);}; #r eval needs ';'

        if ($@) {
            $self->_log_errors(qq{#Gene.TRANSCRIPTSKIPPED# $stable_id\n}
              . qq{$line\nPROBLEM WITH ADDING TRANSCRIPT $stable_id TO GENE\n$line\n$@\n$line});
            
            #y record event and skip
            $trans_notadded++;
            next ENSTRANS;
        }
        else {
            print $log_hndl qq{* Added transcript $stable_id to gene\n};
            #y count it and go on
            $trans_added++;
        }
    } # end of ENSTRANS-level loop - includes ALL transcripts e.g. miRNA, snRNA too 

    #y return the transcript biotypes for later analysis
    return @transcript_biotypes;
} # end of _transcriptLevel_processing_loop sub

sub _cds2exon_phases_translation_processing {

    my ($self, $CDS_refs, $exon_refs, $trans_EnsObj, $gene_strand) = @_;

    #/ do we pass transcript biotype as an arg or just grab it from the object?!?

    my $slices = $self->[2][2][0];
    my $prefix_removed = $self->[2][2][1];

    my $args = $self->[3][1][0];
    my $log_hndl = $self->[3][0][0];
    my $error_hndl = $self->[3][0][1];

    #y Set it as translation object (defines CDS) of transcript
    my $stable_id = $trans_EnsObj->stable_id;

    #y to allow for merging of API transcript and translation code put in simple flag for CDS part
    my $CDS_loop = scalar @{$CDS_refs} > 0 ? 1 : 0;

    #f/ complete twat handling!
    #y are there any actual exons or CDS - i.e. there should be NO 3o features w/o
    #y (1) exons e.g. pseudogene, miRNA...
    #y (2) exons AND CDS - protein_coding

    #/ API won't tolerate exon-less transcripts - fair enough (shnould be at least one can 
    #/ either go to next transcript here OR can check that $transcript->get_all_Exons() exists!?!

    #f no CDS or exons - this is rubbish to be ignored at this level
    if ( (scalar @{$CDS_refs} == 0) && (scalar @{$exon_refs} == 0) ) { # no meaningful offspring
        print $error_hndl qq{#Transcript.NOEXONCDS# $stable_id\n}
          . qq{* Transcript level feature $stable_id has derivative features BUT they are not either exons CDS }
          . qq{features. Skipping transcript.\n};
        $trans_skipped++;
        next ENSTRANS;
        # return;
    }
    #f GTF files - i.e. have ONLY CDS - no exons to map to CDS - bin these!
    #y is it really GTF and has CDS but no exons?!
    # there are sometimes cases of CDS>0 but exon=0 - CDS>exons will be caught by below but not this case
    elsif (scalar @{$exon_refs} == 0) {
        $self->_log_errors(qq{#ProteinCoding.NOEXONS# $stable_id\n}
          . qq{* Transcript $stable_id has CDS but no exons. }
          . qq{This is not a legal Gff3 gene model. Skipping transcript.});
        $trans_skipped++;
        next ENSTRANS;
    }
    #b lazy - can't be bothered to pass biotype as its alreay set with object
    #f already excluded case of no exons - and no CDS is fine with exons - EXCEPT not with protein_coding!?!
    elsif ( (scalar @{$CDS_refs} == 0) && ($trans_EnsObj = q{protein_coding}) ) {
        $self->_log_errors(qq{#ProteinCoding.NOCDS# $stable_id\n}
          . qq{* Transcript $stable_id has biotype 'protein_coding' but has no CDS. }
          . qq{This is not a legal Gff3 gene model. Skipping transcript.});
        $trans_skipped++;
        next ENSTRANS;
        }

    #else {} # no CDS and we wouldn't be here either way

    #y create our translation object if the CDS_loop is flagged - i.e. there are CDS
    my $cds_EnsObj; # just to let it survive till below
    if ($CDS_loop) {

        #y set stable id of translation
        #r/ Make a descion on this - use same name?!?
        my $translation_stable_id = $stable_id;

        $translation_stable_id =~ s/-R(\w)$/-P$1/;
        #$stable_id .= q{_Trans_Temporary};
        #print qq{\ntranslation stable }, $stable_id;
        
        if (
            $cds_EnsObj = Bio::EnsEMBL::Translation->new( # transcript-translation is 1-to-1
                #y assign analysis and version
                %{$self->[3][3][1]},
                -SLICE  => $trans_EnsObj->slice,
            )
        ) { print $log_hndl qq{* Created Bio::EnsEMBL::Translation object for $translation_stable_id\n}; }
        else {
            $self->_log_errors(qq{#Translation.NOTCREATED# $stable_id\n}
              . qq{* Cannot create Bio::EnsEMBL::Translation object for transcript $stable_id});
            $trans_skipped++;
            next ENSTRANS;
        }

        $cds_EnsObj->stable_id($translation_stable_id);

        #y map translation to transcript - i.e. set transcript_id in translation table
        $trans_EnsObj->translation($cds_EnsObj);
    }
 
    #y keep track of the length of the coding region - i.e. should be len%3=0
    my $sumed_cds_length = 0;

    #y create counter for CDS looping within exon loop
    my $cds_num = 0;

    #y create vars to store important exons for later API assignment
    my $exon_EnsObj_transcription_start_Exon;
    my $exon_EnsObj_transcription_end_Exon;
    my $exon_EnsObj_translation_start_Exon;
    my $exon_EnsObj_translation_end_Exon;

    #y loop through exons - additionally extra-loop for matching them to CDS if they exist
    ENSEXON:
    for my $exon_num (0..$#{$exon_refs}) {

        #y increment processed exon counter
        $processed_exons++;

        #y unpack name - redclared with my so is namespace protected
        my $exon_stable_id = $exon_refs->[$exon_num]->{attributes}{ID};

        #y clean off any suffix we may have applied to allow for feature persistence in a hash
        #f if their exons/CDS names are unique that's nothing to do with this - ensembl stores an autoinc as PK so won't throw exception
        $exon_stable_id =~ s/(_nonunique_\d*|_noID_childOf_.*)$// if $exon_stable_id;
	$exon_stable_id = $stable_id.".".($exon_num+1) unless  $exon_stable_id;
	
        #y create ensembl exon object
        my $exon_EnsObj;
        if (
            $exon_EnsObj = Bio::EnsEMBL::Exon->new(
                %{$self->[3][3][1]},
            ) 
        ) { print $log_hndl qq{* Created Bio::EnsEMBL::Exon object for $exon_stable_id\n}; }
        else {
            $self->_log_errors(qq{#Exon.NOTCREATED# $stable_id\n}
              . qq{* Cannot create Bio::EnsEMBL::Exon object for exon $exon_stable_id. Skipping transcript.});
            $trans_skipped++;
            next ENSTRANS;
        }

        #y store start and end exons
        #f API only lets these be declared after all exons have been declared
        $exon_EnsObj_transcription_start_Exon = $exon_EnsObj if ($exon_num == 0);
        $exon_EnsObj_transcription_end_Exon = $exon_EnsObj if ($exon_num == $#{$exon_refs}); 

        #y modify the landmark in place - if they have prefixes
        #f landmark name modification wasn't in place - modifications took place on a separate list of UNIQ landmarks
        $exon_refs->[$exon_num]->{landmark}  =~ s/$prefix_removed// if ($prefix_removed);

        #y unpack the now cleaned landmark
        my $landmark = $exon_refs->[$exon_num]->{landmark};

        #y check the slice exists - really has already been done (i.e. wouldn't have arrive here if this was a prob
        if (!exists $slices->{$landmark}) {
            $self->_log_errors(qq{#Slice.NOSLICE# $stable_id\n}
              . qq{This is very strange. The slice doesn\x27t exist. Skipping transcript.});
            $trans_skipped++;
            next ENSTRANS;
        }

        $exon_EnsObj->stable_id($exon_stable_id);

        #y set other exon details
        $exon_EnsObj->start($exon_refs->[$exon_num]->{start});
        $exon_EnsObj->end($exon_refs->[$exon_num]->{end});
        $exon_EnsObj->strand($exon_refs->[$exon_num]->{strand});

        #y assign the slice to the exon so that ensembl knows to which seqid to hand the exon
        $exon_EnsObj->slice($slices->{$landmark});

        #y give phases to exons if there are no CDS features - i.e. -1 in ensembl
        if (!$CDS_loop) {
        #$exon_EnsObj->phase(-1) if (!$CDS_loop);
            $exon_EnsObj->phase(-1); 
            $exon_EnsObj->end_phase(-1);
        }

        print qq{\n\n\ncalling add_Exon\n\n\n} if ($args->{diagnostic});

        #y assing the exon to the transcript - i.e. via exon_transcript many-too-many cross-ref table
        eval { $trans_EnsObj->add_Exon($exon_EnsObj); }; #r eval needs ';'
      
        if ($@) {
            my $line = q{====================================================================================================};
            $self->_log_errors(qq{#EXON.NOTADDED# $stable_id\n}
              . qq{$line\nPROBLEM WITH ADDING EXON $stable_id TO TRANSCRIPT (do the exon overlap?!?)\n$line\n$@\n$line});
            $trans_notadded++;
            next ENSTRANS;
        }
        #r put in gene saved to log file?!?
        else { print $log_hndl qq{* Added exon $exon_stable_id to transcript.\n}; }

        #y/ match each CDS to an exon - i.e. if this doesn't match the exon on non-coding - 
        #y/ technically we're matching exons to CDS but there may be CDS-less exons
        if ( 
            $CDS_loop && ($exon_refs->[$exon_num]->{start} <=  $CDS_refs->[$cds_num]->{end}) && # bound cds to left end of exon
            ($CDS_refs->[$cds_num]->{start}   <=  $exon_refs->[$exon_num]->{end}) 
           ) {

            #y increment processed cds counter
            $processed_cds++;

            #y keep check of the sumed CDS length for later
            $sumed_cds_length += ($CDS_refs->[$cds_num]->{end} - $CDS_refs->[$cds_num]->{start} +1 );

            #y check appropriate settings for phase calculations
            my $x = $gene_strand == 1 ? 'start' : 'end';
            my $y = $gene_strand == 1 ? 'end'   : 'start';

            #g check that cds/exons start/stop etc in correct places relative to one another
            $self->_cdsExon_overlap_checks($cds_num, $exon_num, $CDS_refs, $exon_refs, $x, $y, $trans_EnsObj->stable_id)
              if (!$args->{nooverlap});

            #y assign ensembl start of exon phase (just 3 - Gffphase)%3
            my $ensembl_phase = 
              #w if starts don't overlap means its first CDS so starting phase of exon is -1
              $exon_refs->[$exon_num]->{$x} == $CDS_refs->[$cds_num]->{$x} 
              ? (3-$CDS_refs->[$cds_num]->{phase})%3 : -1;
            $exon_EnsObj->phase($ensembl_phase);

            #y assign ensembl end of exon phase
            $exon_EnsObj->end_phase(
                $exon_refs->[$exon_num]->{$y} == $CDS_refs->[$cds_num]->{$y} 
                ? (
                #w if they have the same end and ensembl start phase is -1 then it must be the first CDS ->0
                (($ensembl_phase == -1 ? 0 
                #w the ensembl start is not -1 so we calculate the phase properly (ens_start_phase + length)%3
                : $ensembl_phase)
                + ($CDS_refs->[$cds_num]->{end}-$CDS_refs->[$cds_num]->{start}+1) ) % 3
                    )
                #w they don't end so it so this exon must be partially UTR at this end -> -1
                : -1
            );

            #y store first and last coding exon objects for later
            $exon_EnsObj_translation_start_Exon = $exon_EnsObj if ($cds_num == 0);
            $exon_EnsObj_translation_end_Exon = $exon_EnsObj if ($cds_num == $#{$CDS_refs});

        #b junk
        # TWAT YOU CAN'T PUT THIS CODE AFTER THE INCREMENTER - TAHT MUST BE LAST - just printing for diff
        print STDOUT qq{\n-----------------\n| gene strand = $gene_strand},
              qq{\n-----------------\n| exon iteration $exon_num},
                                 qq{\n| cds iteraction: $cds_num},
              qq{\n-----------------\n| cds phase: [}, $CDS_refs->[$cds_num]->{phase},
                                qq{]\n| ens phase: [}, $ensembl_phase,
                                qq{]\n| ens end phase: }, (
            $exon_refs->[$exon_num]->{$y} == $CDS_refs->[$cds_num]->{$y}
              ? (
              #w if they have the same end and ensembl start phase is -1 then it must be the first CDS ->0
              (($ensembl_phase == -1 ? 0 
              #w the ensembl start is not -1 so we calculate the phase properly (ens_start_phase + length)%3
              : $ensembl_phase)
              + ($CDS_refs->[$cds_num]->{end}-$CDS_refs->[$cds_num]->{start}+1) ) % 3
                )
            #w they don't end so it so this exon must be partially UTR at this end -> -1
              : -1),
              qq{\n-----------------\n} if $args->{diagnostic};

            #y absolute last thing we do CDS-exon matcher is increment CDS for next exon match 
            $cds_num++ if ($cds_num != $#{$CDS_refs}); #w try/do next cds

        } #w end of if clause on overlapping cds/exons 
        else {

            #y doesn't have a matching CDS so its non-coding
            $exon_EnsObj->phase(-1);
            $exon_EnsObj->end_phase(-1);
        }

    } #w end exon loop

    #y now that all exons have been allocated assign start and end exons
    $trans_EnsObj->start_Exon($exon_EnsObj_transcription_start_Exon);
    $trans_EnsObj->end_Exon($exon_EnsObj_transcription_end_Exon);

    #y cds QC and saving EnsTranslation details
    if ($CDS_loop) {

        #y check for clearly strange genes and warn
        if ($sumed_cds_length % 3 != 0) {
            $self->_log_errors(qq{#CodingSeq.INVALID# $stable_id\n}
              . qq{* The modulus 3 of the coding sequence length for transcript $stable_id }
              . q{is not 0. This is not a valid coding sequence. Allowing, but address this issue!});
        }

        #y check for CDS that haven't been matched to a gene
        if ($cds_num != $#{$CDS_refs}) {
            my $mapped_cds_num = $cds_num + 1; # i.e. adjust for 0..-1index
            my $n = scalar @{$CDS_refs};
            $self->_log_errors(qq{#CodingSeq.EXCESSCDS# $stable_id\n}
              . qq{* Unable to allocate all CDS regions to exons for transcript $stable_id }
              . qq{ (mapped $mapped_cds_num CDS out a total of $n. Skipping transcript).});
            $trans_skipped++;
            next ENSTRANS;
        }
        else {

            #y set first translated exon - i.e. first exon to map to CDS (arrays were ordered by strand orientation already)
            $cds_EnsObj->start_Exon($exon_EnsObj_translation_start_Exon);

            #y set translation start site - clearly strand dependent
            $cds_EnsObj->start(_exon_coord(
              $exon_EnsObj_translation_start_Exon, 
              ($gene_strand == 1 ? $CDS_refs->[0]->{start} : $CDS_refs->[0]->{end}))
            );

            #y put in last translated exon
            $cds_EnsObj->end_Exon($exon_EnsObj_translation_end_Exon);

            #y set translation end site
            $cds_EnsObj->end(_exon_coord(
             $exon_EnsObj_translation_end_Exon, 
             ($gene_strand == 1 ? $CDS_refs->[$#{$CDS_refs}]->{end} : $CDS_refs->[$#{$CDS_refs}]->{start}))
            );
        }
    }
    return;
}

#y///////////////////////// data checks for commitment ////////////////////////////////////////////

#y probably want to re-think strategy about what biotype mixes are (not) allowed
sub _geneLevel_biotype_checks {

    my ($self, $list, $gene_id) = @_;

    my $error_hndl = $self->[3][0][1];

    my %distinct;
    @distinct{@{$list}} = (); # assign slice to empty list
   
    my @types = (keys %distinct);
    print qq{\nDistinct transcript level biotypes: @types/$gene_id} if $self->[3][1][0]{diagnostic};

    #f this really needs to be an arg
    my $regexp = qr{lincRNA|miRNA|miRNA_pseudogene|misc_RNA|misc_RNA_pseudogene|Mt_rRNA|Mt_tRNA|Mt_tRNA_pseudogene|ncRNA|processed_transcript|protein_coding|pseudogene|rRNA|rRNA_pseudogene|scRNA_pseudogene|snlRNA|snoRNA|snoRNA_pseudogene|snRNA|snRNA_pseudogene|tRNA|tRNA_pseudogene};

    #y already converted mRNA to protein_coding in transcript_processing
    if (scalar @types == 1) {
        if ($types[0] !~ /^$regexp$/) {
            print $error_hndl qq{#Biotype.UNKNOWN# $gene_id\n}
              . qq{* Unknown biotype $types[0], but allowing.\n}; 
        }
        #else { return $types[0]; } # log message?!?
        return $types[0];
    }
    else {
        _log_errors(qq{#Biotypes.DISCORDANT# $gene_id\n}
          . qq{* Biotypes of transcripts for gene $gene_id are discordant! Skipping gene.});
        $gene_skipped++;
        next ENSGENE;
    }
    return;
}

sub _check_cds_phase {

    my ($self, $CDS_refs) = @_;
    my $phase_minusOne;
    my $length_minusOne;

    my $error_hndl = $self->[3][0][1];

    my @temp;
    CDSPHASE:
    for my $i (0..$#{$CDS_refs}) {

        #o ensembl phase = (length+ensphase)%3 = (length-gffphase)%3
        #o i.e. ens-phase = (3-gffphase)%3
        #w done crudely - but if its the first CDS it has phase=0 else...
        my $phase = $i == 0 ? 0 : ((3 -(($length_minusOne%3)-$phase_minusOne) ) % 3);

        #o already checked for acceptable values of strand so no probs
        
        my $length = $CDS_refs->[$i]->{end} - $CDS_refs->[$i]->{start} +1;

        my @parents = @{$CDS_refs->[$i]->{parents}};
        my $cds_numer = $i+1;
        if ($CDS_refs->[$i]->{phase} eq q{.}) {
        #if ($CDS_refs->[$i]->{phase} eq q{.} || $CDS_refs->[$i]->{phase} eq q{.} ) {

            print $error_hndl qq{#Phase.NONEILEGAL#\n* No phase/Ilegal value for phase data for CDS no. $cds_numer with } 
              . qq{parent(s): @parents. Using generated value = $phase.\n};

              #y put in correct value
              $CDS_refs->[$i]->{phase} = $phase;
              #print STDOUT qq{* No phase/Ilegal value for phase data for CDS no. $cds_numer with parents: },
              #@{$CDS_refs->[$i]->{parents}}, qq{. Using generated value.\n};
        }
        elsif ($CDS_refs->[$i]->{phase} != $phase) {

            print $error_hndl qq{#Phase.ERROR#\n* The phase data for CDS no. $cds_numer with parent(s): }
              . qq{@parents is incorrect. Using corrected value = $phase.\n};

              #y put in correct value
              $CDS_refs->[$i]->{phase} = $phase;
        }

        #y store info for next iteration
        $length_minusOne = $length;
        $phase_minusOne = $phase;
    }
    return;
}

sub _cdsExon_overlap_checks {

    my ($self, $cds_num, $exon_num, $CDS_refs, $exon_refs, $x, $y, $stable_id) = @_;

    my $error_hndl = $self->[3][0][1];

    #/ the first CDS in list is always going to be the first CDS of the gene as we have ordered according to strand
    #o thus is corresponds to $x
    if ( ($cds_num == 0) && ($CDS_refs->[$cds_num]->{$x} == $exon_refs->[$exon_num]->{$x}) ) { 
        # if the starts are the same - there's a missing UTR/promoter
        #$self->_log_errors(qq{#UTR.THREEPRIME# $stable_id\n}
        print $error_hndl qq{#UTR.THREEPRIME# $stable_id\n}
          . qq{* First exon and CDS for transcript \x27$stable_id}
          . q{' overlap - is it missing 3' UTR?};
    }
    elsif ( ($cds_num == $#{$CDS_refs}) && ($CDS_refs->[$cds_num]->{$y} == $exon_refs->[$exon_num]->{$y}) ) { 
        # if the starts are the same - there's a missing UTR/promoter
        #$self->_log_errors(qq{#UTR.FIVEPRIME# $stable_id\n}
        print $error_hndl qq{#UTR.FIVEPRIME# $stable_id\n}
          . qq{* Last exon and CDS for transcript \x27$stable_id}
          . q{' overlap - is it missing a 5' UTR?};
    }
    #o/ as we haven't checked ONLY for first and last we haven't discounted #those as options and need to put them in
    elsif ( 
      #o select intermediate iterations
      ( $cds_num != 0 && $cds_num != $#{$CDS_refs} ) &&
      #o select CDS/exon pairs that don't overlap at start
      ( ( ($CDS_refs->[$cds_num]->{start} != $exon_refs->[$exon_num]->{start}) ) ||
      #o select CDS/exon pairs that don't overlap at end
      ($CDS_refs->[$cds_num]->{end} != $exon_refs->[$exon_num]->{end}) ) ) {
        $self->_log_errors(qq{#CodingSeq.CDSEXON# $stable_id\n}
          . qq{* Internal CDS and exons for transcript \x27$stable_id}
          . q{' do not all overlap - there is a problem with this gene model. Ignoring error.});
    }
    return;
}

sub _check_CodingSeq_for_stop_codons {
    #y SCALAR REF - in case the DNA sequence is long don't pass data directly but by reference
    my $dna_ref = shift;
    my $dna = ${$dna_ref};
    my @stops;
    while ($dna =~ /\G(\w{3})*?([Tt][Aa][Gg]|[Tt][Aa][Aa]|[Tt][Ga][Aa])/g) {push @stops, [uc($2) ,(pos $dna)-2];}
    return \@stops;
}

sub _check_CodingSeq_starts_with_ATG {
    my $dna = ${$_[0]};
    return $dna =~ /^[Aa][Tt][Gg]/ ? 1 : 0;
}

#y///////////////////////// summary output ////////////////////////////////////////////////////////

sub _parsing_summary {
    
    my ($self, $ID_hash,$diag,$log_hndl,$report_limit) = @_;

    my %ID_hashOhash = %{$ID_hash};
    my %diagnostic = %{$diag};

    my $line = q{====================================================================================================};
    $self->_log_errors(qq{\n$line\nAnalysing data-structure.\n$line\n}, 1); # turn setting none-error mode with 1

    $self->_log_errors( sprintf(qq{* Processed %d features.},$diagnostic{feature_count}), 1);
    $self->_log_errors( sprintf(qq{* There were %d features with unique ids.},$diagnostic{count_uniqIDs}), 1);
    $self->_log_errors( sprintf(qq{* There were %d features with non-unique ids.},$diagnostic{_count_nonuniqueID_}), 1);
    $self->_log_errors( sprintf(qq{* There were %d features without IDs but that have parents (will be retained).},
      $diagnostic{_child_noID_}), 1);
    $self->_log_errors( sprintf(qq{* There were %d features without IDs or Parents (were ignored).},
      $diagnostic{_count_noID_noParent_}), 1);

    $diagnostic{total} = $diagnostic{count_uniqIDs}+$diagnostic{_child_noID_};

    $self->_log_errors( sprintf(qq{* We will process a total of %d features.}, 
      $diagnostic{total}), 1);

    $diagnostic{hashOhash_keys} = [keys %ID_hashOhash];
    my $key_number = scalar @{$diagnostic{hashOhash_keys}};

    $self->_log_errors(qq{\nData parsed.\n});
    $self->_log_errors( sprintf (qq{* Perl structure has %d features IDs},$key_number), 1);

    #/ key_number = diag{total}. 
    #b there too few IDs relative to feature lines in the gff3 file
    die qq{\n\nWTF!?!} if ( ($diagnostic{total} > $diagnostic{_count_nonuniqueID_} + $key_number) 
      || ($diagnostic{feature_count} > $diagnostic{_count_nonuniqueID_} + $key_number) );

    #b too many IDs relative to feature line number in gff3 file - i.e. inconsistent naming etc.
    $self->_gff3ID_inconsistency_check(\%diagnostic) 
      if ( ($diagnostic{total} < $key_number)
        || ($diagnostic{feature_count} < $key_number) );

    return;
}

sub _commitment_summary {

    my $self = shift;

    my $line = q{====================================================================================================};
    $self->_log_errors(qq{\n$line\nGene commitment summary.\n$line\n}, 1); # turn setting none-error mode with 1

    #y find childless/parentless features
    my @childless = (keys %{$self->[1]});

    #y So this is compltely unecessary - can do this differently and not store [4]
    $self->_log_errors(sprintf(qq{* %d features make up native perl data structure.},
      &_sum($processed_genes,$processed_trans,$processed_exons,   
        $processed_cds,$ignored_utrs,$ignored_polys,$ignored_others,scalar @childless)), 1); # using regular log file
    
    $self->_log_errors(sprintf(qq{* %d features form part of gene-transcript-cds/exon }
      . q{hierachical relationships.}, &_sum($processed_genes,$processed_trans,
      $processed_exons,$processed_cds,$ignored_utrs)), 1);

    $self->_log_errors(sprintf(qq{* %d orphaned (parentless/childless) features.}, 
      scalar @childless), 1);

    $self->_log_errors(sprintf(qq{* Ignored %d UTR features.}, $ignored_utrs), 1); # [-1] really 11 is safer in case we add to it

    $self->_log_errors(sprintf(qq{* Ignored %d polypeptide features.}, $ignored_polys), 1);

    $self->_log_errors(sprintf(qq{* Other ignored %d features.}, $ignored_others), 1);
    
    $self->_log_errors(sprintf(qq{* Processed %d features for EnsEMBL commitment.}, &_sum($processed_genes,
      $processed_trans,$processed_exons,$processed_cds)), 1);

    $self->_log_errors(sprintf(qq{\n* Processed %d gene level features.\n* Processed %d transcript level features.}
      . qq{\n* Processed %d exon level features.\n* Processed %d cds level features.}, 
      $processed_genes,$processed_trans,$processed_exons,$processed_cds), 1);

    $self->_log_errors(sprintf(qq{\n* Saved %d genes.\n* Failed to save %d genes.\n* Skipped %d genes.}, 
      $gene_save, $gene_nosave, $gene_skipped), 1);

    $self->_log_errors(sprintf(qq{\n* Added %d transcripts.\n* Failed to add %d transcripts.\n* Skipped %d transcripts.\n}, 
      $trans_added, $trans_notadded, $trans_skipped), 1);

    #print qq{\nList of orphaned features: @childless};
    return;
}

#y///////////////////////// misc. /////////////////////////////////////////////////////////////////

sub _log_errors {

    my ($self, $string, $mode) = @_;

    $mode ||= 0;

    die qq{\nnot called _log_errors correctly} if ($mode != 0 && $mode != 1);
    
    # my $args = $self->[3][1][0];

    my $hndl;
    if (!$mode) {
        $hndl = $self->[3][0][1];
    }
    else {
        $hndl = $self->[3][0][0];
    }

    # print to standard out
    print STDOUT qq{$string\n};

    # print to errorlog/log
    print $hndl qq{$string\n};

    return;

}

sub _sum {
    my $s = 0;
    $s += $_ for (@_);
    return $s;
}

sub _dbi_unwrap {
    my $dbi_data = shift;
    return [map { $_->[0] } @{$dbi_data}];
}

sub _exon_coord {
    my ($exon, $coord) = @_;

    if ($exon->strand == 1) {
        my $start = $exon->start;
        return $coord - $start + 1;
    } 
    else {
        my $end = $exon->end;
        return $end - $coord + 1;
    }    
}

sub _find_array_differences {
    my ($a, $b) = @_;
    my @a = @{$a};
    my @b = @{$b};
    my %a = map { $_ => 1 } @a; 
    my %b = map { $_ => 1 } @b;
    my @d = grep( !defined $a{$_}, @b);
    return \@d;
}

sub _check_for_prefixes {
    my ($a_ref, $min, $max) = @_; 
    my $string = join (q{ }, @{$a_ref});
    my %prefixes;
    for my $i ($min..$max) {
        #/ twat forgot to anchor it at boundary
        while ($string =~ /\s(\S{$i})\S/ig) { $prefixes{$i}->{$1}++; } 
    }
    return \%prefixes
}

sub _external_call { 
    my $thing = shift;

    print qq{\ntesting};
    open my $hndl, q{<}, q{rar.pl};
    my $code = do { local $/; <$hndl> };
    close $hndl;
    my $code_ref = eval $code;
    $code_ref->($thing);
    print qq{\ntest over\n};
    return;
}

#=fs map CDS to exons to allow for CDS phase conversion to Ensembl exon phase (i.e. protein coding is done diff)

# check if the cds overlaps at any point with the exon - thus is its
# corresponding cds (need at any point to pick out 3' and 5' with diff start/stop!?!

#f/ check if there is any overlap whatsoever between this exon and the first cds
# this is done by making sure that the exon_end is >= (the = part is vital as it
# fixes the exon as at least overlapping) CDS_start (aka CDS_start <= exon_end)
#
#  X<--  s-EXONEXON-e  ---> can move forward must not any further back and still be true
#                   |
#                   s-CDSCDS-e
#
# exon_start <= than CDS_end aka CDS_end >= exon_start
#
#      <---    s-EXONEXON-e   --->X CANNOT MOVE FURTHER IN THIS DIRECTION AND BE TRUE
#              |
#     s-CDSCDS-e
#
#f/ so for any CDS that has BOTH these requirements TRUE must be bounded by exon overlap
#
#
#    X<---  s-cdscdscdscds-e                   s-cdscdscdscds-e                 `  s-cdscdscdscds-e --->X
#                          |             --->                    ---->             |
#                          s-exonexon-e         s-exonexon-e            s-exonexon-e
#
#o/ the mo-fo is constrained... - doesn't matter which way around you do those conditions
#
#y BUT as its a linear scale l->r generally its best to use < as the elements of
#y the condition have the same order as on the scale

#=fe

1;

#//////////////////////////////////////////////////////////////////////////////////////////////////

package main; # {{{1
use lib '/usr/local/ensembl-live/ensembl/modules/';  #'/home/dsth/ensembl_source/ensembl_56/modules/';
use strict;
use warnings;
use Carp;
use Getopt::Long;
#use Data::TreeDraw;
use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis;

#use Bio::EnsEMBL::Registry;
#print qq{\nthis is API version: }, Bio::EnsEMBL::Registry::software_version;

# options lexicals {{{2
my $db = 'dsth_testing'; #my $host = 'mysql.ebi.ac.uk';
my $host = 'localhost';
my $port = '4157';
my $db_id = qq{DBI:mysql:database=$db;host=$host;port=$port}; #my $user = 'anonymous';
my $user = 'dsth';
my $pw = q{};

my %arguments;
$arguments{log} = './log_gff3';
$arguments{gff3file} = 'Anopheles_gambiae_PEST_AgamP3.5.gff3-truncate';
$arguments{report_limit} = 100;
$arguments{noforceunique} = 0;
$arguments{validate} = 0;
$arguments{version} = 1;          
$arguments{nooverlap} = 0;

my $verbose = 1;          
my $external = 0;          
my $silent = 0;
my $source = 'gff3';          
my $gene_prefixes = 1; 
my $analyse_comments = 0;
my $nosave = 0;
my $cs_name = q{};

$arguments{diagnostic} = 0; # turn off all the crap 
open my $diag, '>', 'diagfile' or die if ($arguments{diagnostic});

#y get options #{{{2
&GetOptions( 
    # options
    '-noforce-unique'   => \$arguments{noforceunique},  
    '-notypecheck'      => \$analyse_comments,
    '-ignoregeneprefix' => \$gene_prefixes,
    '-validate'         => \$arguments{validate},
    # args
    '-coordsystem:s'    => \$cs_name,
    '-dbname:s'         => \$db,
    '-host:s'           => \$host,
    '-port:i'           => \$port,
    '-log:s'            => \$arguments{log},
    '-user:s'           => \$user,
    '-password:s'       => \$pw,
    '-version:s'        => \$arguments{version},
    '-source:s'         => \$source,
    '-file:s'           => \$arguments{gff3file},   # :s means its a string - :i is integer - presumably dies otherwise
    '-verbose'          => \$verbose, 
    '-limit:i'          => \$arguments{report_limit}, 
);

#y the actual method calls {{{2
my $vb_Obj = VectorBase::Gff3->new($db, $host, $port, $user, $pw, \%arguments);

#y run the parser
$vb_Obj->_parseGff3File();

##### here's the fucking place they get isolated!?!
#y separate gene-strucutres and isolated features into [0] and [1]
$vb_Obj->_separate_tripartite_and_singleFeatures();
#my %genes; #my %genes = {}; # idit - adicione uma coisa vazia

$vb_Obj->_gene_prefix_detection if !$gene_prefixes;

$vb_Obj->_analyse_comments_column if !$analyse_comments;

#y check the seqids against CSs and grab slices - needs to be broken into MANY subs
if ($cs_name) { $vb_Obj->_fetch_slices_quick($cs_name); }
else { $vb_Obj->_fetch_slices_with_checks($db); }

#=fs Ad-hoc manipulation
#&external::_external(\%just_genes);
&external($vb_Obj->[0]) if $external;
#&external::_external(\%{$vb_Obj->[0]});
#sub {                                                                                
# print qq{\ntesting};     
# my $genelevel_href = shift;
# print qq{\ngene name: }, $_ for (keys %{$genelevel_href});
# print qq{\ntest over\n};
# };
#=fe

draw($vb_Obj->[0]{'vectorbase|AGAP000017'}) if ($arguments{diagnostic});

$vb_Obj->_geneLevel_processing_loop();

$vb_Obj->_commitment_summary if !$silent;

1;

