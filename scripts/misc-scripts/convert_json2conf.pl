#!/bin/env perl

use JSON::Parse ':all';

my @init_keys = qw(
	menu_key 
	menu_name
	submenu_key
	submenu_name	
	source_name
	description
	source_url
	source_type
	maxHeightPixels
	display
);

my %init_config = (
	menu_key => 'maizecode',
	menu_name => 'MaizeCODE',
	submenu_key => 'RNAseq',
	submenu_name => 'RNAseq',
	source_type => 'bigwig',
	maxHeightPixels => 50,
	display => 'off',
);

my %jkeys_mapping = qw(
	key	source_name
	urlTemplate	source_url
	onClick		description
);

my $json_file = shift @ARGV;
warn ("json file is $json_file\n");

my $p = json_file_to_perl ($json_file);

my @track_names;
for my $atrack ( @{$p->{tracks}}) {

	print "\n";
	
	#map{ print "$_ => $atrack->{$_}\n";} keys %$atrack;

	map{ $init_config{$jkeys_mapping{$_}} = $atrack->{$_} } keys %jkeys_mapping;

	next unless ( $init_config{source_name} && $init_config{source_url} =~ /^http/  );

	$init_config{description} = $init_config{source_name}. sprintf(" see <a href='%s'>workflow</a>", $init_config{description}->{url});
warn("description=$init_config{description}\n");

	#for my $track_key (keys %$atrack) {
	#	if( $track_key eq 'key'){
	#		$init_config{source_name} = $atrack->{$track_key};
	#		}elsif( $track_key eq 'onClick' ){
	#		$init_config{description} = $atrack->{$track_key}->{url};
	#	}elsif( $track_key eq 'urlTemplate' ){
	#		$init_config{description} = $atrack->{$track_key};
	#	}
	#}

	my $track_name = $init_config{source_name} ;
	$track_name =~ s/\s+/_/g;

	printf "[%s]\n", $track_name;
	push @track_names, $track_name;
	
	for my $ik (@init_keys){
		print "$ik = $init_config{$ik}\n";
	}

	
}


print "\n";
map{ print "$_ = functional\n" } @track_names;

__END__


menu_key = maizecode
menu_name = MaizeCODE
submenu_key = RNAseq
submenu_name = RNAseq
source_name = leaf_lower_whorl_vegetative_rep1 forward (RNAseq)
description = leaf_lower_whorl_vegetative_rep1 forward (RNAseq), see <a href='https://www.sciapps.org/?wf_id=7c678f53-613e-4750-9435-4c9459d1fd36'>workflow</a>
source_url = https://data.cyverse.org/dav-anon/iplant/home/maizecode/sci_data/results/MCrna-0.0.1_b172adf1-59b7-40f8-b649-4194730511d3/sig_leaf_lower_whorl_vegetative_rep1_R1.bw
source_type = bigwig
maxHeightPixels = 50
display = off

