#!/usr/local/bin/perl
###############################################################################
#
#   Name:           EnsEMBL::HTML::StaticTemplates.pm
#
#   Description:    Populates templates for static content.
#                   Run at server startup
#
###############################################################################

package EnsEMBL::Web::Document::ColourChanger;
use strict;
use warnings;
no warnings "uninitialized";

use GD;
use GD::Text;
use File::Path;

my $TRIANGLE_Y_OFFSET = 3;
my $TRIANGLE_X_OFFSET = 3;
my $TRIANGLE_WIDTH    = 6;
my $TRIANGLE_HEIGHT   = 6;

#----------------------------------------------------------------------
sub change_ddmenu_colours {
    my $species_defs = shift;
    my $colours      = $species_defs->ENSEMBL_STYLE;

    # my $img_directory = $species_defs->ENSEMBL_SERVERROOT . '/htdocs/img';
    my $img_directory = $SiteDefs::ENSEMBL_HTDOCS_DIRS[0] . '/img';
    my %colours
        = qw(170 BACKGROUND1 187 BACKGROUND2 204 BACKGROUND3 221 BACKGROUND4 238 BACKGROUND5);
    my $replace = {
        map {
            $_ =>
                join('', map { chr($_) } hex2rgb($colours->{ $colours{$_} }))
            } keys %colours
    };
    if (opendir DH, "$img_directory/templates") {
        while (my $file = readdir DH) {
            next if $file =~ /zoom\d\.gif/;
            next unless $file =~ /\.gif/;
            local $/ = undef;
            open I, "$img_directory/templates/$file";
            my $X = <I>;
            close I;
            my $flags = ord(substr($X, 10, 1));
            next unless $flags & 128;    # No global colourmap arg!!
            my $colourcount = 2 << ($flags & 7);

            foreach my $N (1 .. $colourcount) {
                my $index  = $N * 3 + 10;
                my $COL    = ord(substr($X, $index, 1));
                my $NEWCOL = $replace->{$COL};
                substr($X, $index, 3) = $NEWCOL if $NEWCOL;
            }
            my $dir = $file =~ /^y/ ? 'dd_menus/' : '';
            my $dd_menu_path = "$img_directory/$dir";
            if (! -d $dd_menu_path) {
                File::Path::mkpath($dd_menu_path);
            }
            if (open O, ">$dd_menu_path$file") {
                print O $X;
                close O;
            } else {
                warn "Cannot write to '$dd_menu_path$file': $!";
            }
        }
    }

    # createGraphicsButtons($species_defs);
    generate_buttons($species_defs, "$img_directory/dd_menus");
}

sub hex2rgb {
    my $hex = shift;
    my ($r, $g, $b) = $hex =~ /(..)(..)(..)/;
    if (defined($r)) {
        return (hex($r), hex($g), hex($b));
    }
}

sub change_zoom_colours {
    my $species_defs  = shift;
    my $colours       = $species_defs->ENSEMBL_STYLE;
    my $img_directory = $species_defs->ENSEMBL_SERVERROOT . '/htdocs/img';
    my @C             = hex2rgb($colours->{'HEADING'});
    my @O             = hex2rgb($colours->{'SPECIESNAME'});
    my @B             = hex2rgb($colours->{'BACKGROUND5'});
    foreach my $i (1 .. 8) {
        local $/ = undef;
        open I, "$img_directory/templates/zoom$i.gif" || next;
        my $Y = my $X = <I>;
        close I;
        my $flags = ord(substr($X, 10, 1));
        next unless $flags & 128;    # No global colourmap arg!!
        my $colourcount = 2 << ($flags & 7);
        my $left        = 80;
        my $right       = 255 - $left;

        foreach my $N (1 .. $colourcount) {
            my $index = $N * 3 + 10;
            my $COL   = (
                ord(substr($X, $index, 1)) + ord(substr($X, $index + 1, 1)) +
                    ord(substr($X, $index + 2, 1))) / 3;
            foreach my $P (0 .. 2) {
                my $N_COL =
                      $COL < $left
                    ? $COL * $C[$P] / $left
                    : $C[$P] + ($COL - $left) / $right * ($B[$P] - $C[$P]);
                substr($X, $index + $P, 1) = chr($N_COL);
            }
        }
        if (open O, ">$img_directory/buttons/zoom$i.gif") {
            print O $X;
            close O;
        }
        foreach my $N (1 .. $colourcount) {
            my $index = $N * 3 + 10;
            my $COL   = (
                ord(substr($Y, $index, 1)) + ord(substr($Y, $index + 1, 1)) +
                    ord(substr($Y, $index + 2, 1))) / 3;
            foreach my $P (0 .. 2) {
                my $N_COL =
                      $COL < $left
                    ? $COL * $C[$P] / $left
                    : $O[$P] + ($COL - $left) / $right * ($B[$P] - $O[$P]);
                substr($Y, $index + $P, 1) = chr($N_COL);
            }
        }
        if (open O, ">$img_directory/buttons/zoom${i}on.gif") {
            print O $Y;
            close O;
        }
        close I;
    }

}

sub change_CSS_colours {
    my $species_defs = shift;
print STDERR "Changing CSS Colours\n";
    local $/ = undef;
    my %colours = %{ $species_defs->ENSEMBL_STYLE || {} };
    foreach (keys %colours) {
        $colours{$_} =~ s/^([0-9A-F]{6})$/#$1/i;
    }
    my $css_directory = $species_defs->ENSEMBL_SERVERROOT . '/htdocs/css';
    if (opendir DH, $css_directory) {
        while (my $file = readdir DH) {
            if ($file =~ /^(.*)-tmpl$/) {
                open I, "$css_directory/$file" || (
                    warn "Cannot read from '$css_directory/$file'" && next);
                if (open O, ">$css_directory/$1") {
                    local ($/) = undef;
                    my $T = <I>;
                    $T
                        =~ s/\[\[(\w+)\]\]/$colours{$1}||"\/* ARG MISSING DEFINITION $1 *\/"/eg;
                    print O $T;
                    close O;
                } else {
                    warn "Cannot write to '$css_directory/$file': $!";
                }
                close I;
            }
        }
    }
}

sub createGraphicsButtons {
    my ($SPECIES_DEFS) = @_;
    my $species = $SiteDefs::ENSEMBL_PERL_SPECIES;
    my $colours = $SPECIES_DEFS->get_config($species, 'ENSEMBL_COLOURS');
    my $map0 = quotemeta(pack('H6', 'ffdf27'));
    my $col0 = pack('H6', $colours->{'background0'});
    my $map1 = quotemeta(pack('H6', 'ffffe7'));
    my $col1 = pack('H6', $colours->{'background1'});
    my $map2 = quotemeta(pack('H6', 'ffffdd'));
    my $col2 = pack('H6', $colours->{'background2'});
    my $map3 = quotemeta(pack('H6', 'ffffcc'));
    my $col3 = pack('H6', $colours->{'background3'});
    my $map4 = quotemeta(pack('H6', '999900'));
    my $col4 = pack('H6', $colours->{'background4'});

    foreach my $dir (qw(buttons dd_menus)) {
        my $T = $SiteDefs::ENSEMBL_HTDOCS_DIRS[0] . "/templates/$dir";
        my $D = $SiteDefs::ENSEMBL_HTDOCS_DIRS[0] . "/img/$dir";
        opendir(DIR, "$T") or next;
        foreach (readdir(DIR)) {
            next unless /\.gif$/;
            my $F = "$T/$_";
            print STDERR "FILE: $F\n";
            next unless -f $F;
            my $l = -s $F;
            open I, $F or (warn $! && next);
            my $gif;
            read I, $gif, $l;
            $gif =~ s/$map1/$col0/;
            $gif =~ s/$map2/$col1/;
            $gif =~ s/$map3/$col2/;
            $gif =~ s/$map4/$col3/;
            $gif =~ s/$map0/$col4/;
            open(O, ">$D/$_") || (warn $! && next);
            binmode(O);
            print O $gif;
            close O;
        }
    }
}

sub generate_buttons {
    my ($species_defs, $directory) = @_;
    my %BUTTONS = (
        'features'       => 'Features',
        'options'        => 'Options',
        'gss_features'   => 'GSSs',
        'est_features'   => 'ESTs',
        'fst_features'   => 'FSTs',
        'imgsize'        => 'Image',
        'array_features' => 'Arrays',
        'markers'        => 'Markers',
        'help'           => 'Help',
        'repeats'        => 'Repeats',
        'mdr_features'   => 'MDRs',
        'dassource'      => 'DAS',
        'exportas'       => 'Export',
        'compara'        => 'Comparative',
        'gene_options'   => 'Genes',
    );
    my $HORIZONTAL_PAD = 5;
    my $VERTICAL_PAD   = 6;
    my $RAISED_X_OFFSET = 2;
    my $RAISED_Y_OFFSET = 2;

    my $style = $species_defs->ENSEMBL_STYLE;

    my $fontname
        = "$style->{'GRAPHIC_TTF_PATH'}/$style->{'GRAPHIC_FONT'}.ttf";
    my $fontsize = $style->{'GRAPHIC_FONTSIZE'};

    for my $key (keys %BUTTONS) {
print STDERR "Creating button for $key\n";

        my $image_filename = "$directory/y-$key.gif";
        my $raised_filename = "$directory/y-$key-o.gif";
        my $gd_text  = GD::Text->new();
        $gd_text->set_font($fontname, $fontsize);
        $gd_text->set_text($BUTTONS{$key});

        my ($width, $height) = $gd_text->get('width', 'height');
        my $image_width  = $width + 3 * $HORIZONTAL_PAD + $TRIANGLE_WIDTH + 1;
        my $image_height = $height + 2 * $VERTICAL_PAD + 1;
        my $image = GD::Image->new($image_width, $image_height);
        my $background4
            = $image->colorAllocate(hex2rgb($style->{'BACKGROUND4'}));
        my $background1
            = $image->colorAllocate(hex2rgb($style->{'BACKGROUND1'}));
        my $background2
            = $image->colorAllocate(hex2rgb($style->{'BACKGROUND2'}));
        my $white = $image->colorAllocate(255, 255, 255);
        my $black = $image->colorAllocate(0,   0,   0);
        $image->transparent($white);
        $image->interlaced('true');

        $image->stringTTF(
            $black, $fontname, $fontsize, 0,
            $HORIZONTAL_PAD + 1,
            int($VERTICAL_PAD / 2) + $height + 1,
            $BUTTONS{$key}
        );

        $image->setAntiAliased($black);
        $image->filledPolygon(
            pulldown_arrow($HORIZONTAL_PAD + $width, $VERTICAL_PAD), $black);

        $image->line(0, $image_height - 2,
            $image_width, $image_height - 2, $background1);

        save_image($image, $image_filename);
        
        $image->line(
                $RAISED_X_OFFSET, $RAISED_Y_OFFSET,
                $image_width - $RAISED_X_OFFSET,
                $RAISED_Y_OFFSET, $background2);
        $image->line($RAISED_X_OFFSET, $RAISED_Y_OFFSET,
                $RAISED_X_OFFSET,
                $image_height - 2*$RAISED_X_OFFSET, $background2);
        $image->line(
                $RAISED_X_OFFSET, $image_height - 2*$RAISED_Y_OFFSET,
                $image_width - $RAISED_X_OFFSET,
                $image_height - 2*$RAISED_Y_OFFSET, $background1);
        $image->line($image_width - $RAISED_X_OFFSET, $RAISED_Y_OFFSET,
                $image_width - $RAISED_X_OFFSET,
                $image_height - 2*$RAISED_X_OFFSET, $background1);
        
        save_image($image, $raised_filename);
    }
}

sub save_image {
    my ($image, $filename) = @_;
    open(my $fh, '>', $filename) || (print(STDERR "Cannot save $filename: $!\n") && return);
    binmode $fh;
    print $fh $image->gif;
    close $fh;
}


sub pulldown_arrow {
    my ($x_start, $y_start) = @_;

    $x_start += $TRIANGLE_X_OFFSET;
    $y_start += $TRIANGLE_Y_OFFSET;

    my $triangle = GD::Polygon->new();
    $triangle->addPt($x_start, $y_start);
    $triangle->addPt($x_start + $TRIANGLE_WIDTH, $y_start);
    $triangle->addPt($x_start + $TRIANGLE_WIDTH / 2,
        $y_start + $TRIANGLE_HEIGHT);

    return $triangle;
}

1;
