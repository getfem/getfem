#!/usr/bin/perl

my $RCS_Id = '$Id: eps2png.pl,v 2.6 2008/03/27 15:07:11 jv Exp $ ';

# Author          : Johan Vromans
# Created On      : Tue Sep 15 15:59:04 1992
# Last Modified By: Johan Vromans
# Last Modified On: Thu Mar 27 16:06:22 2008
# Update Count    : 165
# Status          : Okay

################ Common stuff ################

use strict;
use Getopt::Long 2.1;

my $my_package = "Sciurix";
my ($my_name, $my_version) = $RCS_Id =~ /: (.+).pl,v ([\d.]+)/;
$my_version .= '*' if length('$Locker:  $ ') > 12;

use vars qw($VERSION);
( $VERSION ) = '$Revision: 2.6 $ ' =~ /\$Revision:\s+([^\s]+)/;

################ Program parameters ################

# Some GhostScript programs can produce GIF directly.
# If not, we need the PBM package for the conversion.
# NOTE: This will be changed upon install.
my $use_pbm = 0;

my $res = 82;			# default resolution
my $scale = 1;			# default scaling
my $mono = 0;			# produce BW images if non-zero
my $format;			# output format
my $gs_format;			# GS output type
my $output;			# output, defaults to STDOUT
my $antialias = 4;              # antialiasing
my $width;			# desired widht
my $height;			# desired height

my ($verbose,$trace,$test,$debug) = (0,0,0,0);
handle_options ();
unless ( defined $format ) {
    if ( $0 =~ /2(gif|jpg|png)$/ ) {
	set_out_type ($1);
    }
    else {
	set_out_type ('png') unless defined $format;
    }
}
print STDERR ("Producing $format ($gs_format) image.\n") if $verbose;

$trace |= $test | $debug;
$verbose |= $trace;

################ Presets ################

################ The Process ################

my $eps_file;
my $err = 0;

foreach $eps_file ( @ARGV ) {

    unless ( open (EPS, $eps_file) ) {
	print STDERR ("Cannot open $eps_file [$!], skipped\n");
	$err++;
	next;
    }

    my $line = <EPS>;
    unless ( $line =~ /^%!PS-Adobe.*EPSF-/ ) {
	print STDERR ("Not EPS file: $eps_file, skipped\n");
	$err++;
	next;
    }

    my $ps = "";		# PostScript input data
    my $xscale;
    my $yscale;
    my $gotbb;

    # Prevent derived values from propagating.
    my $width = $width;
    my $height = $height;

    while ( $line = <EPS> ) {

	# Search for BoundingBox.
	if ( $line =~ /^%%BoundingBox:\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/i ) {
	    $gotbb++;
	    print STDERR ("$eps_file: x0=$1, y0=$2, w=", $3-$1, ", h=", $4-$2)
		if $verbose;

	    if ( defined $width ) {
		$res = 72;
		$xscale = $width / ($3 - $1);
		if ( defined $height ) {
		    $yscale = $height / ($4 - $2);
		}
		else {
		    $yscale = $xscale;
		    $height = ($4 - $2) * $yscale;
		}
	    }
	    elsif ( defined $height ) {
		$res = 72;
		$yscale = $height / ($4 - $2);
		if ( defined $width ) {
		    $xscale = $width / ($3 - $1);
		}
		else {
		    $xscale = $yscale;
		    $width = ($3 - $1) * $xscale;
		}
	    }
	    unless ( defined $xscale ) {
		$xscale = $yscale = $scale;
		# Calculate actual width.
		$width  = $3 - $1;
		$height = $4 - $2;
		# Normal PostScript resolution is 72.
		$width  *= $res/72 * $xscale;
		$height *= $res/72 * $yscale;
		# Round up.
		$width  = int ($width + 0.5) + 1;
		$height = int ($height + 0.5) + 1;
	    }
	    print STDERR (", width=$width, height=$height\n") if $verbose;

	    # Scale.
	    $ps .= "$xscale $yscale scale\n"
	      if $xscale != 1 || $yscale != 1;

	    # Create PostScript code to translate coordinates.
	    $ps .= (0-$1) . " " . (0-$2) . " translate\n"
	      unless $1 == 0 && $2 == 0;

	    # Include the image, show and quit.
	    $ps .= "($eps_file) run\n".
	      "showpage\n".
		"quit\n";

	    last;
	}
	elsif ( $line =~ /^%%EndComments/i ) {
	    last;
	}
    }
    close (EPS);

    unless ( $gotbb ) {
	print STDERR ("No bounding box in $eps_file\n");
	$err++;
	return;
    }

    my $out_file;		# output file
    my $pbm_file;		# temporary file for PBM conversion

    # Note the temporary PBM file is created where the output file is
    # located, since that will guarantee accessibility (and a valid
    # filename).
    if ( defined $output ) {
	$out_file = $output;
	$pbm_file = $output.".ppm";
    }
    elsif ( $eps_file =~ /^(.+).epsf?$/i ) {
	$out_file = "$1.$format";
	$pbm_file = $1.".ppm";
    }
    else {
	$out_file = $eps_file . ".$format";
	$pbm_file = $eps_file . ".ppm";
    }
    print STDERR ("Creating $out_file\n") if $verbose;

    my $gs0 = "gs -q -dNOPAUSE -r$res -g${width}x$height";
    my $gs1 = "-";
    $gs0 .= " -dTextAlphaBits=$antialias -dGraphicsAlphaBits=$antialias"
      if $antialias;
    if ( $format eq 'png' ) {
	mysystem ("$gs0 -sDEVICE=". ($mono ? "pngmono" : $gs_format).
		  " -sOutputFile=$out_file $gs1", $ps);
    }
    elsif ( $format eq 'jpg' ) {
	mysystem ("$gs0 -sDEVICE=". ($mono ? "jpeggray" : $gs_format).
		  " -sOutputFile=$out_file $gs1", $ps);
    }
    elsif ( $format eq 'gif' ) {
	if ( $use_pbm ) {
	    # Convert to PPM and use some of the PBM converters.
	    mysystem ("$gs0 -sDEVICE=". ($mono ? "pbm" : "ppm").
		      " -sOutputFile=$pbm_file $gs1", $ps);
	    # mysystem ("pnmcrop $pbm_file | ppmtogif > $out_file");
	    mysystem ("ppmtogif $pbm_file > $out_file");
	    unlink ($pbm_file);
	}
	else {
	    # GhostScript has GIF drivers built-in.
	    mysystem ("$gs0 -sDEVICE=". ($mono ? "gifmono" : "gif8").
		      " -sOutputFile=$out_file $gs1", $ps);
	}
    }
    else {
	print STDERR ("ASSERT ERROR: Unhandled output type: $format\n");
	exit (1);
    }

    unless ( -s $out_file ) {
	print STDERR ("Problem creating $out_file for $eps_file\n");
	$err++;
    }

}

exit 1 if $err;

################ Subroutines ################

sub mysystem {
    my ($cmd, $data) = @_;
    print STDERR ("+ $cmd\n") if $trace;
    if ( $data ) {
	if ( $trace ) {
	    my $dp = ">> " . $data;
	    $dp =~ s/\n(.)/\n>> $1/g;
	    print STDERR ("$dp");
	}
	open (CMD, "|$cmd") or die ("cmd: $!\n");
	print CMD $data;
	close CMD or die ("cmd close: $!\n");
    }
    else {
	system ($cmd);
    }
}

sub set_out_type {
    my ($opt) = lc (shift (@_));
    if ( $opt =~ /^png(mono|gray|16|256|16m|alpha)?$/ ) {
	$format = 'png';
	$gs_format = $format.(defined $1 ? $1 : '16m');
    }
    elsif ( $opt =~ /^gif(mono)?$/ ) {
	$format = 'gif';
	$gs_format = $format.(defined $1 ? $1 : '');
    }
    elsif ( $opt =~ /^(jpg|jpeg)(gray)?$/ ) {
	$format = 'jpg';
	$gs_format = 'jpeg'.(defined $2 ? $2 : '');
    }
    else {
	print STDERR ("ASSERT ERROR: Invalid value to set_out_type: $opt\n");
	exit (1);
    }
}

sub handle_options {
    my  ($help) = 0;		# handled locally
    my ($ident) = 0;		# handled locally

    # Process options.
    if ( @ARGV > 0 && $ARGV[0] =~ /^[-+]/ ) {
	usage () 
	  unless GetOptions ('ident'	   => \$ident,
			     'verbose'	   => \$verbose,
			     'antialias|aa=i'   => \$antialias,
			     'noantialias|noaa' => sub { $antialias = 0 },
			     'scale=f'     => \$scale,
			     'width=i'	   => \$width,
			     'height=i'	   => \$height,
			     'output=s'    => \$output,
			     'png'	   => \&set_out_type,
			     'pngmono'	   => \&set_out_type,
			     'pnggray'	   => \&set_out_type,
			     'png16'	   => \&set_out_type,
			     'png256'	   => \&set_out_type,
			     'png16m'	   => \&set_out_type,
			     'pngalpha'	   => \&set_out_type,
			     'jpg'	   => \&set_out_type,
			     'jpggray'	   => \&set_out_type,
			     'jpeg'	   => \&set_out_type,
			     'jpeggray'	   => \&set_out_type,
			     'gif'	   => \&set_out_type,
			     'gifmono'	   => \&set_out_type,
			     'mono!'	   => \$mono,
			     'resolution=i' => \$res,
			     'pbm!'	   => \$use_pbm,
			     'trace'	   => \$trace,
			     'help'	   => \$help,
			     'debug'	   => \$debug)
	    && !$help;
    }
    print STDERR ("This is $my_package [$my_name $my_version]\n")
	if $ident;
    die ("Only one file argument is allowed when -output is used\n")
      if @ARGV > 1 && defined $output;
    die ("At least one input file name must be specified\n")
      unless @ARGV;
    die ("Antialias value must be 0, 1, 2, 4, or 8\n")
      unless "$antialias" =~ /^[01248]$/;
}

sub usage {
    print STDERR <<EndOfUsage;
This is $my_package [$my_name $my_version]
Usage: $0 [options] file [...]

    -png -pngmono -pnggray -png16 -png256 -png16m -pngalpha
                        produce PNG image
    -jpg -jpggray -jpeg -jpeggray
                        produce JPG image
    -gif -gifmono       produce GIF image
    -[no]mono		monochrome/colour rendition
    -width XXX		desired with
    -height XXX		desired height
    -resolution XXX	resolution (default = $res)
    -scale XXX		scaling factor
    -antialias XX	antialias factor (must be 0, 1, 2, 4 or 8; default: 4)
    -noantialias	no antialiasing (same as -antialias 0)
    -[no]pbm		GIF only: [do not] convert via pbm format
    -output XXX		output to this file (only one input file)
    -help		this message
    -ident		show identification
    -verbose		verbose information
EndOfUsage
    exit 1;
}

# For install testing
1;

__END__

=pod

=head1 NAME

eps2png - convert EPS files to PNG, JPG or GIF images

=head1 SYNOPSIS

    eps2png [ options ] files ...
    eps2gif [ options ] files ...
    eps2jpg [ options ] files ...

=head1 DESCRIPTION

Converts files from EPS format (Encapsulated PostScript) to some
popular image formats.

If installed as C<eps2png> (the default), it produces PNG images by
default. Likewise, C<eps2gif> defaults to GIF images and C<eps2jpg>
defaults to JPG. Note that the normal installation procedure will
I<only> install C<eps2png>.

It uses GhostScript to produce the images. Since modern GhostScript
programs do not support GIF anymore, GIF images are produced via the
Portable PixMap converters (PBM-package). In this case, a temporary
file is created, named after the output file, with the extension
replaced by ".ppm". It is deleted upon completion.

=head1 ARGUMENTS

B<eps2png> always requires at least one argument: the name of the EPS
file to be converted. It is possible to specify more than one file
name. This will cause all named files to be converted into separate
files, e.g., "C<sample.eps>" will be converted to "C<sample.png>" and
so on.

=over 4

=item B<-png -pngmono -pnggray -png16 -png256 -png16m -pngalpha>

Each of these options will instruct Ghostscript to use the
corresponding bitmap generator, and supply a C<.png> default
extension for output files.

=item B<-jpg -jpggray -jpeg -jpeggray>

Same, but with a C<.jpg> default extension for output files.

=item B<-gif -gifmono>

Same, but with a C<.gif> default extension for output files.

Note: Since modern Ghostscript versions no longer support the GIF
format due to copyright restrictions, B<eps2png> will request
Ghostscript to produce a Portable Bitmap File (.ppm or .pbm) instead
and run the B<ppmtogif> converter to produce the actual GIF file.

=item B<-mono>

This option will select monochrome (BW or gray) output. It forces the
Ghostscript driver to C<pngmono>, C<jpeggray>, C<pbm>, or C<gifmono>.

=item B<-nomono>

Produces colour images. This is the default.

=item B<-width> I<NN>

The desired width of the output image.

If B<-height> is not specified, the image will be scaled proportionally.

=item B<-height> I<NN>

The desired height of the output image.

If B<-width> is not specified, the image will be scaled proportionally.

=item B<-resolution> I<NN>

Specifies the resolution for the output image. This is the width, in
pixels, of the bitmap image for an EPS image of one inch wide (72
PostScript points).

Note that for best results, use the B<-width> and B<-height> options
instead.

Default value is 82, which causes the converted image to be of more
or less the same size as the EPS image. On my screen, that is.

=item B<-scale> I<NN>

Specify a scaling factor. This may be a fractional number.

For a one-inch EPS image, the resultant bitmap image will be
I<scale> times I<resolution>.

Note that for best results, use the B<-width> and B<-height> options
instead.

=item B<-antialias> I<NN>

Sets the antialiasing depth. I<NN> must be 0 (no antialiasing), 1, 2,
4, or 8. Default value is 4.

=item B<-noantialias>

Sets the antialiasing depth to 0.

=item B<-pbm>

Forces GIF conversion through the PBM converters.

=item B<-nopbm>

Forces GIF conversion through Ghostscript.

=item B<-output> I<XXX>

Stores the output in this file. Only one input file may be supplied if
this option is specified.

=item B<-help>

Prints a help message and exits.

=item B<-ident>

Prints the program version before doing anything else.

=item B<-verbose>

Provides more verbose information.

=back

=head1 AUTHOR

Johan Vromans, <jvromans@squirrel.nl>.

=head1 BUGS

GhostScript and, if required, the PBM package, need to be installed and
accessible through the user's C<PATH>.

GhostScript is assumed to be capable of handling all the image types
listed above.

The EPS should be well-behaving.

=head1 COPYRIGHT AND DISCLAIMER

This program is Copyright 1994,2008 by Johan Vromans.
This program is free software; you can redistribute it and/or
modify it under the terms of the Perl Artistic License or the
GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

If you do not have a copy of the GNU General Public License write to
the Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
MA 02139, USA.

=cut
