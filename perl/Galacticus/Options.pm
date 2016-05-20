# Parse command line options.

package Options;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Scalar::Util 'reftype';
use XML::Simple;
require List::ExtraUtils;

sub Config {
    # Return a section of the Galacticus config file for the current host.
    my $section          = shift();
    my $xml              = new XML::Simple;
    my $configSection;
    if ( -e $galacticusPath."galacticusConfig.xml" ) {
	my $galacticusConfig = $xml->XMLin($galacticusPath."galacticusConfig.xml", KeyAttr => 0);   
	foreach ( &ExtraUtils::as_array($galacticusConfig->{$section}->{'host'}) ) {
	    if ( $_->{'name'} eq $ENV{'HOSTNAME'} || $_->{'name'} eq "default" ) {
		$configSection = $_;
		last;
	    }
	}
    }
    return $configSection;
}

sub Parse_Options {
    # Parses command line options of the form:
    #   --optionName value
    # and places them into a hash (which can have default values preloaded into it). Supports double
    # quotes in the "value" of options to allow for values that contain spaces,
    my @arguments = @{$_[0]};
    my $options   =   $_[1];
    my $iArg = -1;
    my %argumentsSeen;
    while ( $iArg < scalar(@arguments)-1 ) {
	++$iArg;
	if ( $arguments[$iArg] =~ m/^\-\-(.*)/ ) {
	    my $argument = $1;
	    my $value = $arguments[$iArg+1];
	    ++$iArg;
	    if ( $arguments[$iArg] =~ m/^\"/ ) {
		until ( $value =~ m/\"$/ ) {
		    ++$iArg;
		    $value .= $arguments[$iArg];
		}
		$value =~ s/^\"(.*)\"$/$1/;
	    }
	    if ( exists($argumentsSeen{$argument}) ) {
		if ( reftype($options->{$argument}) && reftype($options->{$argument}) eq "ARRAY" ) {
		    push(@{$options->{$argument}},$value);
		} else {
		    $options->{$argument} = [ $options->{$argument}, $value ];
		}
	    } else {
		$options->{$argument} = $value;
	    }
	    $argumentsSeen{$argument} = 1;
	}
    }
    return %argumentsSeen;
}

1;
