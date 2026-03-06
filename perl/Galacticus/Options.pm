# Parse command line options.

package Galacticus::Options;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Scalar::Util 'reftype';
use XML::Simple;
use List::ExtraUtils;

sub LoadConfig {
    # Load and return a Galacticus config file.
    my $config;
    if ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/galacticusConfig.xml" ) {
	my $xml = new XML::Simple;
	$config = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/galacticusConfig.xml" );
    } elsif ( -e $ENV{'HOME'}."/.galacticusConfig.xml" ) {
	my $xml = new XML::Simple;
	$config = $xml->XMLin($ENV{'HOME'                }."/.galacticusConfig.xml");
    }
}

sub Config {
    # Return a section of the Galacticus config file for the current host.
    my $section          = shift();
    my $xml              = new XML::Simple;
    my $configSection;
    if ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/galacticusConfig.xml" ) {
	my $galacticusConfig = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/galacticusConfig.xml", KeyAttr => 0);   
	foreach ( &List::ExtraUtils::as_array($galacticusConfig->{$section}->{'host'}) ) {
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
	die("unexpected dash in argument '".$arguments[$iArg]."'? can happen as a result of copy-and-paste - try changing to regular dashes")
	    if ( $arguments[$iArg] =~ m/^\-?[––—‒‑‐]/ );
	if ( $arguments[$iArg] =~ m/^\-\-(.*)/ ) {
	    my $argument = $1;
	    die("Galacticus::Options::Parse_Options: missing final value after option '--".$argument."'")
		unless ( $iArg+1 < scalar(@arguments) );
	    die("Galacticus::Options::Parse_Options: missing value after option '--".$argument."'")
		if ( $arguments[$iArg+1] =~ m/^\-\-/ );
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

sub Serialize_Options {
    # Serialize options back to command line arguments.
    my %options = %{shift()};
    my $commandLine;
    foreach my $name ( sort(keys(%options)) ) {
	$commandLine .= " ".join(" ",map {"--".$name." ".$_} &List::ExtraUtils::as_array($options{$name}));
    }
    return $commandLine;
}

1;
