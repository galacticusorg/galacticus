#!/usr/bin/env perl
use File::Find;
use XML::Simple;
use Data::Dumper;

# Sort models by goodness of fit.
# Andrew Benson (2-June-2010)

# Specify directory containing models.
@modelDirectories = ( "models" );

# Specify list of parameters to include in the output.
@parametersToShow = ( "starFormationSpheroidEfficiency" );

# Specify value on which to sort - this will be treated as a regular expression
$sortName = "black hole vs\\. bulge mass relation";

# Specify whether or not to use reduced chi^2 for sort.
$useReducedChi2 = 0;

# Scan for model directories.
find(\&processFile, @modelDirectories);

# Output header.
print "# Model\tFit";
foreach $parameter ( @parametersToShow ) {
    print "\t".$parameter;
}
print "\n";

# Output sorted list of models.
foreach $value (sort {$results{$a}->{'result'} <=> $results{$b}->{'result'} } keys %results ) {
    print $value."\t".$results{$value}->{'result'};
    foreach $parameter ( @parametersToShow ) {
	print "\t".$results{$value}->{'parameters'}->{$parameter};
    }
    print "\n";
}

exit;

sub processFile {

    $fileName = $_;
    chomp($fileName);

    # Extract fit value.
    if ( $fileName =~ m/^galacticusFits\.xml(\.bz2)??/ ) {
	$fitFile = "galacticusFits.xml";
	if ( -e $fitFile.".bz2" ) {
	    system("bunzip2 ".$fitFile.".bz2");
	    $reCompressFits = 1;
	} else {
	    $reCompressFits = 0;
	}
	if ( -e $fitFile ) {
	    $xml = new XML::Simple;
	    $data = $xml->XMLin($fitFile,ForceArray => 1);
	    foreach $fit ( @{$data->{'galacticusFit'}} ) {
		if ( ${$fit->{'name'}}[0] =~ m/$sortName/ ) {
		    $modelName = $File::Find::dir;
		    unless ( ${$fit->{'chiSquared'}}[0] eq "nan" ) {
			$results{$modelName}->{'result'} = ${$fit->{'chiSquared'}}[0];
			if ( $useReducedChi2== 1 ) {
			    $results{$modelName}->{'result'} /= ${$fit->{'degreesOfFreedom'}}[0];
			}
		    }
		}
	    }
	    system("bzip2 ".$fitFile) if ( $reCompressFits == 1 );
	}
    }

    # Extract any parameters.
    if ( $fileName =~ m/^newParameters\.xml(\.bz2)??/ ) {
	# Include parameters that must be shown.
	$modelName = $File::Find::dir;	
	if ( $#parametersToShow >= 0 ) {
	    if ( -e "newParameters.xml.bz2" ) {
		system("bunzip2 newParameters.xml.bz2");
		$reCompressParameters = 1;
	    } else {
		$reCompressParameters = 0;
	    }
	    $xml = new XML::Simple;
	    $data = $xml->XMLin("newParameters.xml");
	    foreach $parameter ( @parametersToShow ) {
		$results{$modelName}->{'parameters'}->{$parameter} = $data->{'parameter'}->{$parameter}->{'value'};
	    }
	    if ( $results{$modelName}->{'result'} eq "" ) {$results{$modelName}->{'result'} = "undefined"};
	}
    }
}

