#!/usr/bin/env perl
use File::Find;
use XML::Simple;
use Data::Dumper;
use Date::Parse;
use Text::Table;

# Sort models by goodness of fit.
# Andrew Benson (2-June-2010)

# Set default list of model directories to scan.
$parameters{'modelDirectories'} = "projects/parameterSearch/models/model0013";

# Set default list of parameters to show.
$parameters{'parametersToShow' } = "starFormationDiskEfficiency";

# Set default fit value on which to sort.
$parameters{'sortOn'} = [
			 {
			     name => "Volume averaged star formation rate history.",
			     weight => 10.0
			     },
			 {
			     name => "Li and White (2009) stellar mass function",
			     weight => 0.25
			     },
			 {
			     name => "Feoli & Mancini (2009) black hole vs. bulge mass relation",
			     weight =>0.001
			     },
			 {
			     name => "Pizagno et al. (2007) SDSS Tully-Fisher relation",
			     weight => 0.001
			     },
			 {
			     name => "Dejong & Lacey (2000) disk scale length distributions",
			     weight => 0.075
			     },
			 {
			     name => "Weinmann et al. (2006) SDSS galaxy colors",
			     weight => 0.033
			     },
			 ];

# Set default sort methods.
$parameters{'useReducedChi2'} = 1;
$parameters{'sortOnRating'}   = 1;

# Parse command line parameters.
for($iParameter=0;$iParameter<=$#ARGV;++$iParameter) {
    $parameter = $ARGV[$iParameter];
    if ( $parameter =~ m/\"/ ) {
	do {
	    ++$iParameter;
	    $parameter .= $ARGV[$iParameter];
	} until ( $parameter =~ m/\"$/ );
    }
    if ( $parameter =~ m/(\S+)=(.*)/ ) {
	$parameterName  = $1;
	$parameterValue = $2;
	$parameters{$parameterName} = $parameterValue;
    }
}

# Specify directory containing models.
@modelDirectories = split(/,/,$parameters{'modelDirectories'});

# Specify list of parameters to include in the output.
@parametersToShow = split(/,/,$parameters{'parametersToShow'});

# Specify minimum date/time to include. (Applies only to archived models.)
if ( exists($parameters{'modelsRunAfter'}) ) {
    $timeMinimum = str2time($parameters{'modelsRunAfter'});
    print "Sort includes models run after: ".$parameters{'modelsRunAfter'}."\n";
} else {
    $timeMinimum = 0;
}

# Specify value on which to sort - this will be treated as a regular expression
$sortOn = $parameters{'sortOn'};
print "Sorting on fit:\n";
print Dumper($parameters{'sortOn'});

# Specify whether or not to use reduced chi^2 for sort.
$useReducedChi2 = $parameters{'useReducedChi2'};
if ( $useReducedChi2 ) {
    print "Sorting using: reduced chi^2\n";
} else {
    print "Sorting using: chi^2\n";
}

# Specify whether or not to sort on rating.
$sortOnRating = $parameters{'sortOnRating'};
if ( $sortOnRating ) {
    print "Sorting using: ratings then chi^2\n";
} else {
    print "Sorting using: chi^2\n";
}

# Scan for model directories.
find(\&processFile, @modelDirectories);

# Create header.
$table = Text::Table->new();
undef(@rowData);
$rowData[0] = "Model";
$rowData[1] = "Rating";
$rowData[2] = "Chi^2";
push(@rowData,@parametersToShow);
$table->add(@rowData);

# Output sorted list of models.
sub sortOrder {
    $result = 0;
    if ( $sortOnRating == 1 ) {
	if ( $results{$a}->{'rating'} eq "undefined" ) {
	    if ( $results{$b}->{'rating'} eq "undefined" ) {
		$result = 0;
	    } else {
		$result = 1;
	    }
	} elsif ( $results{$b}->{'rating'} eq "undefined" ) {
	    if ( $results{$a}->{'rating'} eq "undefined" ) {
		$result = 0;
	    } else {
		$result = -1;
	    }
	} else {
	    $result = $results{$b}->{'rating'} <=> $results{$a}->{'rating'};
	}
    }
    if ( $result == 0 ) {
	if ( $results{$a}->{'chiSquared'} eq "undefined" ) {
	    if ( $results{$b}->{'chiSquared'} eq "undefined" ) {
		$result = 0;
	    } else {
		$result = 1;
	    }
	} elsif ( $results{$b}->{'chiSquared'} eq "undefined" ) {
	    if ( $results{$a}->{'chiSquared'} eq "undefined" ) {
		$result = 0;
	    } else {
		$result = -1;
	    }
	} else {
	    $result = $results{$a}->{'chiSquared'} <=> $results{$b}->{'chiSquared'};
	}
    }
    return $result;
}
foreach $value ( sort sortOrder keys %results ) {
    undef(@rowData);
    $rowData[0] = $value;
    $rowData[1] = $results{$value}->{'rating'};
    $rowData[2] = $results{$value}->{'chiSquared'};
    foreach $parameter ( @parametersToShow ) {
	push(@rowData,$results{$value}->{'parameters'}->{$parameter});
    }
    $table->add(@rowData);
}
print "\n";
print $table;

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
		foreach my $sortOn ( @{$parameters{'sortOn'}} ) {
		    if ( ${$fit->{'name'}}[0] eq $sortOn->{'name'} ) {
		       $modelName = $File::Find::dir;
		       unless ( ${$fit->{'chiSquared'}}[0] eq "nan" ) {
		          my $chiSquared = ${$fit->{'chiSquared'}}[0];
			  if ( $useReducedChi2== 1 ) {
			      $chiSquared /= ${$fit->{'degreesOfFreedom'}}[0];
	                  }
	                  $results{$modelName}->{'chiSquared'} += $chiSquared*$sortOn->{'weight'};
	               }
	               $results{$modelName}->{'rating'} = ${${$fit->{'rating'}}[0]->{'value'}}[0];
		    }
		}
	    }
	    system("bzip2 ".$fitFile) if ( $reCompressFits == 1 );
	}
	if ( $results{$modelName}->{'chiSquared'} eq "" ) {$results{$modelName}->{'chiSquared'} = "undefined"};
	if ( $results{$modelName}->{'rating'}     eq "" ) {$results{$modelName}->{'rating'}     = "undefined"};
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
	}
	if ( $results{$modelName}->{'chiSquared'} eq "" ) {$results{$modelName}->{'chiSquared'} = "undefined"};
	if ( $results{$modelName}->{'rating'}     eq "" ) {$results{$modelName}->{'rating'}     = "undefined"};
    }

    # Extract fit value.
    if ( $fileName =~ m/^[\d\.]+\.xml(\.bz2)??/ ) {
	($fitFile = $fileName) =~ s/\.bz2$//;
	if ( -e $fitFile.".bz2" ) {
	    system("bunzip2 ".$fitFile.".bz2");
	    $reCompressFits = 1;
	} else {
	    $reCompressFits = 0;
	}
	if ( -e $fitFile ) {
	    $xml = new XML::Simple;
	    $dataFit = $xml->XMLin($fitFile,ForceArray => 1);
	    $xml = new XML::Simple;
	    $data    = $xml->XMLin($fitFile);
            # Check the date.
	    $time = str2time($data->{'version'}->{'runTime'});
	    if ( $time >= $timeMinimum ) {
		$modelName = $fitFile;
		# Extract fit results.
		foreach $fit ( @{${$dataFit->{'galacticusFits'}}[0]->{'galacticusFit'}} ) {
		    if ( ${$fit->{'name'}}[0] =~ m/$sortOn/ ) {
			unless ( ${$fit->{'chiSquared'}}[0] eq "nan" ) {
			    $results{$modelName}->{'chiSquared'} = ${$fit->{'chiSquared'}}[0];
			    if ( $useReducedChi2== 1 ) {
				$results{$modelName}->{'chiSquared'} /= ${$fit->{'degreesOfFreedom'}}[0];
			    }
			}
			$results{$modelName}->{'rating'} = ${${$fit->{'rating'}}[0]->{'value'}}[0];
		    }
		}
		# Extract any parameters.
		if ( $#parametersToShow >= 0 ) {
		    foreach $parameter ( @parametersToShow ) {
			$results{$modelName}->{'parameters'}->{$parameter} =$data->{'parameters'}->{'parameter'}->{$parameter}->{'value'};
		    }
		}
		system("bzip2 ".$fitFile) if ( $reCompressFits == 1 );
		if ( $results{$modelName}->{'chiSquared'} eq "" ) {$results{$modelName}->{'chiSquared'} = "undefined"};
		if ( $results{$modelName}->{'rating'}     eq "" ) {$results{$modelName}->{'rating'}     = "undefined"};
	    }
	}
    }
}

