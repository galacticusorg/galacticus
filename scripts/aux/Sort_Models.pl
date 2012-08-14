#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use XML::Simple;
use Data::Dumper;
use Date::Parse;
use Text::Table;

# Sort models by goodness of fit.
# Andrew Benson (2-June-2010)

# Initialize data structures.
our %parameters;
my  %results;

# Set default list of model directories to scan.
$parameters{'modelDirectories'} = ".";

# Set default list of parameters to show.
$parameters{'parametersToShow' } = "";

# Get the name of the file containing the list of fits to sort on.
die("Usage: Sort_Models.pl <fits> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $fitsFile = $ARGV[0];

# Set default fit value on which to sort.
my $fitsXML = new XML::Simple;
$parameters{'sortOn'} = $fitsXML->XMLin($fitsFile, KeyAttr => "", ForceArray => [ "analysis" ]);

# Set default sort methods.
$parameters{'useReducedChi2'} = 1;
$parameters{'sortOnRating'}   = 1;

# Set default verbosity.
$parameters{'quiet'} = 0;

# Parse command line parameters.
for(my $iParameter=1;$iParameter<scalar(@ARGV);++$iParameter) {
    my $parameter = $ARGV[$iParameter];
    if ( $parameter =~ m/\"/ ) {
	do {
	    ++$iParameter;
	    $parameter .= $ARGV[$iParameter];
	} until ( $parameter =~ m/\"$/ );
    }
    if ( $parameter =~ m/\-\-(\S+)=(.*)/ ) {
	my $parameterName  = $1;
	my $parameterValue = $2;
	$parameters{$parameterName} = $parameterValue;
    }
}

# Specify directory containing models.
my @modelDirectories = split(/,/,$parameters{'modelDirectories'});

# Specify list of parameters to include in the output.
my @parametersToShow = split(/,/,$parameters{'parametersToShow'});

# Specify minimum date/time to include. (Applies only to archived models.)
my $timeMinimum = 0;
if ( exists($parameters{'modelsRunAfter'}) ) {
    $timeMinimum = str2time($parameters{'modelsRunAfter'});
    print "Sort includes models run after: ".$parameters{'modelsRunAfter'}."\n"
	unless ( $parameters{'quiet'} == 1 );
}

# Specify value on which to sort - this will be treated as a regular expression
unless ( $parameters{'quiet'} == 1 ) {
    print "Sorting on fit:\n";
    my $table = Text::Table->new();
    foreach my $fit ( @{$parameters{'sortOn'}->{'analysis'}} ) {
	my @rowData = ( $fit->{'name'}, "[weight: ".$fit->{'weight'}."]" );
	$table->add(@rowData);
    }
    print $table;
}

# Specify whether or not to use reduced chi^2 for sort.
unless ( $parameters{'quiet'} == 1 ) {
    if ( $parameters{'useReducedChi2'} ) {
	print "Sorting using: reduced chi^2\n";
    } else {
	print "Sorting using: chi^2\n";
    }
}

# Specify whether or not to sort on rating.
unless ( $parameters{'quiet'} == 1 ) {
    if ( $parameters{'sortOnRating'} ) {
	print "Sorting using: ratings then chi^2\n";
    } else {
	print "Sorting using: chi^2\n";
    }
}

# Scan for model directories.
find(\&processFile, @modelDirectories);

# Create header.
my $table = Text::Table->new();
my @rowData;
undef(@rowData);
$rowData[0] = "Model";
$rowData[1] = "Rating";
$rowData[2] = "Chi^2";
push(@rowData,@parametersToShow);
$table->add(@rowData);

# Output sorted list of models.
sub sortOrder {
    my $result = 0;
    if ( $parameters{'sortOnRating'} == 1 ) {
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
foreach my $value ( sort sortOrder keys %results ) {
    undef(@rowData);
    $rowData[0] = $value;
    $rowData[1] = $results{$value}->{'rating'};
    $rowData[2] = $results{$value}->{'chiSquared'};
    foreach my $parameter ( @parametersToShow ) {
	push(@rowData,$results{$value}->{'parameters'}->{$parameter});
    }
    $table->add(@rowData);
}
print "\n"
    unless ( $parameters{'quiet'} == 1 );
print $table;

exit;

sub processFile {

    my $fileName = $_;
    chomp($fileName);

     # Extract fit value.
     if ( $fileName =~ m/^galacticusFits\.xml(\.bz2)??/ ) {
	 # Find the model directory.
	 my $modelName = $File::Find::dir;	
	 print "Processing: ".$modelName."\n"
	     unless ( $parameters{'quiet'} == 1 );
	 my $fitFile = "galacticusFits.xml";
	 system("bunzip2 ".$fitFile.".bz2")
	     if ( -e $fitFile.".bz2" );
	 if ( -e $fitFile ) {
	     my $xml = new XML::Simple;
	     my $data = $xml->XMLin($fitFile, KeyAttr => "", ForceArray => [ "galacticusFit" ]);
	     foreach my $fit ( @{$data->{'galacticusFit'}} ) {
		 foreach my $sortOn ( @{$parameters{'sortOn'}->{'analysis'}} ) {
		     if ( $fit->{'name'} eq $sortOn->{'name'} ) {
			 unless ( $fit->{'chiSquared'} eq "nan" ) {
			     my $chiSquared = $fit->{'chiSquared'};
			     $chiSquared /= $fit->{'degreesOfFreedom'}
			        if ( $parameters{'useReducedChi2'} == 1 );
			     $results{$modelName}->{'chiSquared'} += $chiSquared*$sortOn->{'weight'};
			 }
			 $results{$modelName}->{'rating'} = $fit->{'rating'}->{'value'};
		     }
		 }
	     }
	 }
	 $results{$modelName}->{'chiSquared'} = "undefined"
	     unless ( defined($results{$modelName}->{'chiSquared'}) );
	 $results{$modelName}->{'rating'}     = "undefined"
	     unless ( defined($results{$modelName}->{'rating'    }) );
	 
	 # Include parameters that must be shown.
	 if ( scalar(@parametersToShow) > 0 ) {
	     system("bunzip2 newParameters.xml.bz2")
		 if ( -e "newParameters.xml.bz2" );
	     my $xml = new XML::Simple;
	     my $data = $xml->XMLin("newParameters.xml");
	     foreach my $parameter ( @parametersToShow ) {
		 $results{$modelName}->{'parameters'}->{$parameter} = $data->{'parameter'}->{$parameter}->{'value'};
	     }
	 }
     }
}

