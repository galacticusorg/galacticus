#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Data::Dumper;
use XML::Simple;
require Galacticus::Launch::PBS;

# Visualize the posterior distribution from MCMC chains using a triangle arrangement.
# Andrew Benson (12-June-2012)

# Get file name to process.
die("Usage: mcmcVisualizeTriangle.pl fileRoot configFile [options]")
    unless ( scalar(@ARGV) > 1 );
my $fileRoot       = $ARGV[0];
my $configFileName = $ARGV[1];

# Create a hash of named arguments.
my %arguments;
my $iArg = 0;
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	my $argument = $1;
	my $value = $ARGV[$iArg+1];
	++$iArg;
	if ( $ARGV[$iArg] =~ m/^\"/ ) {
	    until ( $value =~ m/\"$/ ) {
		++$iArg;
		$value .= $ARGV[$iArg];
	    }
	   $value =~ s/^\"(.*)\"$/$1/;
	}
	if ( $argument eq "property" ) {
	    my @propertySpecifiers = split(/:/,$value);
	    my %propertyData = 
		(
		 name    => $propertySpecifiers[0],
		 scaling => "linear"
		);
	    for(my $i=1;$i<scalar(@propertySpecifiers);++$i) {
		$propertyData{'scaling'} = "log"
		    if ( $propertySpecifiers[$i] eq "logarithmic" );
		if ( $propertySpecifiers[$i] =~ m/^xLabel=(.*)/ ) {
		    $propertyData{'xLabel'} = $1;
		}
		if ( $propertySpecifiers[$i] =~ m/^zLabel=(.*)/ ) {
		    $propertyData{'zLabel'} = $1;
		}
	    }
	    push(@{$arguments{$argument}},\%propertyData);
	} elsif ( $argument eq "range" ) {
	    push(@{$arguments{$argument}},$value);
	} else {
	    $arguments{$argument} = $value;
	}
    }
}

# Set output file name.
my $outputFileName = "triangle";
$outputFileName = $arguments{'output'}
    if ( exists($arguments{'output'}) );

# Set scale for graphics.
my $scale = 0.4;
$scale = $arguments{'scale'}
    if ( exists($arguments{'scale'}) );
my $textSize = 7;
$textSize = $arguments{'textSize'}
    if ( exists($arguments{'textSize'}) );
my $labelStyle = "normalsize";
$labelStyle = $arguments{'labelStyle'}
    if ( exists($arguments{'labelStyle'}) );

# Set work directory.
my $workDirectory = ".";
$workDirectory = $arguments{'workDirectory'}
    if ( exists($arguments{'workDirectory'}) );

# Construct options to pass to plotting script.
my $ngood = 1000;
$ngood = $arguments{'ngood'}
    if ( exists( $arguments{'ngood'}) );
my $ngrid = 50;
$ngrid = $arguments{'ngrid'}
    if ( exists( $arguments{'ngrid'}) );
my $options = " --title none --ngood ".$ngood." --ngrid ".$ngrid;
$options .= " --useUnconverged ".$arguments{'useUnconverged'}
   if ( exists($arguments{'useUnconverged'}) );

# Open the config file and parse the available parameter names.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFileName, KeyAttr => []);
my @propertyNamesAvailable;
foreach my $parameter ( @{$config->{'parameters'}->{'parameter'}} ) {
    push(@propertyNamesAvailable,$parameter->{'name'})
	 if ( exists($parameter->{'prior'}) );
}
my @propertyNames;
if ( exists($arguments{'property'}) ) {
    @propertyNames = @{$arguments{'property'}};
    die("mcmcVisualizeTriangle.pl: at least 3 propertyNames must be specified")
	if ( scalar(@propertyNames) < 3 );
    foreach my $property ( @propertyNames ) {
	die("Property '".$property->{'name'}."' is not available") unless ( grep {$_ eq $property->{'name'}} @propertyNamesAvailable );
    }
} else {
    @propertyNames = @propertyNamesAvailable;
}

# Determine any ranges to be applied to parameters.
my @ranges;
if ( exists($arguments{'range'}) ) {
    @ranges = @{$arguments{'range'}};
    foreach my $property ( @ranges ) {
	my @range = split(/:/,$property);
	die("Property '".$range[0]."' is not available") unless ( grep {$_ eq $range[0]} @propertyNamesAvailable );
    }
}

# Construct property list.
my @properties;
foreach my $parameter ( @{$config->{'parameters'}->{'parameter'}} ) {
    print $parameter->{'name'}."\t".join(" ",@propertyNames)."\n";
    push(@properties,$parameter)
	if ( grep {$parameter->{'name'} eq $_} @propertyNames );
}

# Loop over parameters.
my @pbsStack;
my $standardWidth;
my $standardHeight;
for(my $i=0;$i<scalar(@properties);++$i) {
    my $command = "constraints/visualization/mcmcVisualize.pl ".$configFileName." ".$fileRoot." --workDirectory ".$workDirectory." --xProperty '".$properties[$i]->{'name'}."' --xScale ".$properties[$i]->{'mapping'}->{'type'}." --textSize ".$textSize." --labelStyle ".$labelStyle." --output ".$outputFileName."_".$i.".pdf --data ".$outputFileName."_".$i.".xml ".$options;
    $command .= " --xLabel '".$properties[$i]->{'xLabel'}."'"
        if ( exists($properties[$i]->{'xLabel'}) );
    $command .= " --zLabel '".$properties[$i]->{'zLabel'}."'"
        if ( exists($properties[$i]->{'zLabel'}) );
    if ( @ranges ) {
	foreach my $range ( @ranges ) {
	    $command .= " --range ".$range;
	}
    }
    # Create PBS job.
    unless ( -e $outputFileName."_".$i.".pdf" ) {	
	my %job =
	    (
	     launchFile => $outputFileName."_".$i.".pbs",
	     label      => "triangle",
	     logFile    => $outputFileName."_".$i.".log",
	     command    => $command,
	     ppn        => 1
	    );
	foreach ( 'walltime', 'memory' ) {
	    $job{$_} = $arguments{$_}
	    if ( exists($arguments{$_}) );
	}
	push(@pbsStack,\%job);
    }
    if ( $i < scalar(@properties)-1 ) { 
	for(my $j=$i+1;$j<scalar(@properties);++$j) {
	    my $command = "constraints/visualization/mcmcVisualize.pl ".$configFileName." ".$fileRoot." --workDirectory ".$workDirectory." --yProperty '".$properties[$i]->{'name'}."' --yScale ".$properties[$i]->{'mapping'}->{'type'}." --xProperty '".$properties[$j]->{'name'}."' --xScale ".$properties[$j]->{'mapping'}->{'type'}." --textSize ".$textSize." --labelStyle ".$labelStyle." --output ".$outputFileName."_".$i."_".$j.".pdf --data ".$outputFileName."_".$i."_".$j.".xml ".$options;
	    $command .= " --xLabel '".$properties[$j]->{'xLabel'}."'"
		if ( exists($properties[$j]->{'xLabel'}) );
	    $command .= " --yLabel '".$properties[$i]->{'xLabel'}."'"
		if ( exists($properties[$i]->{'xLabel'}) );
	    if ( $j == scalar(@properties)-1 ) {
		$command .= " --labels y2 --colorbox 0";
	    } else {
		$command .= " --labels none --colorbox 0";
	    }
	    # Create PBS job.
	    unless ( -e $outputFileName."_".$i."_".$j.".pdf" ) {
	    	my %job =
		    (
		     launchFile => $outputFileName."_".$i."_".$j.".pbs",
		     label      => "triangle",
		     logFile    => $outputFileName."_".$i."_".$j.".log",
		     command    => $command,
		     ppn        => 1
		    );
		foreach ( 'walltime', 'memory' ) {
		    $job{$_} = $arguments{$_}
		    if ( exists($arguments{$_}) );
		}
		push(@pbsStack,\%job);
	    }
	}
    }
}
# Send jobs to PBS.
print Dumper(\@pbsStack);
&PBS::SubmitJobs(\%arguments,@pbsStack);
exit;

# Get output size.
open(iHndl,"pdfinfo ".$outputFileName."_0_1.pdf|");
while ( my $line = <iHndl> ) {
    if ( $line =~ m/Page size:\s*(\d+)\s*x\s*(\d+)\s*pts/ ) {
	$standardWidth  = $1;
	$standardHeight = $2;
    }		    
}
close(iHndl);

# Open output file.
open(oHndl,">".$outputFileName.".tex");

# Loop over parameters.
my $outputDirectoryName = `dirname  $outputFileName`;
my $outputBaseName      = `basename $outputFileName`;
chomp($outputDirectoryName);
chomp($outputBaseName     );
print oHndl "\\newcommand{\\triangledir}{".$outputDirectoryName."}\n";
print oHndl "\\renewcommand{\\arraystretch}{0}\n";
print oHndl "\\begin{tabular}{".("l\@{}" x scalar(@properties))."}\n";
for(my $i=0;$i<scalar(@properties);++$i) {
    print oHndl ("&" x $i);
    my $width;
    my $height;
    open(iHndl,"pdfinfo ".$outputFileName."_".$i.".pdf|");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/Page size:\s*(\d+)\s*x\s*(\d+)\s*pts/ ) {
	    $width  = $1;
	    $height = $2;
	}		    
    }
    close(iHndl);
    my $shiftHorizontal = $scale*($width -$standardWidth );
    my $shiftVertical   = $scale*($height-$standardHeight);
    print oHndl "\\hspace{-".$shiftHorizontal."pt}"
	if ( $i > 0 );
    print oHndl "\\includegraphics[scale=".$scale."]{\\triangledir/".$outputBaseName."_".$i.".pdf}";
    print oHndl "\\vspace{-".$shiftVertical."pt}";
    if ( $i < scalar(@properties)-1 ) { 
	for(my $j=$i+1;$j<scalar(@properties);++$j) {
	    print oHndl "&\\raisebox{".$shiftVertical."pt}{\\includegraphics[scale=".$scale."]{\\triangledir/".$outputBaseName."_".$i."_".$j.".pdf}}";
	}
    }
    print oHndl "\\\\\n";
}

# Close the output file.
print oHndl "\\end{tabular}\n";
close(oHndl);

exit;

