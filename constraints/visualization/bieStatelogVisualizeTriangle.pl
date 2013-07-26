#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Data::Dumper;

# Visualize the posterior distribution from a BIE statelog file using a triangle arrangement.
# Andrew Benson (12-June-2012)

# Get file name to process.
die("Usage: bieStatelogVisualizeTriangle.pl filename [options]") unless ( scalar(@ARGV) > 0 );
my $fileName = $ARGV[0];

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
		    if ( $propertySpecifiers[$i] eq "log" );
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

# Construct options to pass to plotting script.
my $ngood = 1000;
$ngood = $arguments{'ngood'}
    if ( exists( $arguments{'ngood'}) );
my $ngrid = 50;
$ngrid = $arguments{'ngrid'}
    if ( exists( $arguments{'ngrid'}) );
my $options = " --title none --ngood ".$ngood." --ngrid ".$ngrid;

# Open the file and parse the available parameter names.
open(iHndl,$fileName);
my $header = <iHndl>;
close(iHndl);
$header =~ s/^.*\"Prior\"\s*\"//;
$header =~ s/\"\s*$//;
my @propertiesAvailable = split(/\"\s*\"/,$header);
my @properties;
if ( exists($arguments{'property'}) ) {
    @properties = @{$arguments{'property'}};
    die("bieStatelogVisualizeTriangle.pl: at least 3 properties must be specified")
	if ( scalar(@properties) < 3 );
    foreach my $property ( @properties ) {
	die("Property '".$property->{'name'}."' is not available") unless ( grep {$_ eq $property->{'name'}} @propertiesAvailable );
    }
} else {
    @properties = @propertiesAvailable;
}

# Determine any ranges to be applied to parameters.
my @ranges;
if ( exists($arguments{'range'}) ) {
    @ranges = @{$arguments{'range'}};
    foreach my $property ( @ranges ) {
	my @range = split(/:/,$property);
	die("Property '".$range[0]."' is not available") unless ( grep {$_ eq $range[0]} @propertiesAvailable );
    }
}

# Loop over parameters.
my $standardWidth;
my $standardHeight;
for(my $i=0;$i<scalar(@properties);++$i) {
    my $command = "constraints/visualization/bieStatelogVisualize.pl ".$fileName." --xProperty '".$properties[$i]->{'name'}."' --xScale ".$properties[$i]->{'scaling'}." --output ".$outputFileName."_".$i.".pdf".$options;
    $command .= " --xLabel '".$properties[$i]->{'xLabel'}."'"
        if ( exists($properties[$i]->{'xLabel'}) );
    $command .= " --zLabel '".$properties[$i]->{'zLabel'}."'"
        if ( exists($properties[$i]->{'zLabel'}) );
    if ( @ranges ) {
	foreach my $range ( @ranges ) {
	    $command .= " --range ".$range;
	}
    }
    system($command)
	unless ( -e $outputFileName."_".$i.".pdf" );
    if ( $i < scalar(@properties)-1 ) { 
	for(my $j=$i+1;$j<scalar(@properties);++$j) {
	    my $command = "constraints/visualization/bieStatelogVisualize.pl ".$fileName." --yProperty '".$properties[$i]->{'name'}."' --yScale ".$properties[$i]->{'scaling'}." --xProperty '".$properties[$j]->{'name'}."' --xScale ".$properties[$j]->{'scaling'}." --output ".$outputFileName."_".$i."_".$j.".pdf".$options;
	    $command .= " --xLabel '".$properties[$j]->{'xLabel'}."'"
		if ( exists($properties[$j]->{'xLabel'}) );
	    $command .= " --yLabel '".$properties[$i]->{'xLabel'}."'"
		if ( exists($properties[$i]->{'xLabel'}) );
	    if ( $j == scalar(@properties)-1 ) {
		$command .= " --labels y2 --colorbox 0";
	    } else {
		$command .= " --labels none --colorbox 0";
	    }
	    system($command)
		unless ( -e $outputFileName."_".$i."_".$j.".pdf" );
	    if ( $i == 0 && $j == 1 ) {
		open(iHndl,"pdfinfo ".$outputFileName."_".$i."_".$j.".pdf|");
		while ( my $line = <iHndl> ) {
		    if ( $line =~ m/Page size:\s*(\d+)\s*x\s*(\d+)\s*pts/ ) {
			$standardWidth  = $1;
			$standardHeight = $2;
		    }		    
		}
		close(iHndl);
	    }
	}
    }
}

# Open output file.
open(oHndl,">".$outputFileName.".tex");
print oHndl "\\renewcommand{\\arraystretch}{0}\n";
print oHndl "\\begin{tabular}{".("l\@{}" x scalar(@properties))."}\n";

# Loop over parameters.
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
    print oHndl "\\includegraphics[scale=".$scale."]{".$outputFileName."_".$i.".pdf}";
    print oHndl "\\vspace{-".$shiftVertical."pt}";
    if ( $i < scalar(@properties)-1 ) { 
	for(my $j=$i+1;$j<scalar(@properties);++$j) {
	    print oHndl "&\\raisebox{".$shiftVertical."pt}{\\includegraphics[scale=".$scale."]{".$outputFileName."_".$i."_".$j.".pdf}}";
	}
    }
    print oHndl "\\\\\n";
}

# Close the output file.
print oHndl "\\end{tabular}\n";
close(oHndl);

exit;

