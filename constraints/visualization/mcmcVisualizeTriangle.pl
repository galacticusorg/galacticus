#!/usr/bin/env perl
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}.'/perl';
use strict;
use warnings;
use Data::Dumper;
use XML::Simple;
use POSIX;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::Local;
use Galacticus::Launch::PBS;
use Galacticus::Launch::MonolithicPBS;
use Galacticus::Launch::Slurm;

# Visualize the posterior distribution from MCMC chains using a triangle arrangement.
# Andrew Benson (12-June-2012)

# Get file name to process.
die("Usage: mcmcVisualizeTriangle.pl fileRoot parameterFile [options]")
    unless ( scalar(@ARGV) > 1 );
my $fileRoot          = $ARGV[0];
my $parameterFileName = $ARGV[1];

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
	    my @propertySpecifiers = split(/(?<!:):(?!:)/,$value);
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
# Set number of significant figures in labels.
my $sigFigs = 2;
$sigFigs = $arguments{'sigFigs'}
    if ( exists($arguments{'sigFigs'}) );
# Set scale for graphics.
my $scale = 0.4;
$scale = $arguments{'scale'}
    if ( exists($arguments{'scale'}) );
my $textSize = 7;
$textSize = $arguments{'textSize'}
    if ( exists($arguments{'textSize'}) );
my $plotSize = "5cm,5cm";
$plotSize = $arguments{'plotSize'}
    if ( exists($arguments{'plotSize'}) );
my $labelStyle = "normalsize";
$labelStyle = $arguments{'labelStyle'}
    if ( exists($arguments{'labelStyle'}) );
my $drawLabels = "gnuplot";
$drawLabels = $arguments{'drawLabels'}
    if ( exists($arguments{'drawLabels'}) );
my $lineWeight = "2,1";
$lineWeight = $arguments{'lineWeight'}
    if ( exists($arguments{'lineWeight'}) );
my $shortcut = ".";
$shortcut = $arguments{'shortcut'}
    if ( exists($arguments{'shortcut'}) );

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

# Open the parameter file and parse the available parameter names.
my $xml        = new XML::Simple;
my $parameters = $xml->XMLin($parameterFileName);
my @propertyNamesAvailable;
foreach my $parameter ( @{$parameters->{'posteriorSampleSimulation'}->{'modelParameter'}} ) {
    push(@propertyNamesAvailable,$parameter->{'name'}->{'value'})
	 if ( $parameter->{'value'} eq "active" );
}
my @propertyNames;
if ( exists($arguments{'property'}) ) {
    @propertyNames = map {$_->{'name'}} @{$arguments{'property'}};
    die("mcmcVisualizeTriangle.pl: at least 2 propertyNames must be specified")
	if ( scalar(@propertyNames) < 2 );
    foreach my $property ( @propertyNames ) {
	die("Property '".$property."' is not available") unless ( grep {$_ eq $property} @propertyNamesAvailable );
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
foreach my $parameter ( @{$parameters->{'posteriorSampleSimulation'}->{'modelParameter'}} ) {
    if ( grep {$parameter->{'name'}->{'value'} eq $_} @propertyNames ) {
	foreach my $property ( @{$arguments{'property'}} ) {
	    if ( $property->{'name'} eq $parameter->{'name'}->{'value'} ) {
		foreach ( 'xLabel', 'zLabel' ) {
		    $parameter->{$_} = $property->{$_};
		}
	    }
	}
	push(@properties,$parameter);
    }
}

# Loop over parameters.
my $queueManager = &Galacticus::Options::Config('queueManager');
my @jobStack;
my $standardWidth;
my $standardHeight;
for(my $i=0;$i<scalar(@properties);++$i) {
    my $command = "constraints/visualization/mcmcVisualize.pl ".$parameterFileName." ".$fileRoot." --workDirectory ".$workDirectory." --xProperty '".$properties[$i]->{'name'}->{'value'}."' --xScale ".$properties[$i]->{'operatorUnaryMapper'}->{'value'}." --textSize ".$textSize." --plotSize ".$plotSize." --lineWeight ".$lineWeight." --labelStyle ".$labelStyle." --output ".$outputFileName."_".$i.".pdf --data ".$outputFileName."_".$i.".xml ".$options;
    $command .= " --oldChainFormat ".$arguments{'oldChainFormat'}
	if ( exists($arguments{'oldChainFormat'}) );
    $command .= " --showLabels no"
	if ( $drawLabels ne "gnuplot" );
    $command .= " --xLabel '".$properties[$i]->{'xLabel'}."'"
        if ( exists($properties[$i]->{'xLabel'}) );
    $command .= " --zLabel '".$properties[$i]->{'zLabel'}."'"
        if ( exists($properties[$i]->{'zLabel'}) );
    if ( @ranges ) {
	foreach my $range ( @ranges ) {
	    $command .= " --range \"".$range."\"";
	}
    }
    # Create job.
    unless ( -e $outputFileName."_".$i.".pdf" ) {
	my %job =
	    (
	     launchFile => $outputFileName."_".$i.".sh",
	     label      => "triangle",
	     logFile    => $outputFileName."_".$i.".log",
	     command    => $command,
	     ppn        => 1
	    );
	foreach ( 'walltime', 'memory' ) {
	    $job{$_} = $arguments{$_}
	    if ( exists($arguments{$_}) );
	}
	push(@jobStack,\%job);
    }
    if ( $i < scalar(@properties)-1 ) { 
	for(my $j=$i+1;$j<scalar(@properties);++$j) {
	    my $command = "constraints/visualization/mcmcVisualize.pl ".$parameterFileName." ".$fileRoot." --workDirectory ".$workDirectory." --yProperty '".$properties[$i]->{'name'}->{'value'}."' --yScale ".$properties[$i]->{'operatorUnaryMapper'}->{'value'}." --xProperty '".$properties[$j]->{'name'}->{'value'}."' --xScale ".$properties[$j]->{'operatorUnaryMapper'}->{'value'}." --textSize ".$textSize." --plotSize ".$plotSize." --lineWeight ".$lineWeight." --labelStyle ".$labelStyle." --output ".$outputFileName."_".$i."_".$j.".pdf --data ".$outputFileName."_".$i."_".$j.".xml ".$options;
	    $command .= " --oldChainFormat ".$arguments{'oldChainFormat'}
	        if ( exists($arguments{'oldChainFormat'}) );
	    $command .= " --showLabels no"
		if ( $drawLabels ne "gnuplot" );
	    $command .= " --xLabel '".$properties[$j]->{'xLabel'}."'"
		if ( exists($properties[$j]->{'xLabel'}) );
	    $command .= " --yLabel '".$properties[$i]->{'xLabel'}."'"
		if ( exists($properties[$i]->{'xLabel'}) );
	    if ( $j == scalar(@properties)-1 ) {
		$command .= " --labels y2 --colorbox 0";
	    } else {
		$command .= " --labels none --colorbox 0";
	    }
	    if ( @ranges ) {
		foreach my $range ( @ranges ) {
		    $command .= " --range \"".$range."\"";
		}
	    }
	    # Create job.
	    unless ( -e $outputFileName."_".$i."_".$j.".pdf" ) {
	    	my %job =
		    (
		     launchFile => $outputFileName."_".$i."_".$j.".sh",
		     label      => "triangle",
		     logFile    => $outputFileName."_".$i."_".$j.".log",
		     command    => $command,
		     ppn        => 1
		    );
		foreach ( 'walltime', 'memory' ) {
		    $job{$_} = $arguments{$_}
		    if ( exists($arguments{$_}) );
		}
		push(@jobStack,\%job);
	    }
	}
    }
}
# Submit jobs.
die("queue manager does not support submitting job arrays")
    unless ( exists($Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}) );
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%arguments,@jobStack);

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
print oHndl "\\newcommand{\\triangledir}{".$shortcut."}\n";
print oHndl "\\renewcommand{\\arraystretch}{0}\n";
print oHndl "\\setlength{\\tabcolsep}{0pt}\n";
if ( $drawLabels eq "gnuplot" ) {
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
	my $shiftHorizontal = -$scale*($width -$standardWidth );
	my $shiftVertical   = -$scale*($height-$standardHeight);
	print oHndl "\\hspace{".$shiftHorizontal."pt}"
	    if ( $i > 0 );
	print oHndl "\\includegraphics[scale=".$scale."]{\\triangledir/".$outputBaseName."_".$i.".pdf}";
	print oHndl "\\vspace{".$shiftVertical."pt}";
	if ( $i < scalar(@properties)-1 ) { 
	    for(my $j=$i+1;$j<scalar(@properties);++$j) {
		print oHndl "&\\raisebox{".$shiftVertical."pt}{\\includegraphics[scale=".$scale."]{\\triangledir/".$outputBaseName."_".$i."_".$j.".pdf}}";
	    }
	}
	print oHndl "\\\\\n";
    }
} elsif ( $drawLabels eq "latex" ) {
    my $width;
    my $height;
    print oHndl "\\begin{tabular}{".("l\@{}c\@{}r\@{}" x scalar(@properties))."}\n";
    for(my $i=0;$i<=scalar(@properties);++$i) {
	unless ( defined($width) ) {
	    open(iHndl,"pdfinfo ".$outputFileName."_".$i.".pdf|");
	    while ( my $line = <iHndl> ) {
		if ( $line =~ m/Page size:\s*(\d+)\s*x\s*(\d+)\s*pts/ ) {
		    $width  = $1;
		    $height = $2;
		}
	    }
	}
	if ( $i > 1 ) {
	    print oHndl ("&&&" x ($i-1));
	}
	if ( $i > 0 ) {
	    my $xml  = new XML::Simple();
	    my $data = $xml->XMLin($outputDirectoryName."/".$outputBaseName."_".($i-1).".xml");
	    my $rangeMinimum = &latexFormat($data->{'x'}->[ 0],$sigFigs);
	    my $rangeMaximum = &latexFormat($data->{'x'}->[-1],$sigFigs);
	    my @columns;
	    my $label = exists($properties[$i-1]->{'label'}) ? "\$".$properties[$i-1]->{'label'}."\$" : $properties[$i-1]->{'xLabel'};
	    my $j = -1;
	    foreach ( $rangeMinimum, $label, $rangeMaximum ) {
		++$j;
		my $content;
		if ( $j == 0 ) {
		    $content .= "\\multicolumn{1}{p{".($width*$scale/3.25)."pt}}{";
		} elsif ( $j == 1 ) {
		    $content .= "\\multicolumn{1}{p{".($width*$scale/3.25)."pt}}{\\hspace{".($width*$scale/3.1/2.0)."pt}";
		} elsif ( $j == 2 ) {
		    $content .= "\\multicolumn{1}{r}{";
		}
		$content .= "\\raisebox{";
		if ( $i == scalar(@properties) ) {
		    $content .= "0pt";
		} else {		
		    $content .= ($height*$scale)."pt-\\widthof{\\".$labelStyle." x}";
		}
		$content .= "-\\widthof{\\".$labelStyle." ".$_."}}[0pt][0pt]{\\rotatebox{90}{\\".$labelStyle." ".$_."}".($j == 2 ? "\\hspace{3pt}" : "")."}";
		$content .= "}";
		push(@columns,$content);	    
	    }
	    print oHndl join("&",@columns);
	    print oHndl "&"
		if ( $i < scalar(@properties) );
	}
	if ( $i < scalar(@properties) ) {
	    my $xml          = new XML::Simple();
	    my $data         = $xml->XMLin($outputDirectoryName."/".$outputBaseName."_".$i.".xml");
	    my $centralValue = +$data->{'modeValueX'}                   ;
	    my $upperOffset  = +$data->{'spanUpperX'}->[0]-$centralValue;
	    my $lowerOffset  = -$data->{'spanLowerX'}->[0]+$centralValue;
	    my $label        = exists($properties[$i]->{'label'}) ? "\$".$properties[$i]->{'label'}."\$" : $properties[$i]->{'xLabel'};
	    open(my $constraint,">".$outputDirectoryName."/".$outputBaseName."_".$i.".tex");
	    print $constraint $label." \$=".&latexFormatErrors($centralValue,$lowerOffset,$upperOffset,$sigFigs,mathMode => 1)."\$";
	    close($constraint);
	    system("pdfcrop --margins \"-1 -1 -1 -1\" ".$outputDirectoryName."/".$outputBaseName."_".$i.".pdf ".$outputDirectoryName."/".$outputBaseName."_".$i."_cropped.pdf")
		unless ( -e $outputDirectoryName."/".$outputBaseName."_".$i."_cropped.pdf" );
	    print oHndl "\\multicolumn{3}{c}{\\includegraphics[scale=".$scale."]{\\triangledir/".$outputBaseName."_".$i."_cropped.pdf}}";
	    if ( $i < scalar(@properties)-1 ) { 
		for(my $j=$i+1;$j<scalar(@properties);++$j) {
		    system("pdfcrop --margins \"-1 -1 -1 -1\" ".$outputDirectoryName."/".$outputBaseName."_".$i."_".$j.".pdf ".$outputDirectoryName."/".$outputBaseName."_".$i."_".$j."_cropped.pdf")
			unless ( -e $outputDirectoryName."/".$outputBaseName."_".$i."_".$j."_cropped.pdf" );
		    print oHndl "&\\multicolumn{3}{c}{\\includegraphics[scale=".$scale."]{\\triangledir/".$outputBaseName."_".$i."_".$j."_cropped.pdf}}";
		}
	    }
	}
	print oHndl "\\\\\n";
    }
}
# Close the output file.
print oHndl "\\end{tabular}\n";
close(oHndl);

exit;

sub latexFormat {
    my $value   = shift();
    my $sigFigs = shift();
    my %options = @_;
    # Handle zero.
    return "0.0"
	if ( $value == 0.0 );
    # Handle negative values.
    my $isNegative = $value < 0.0 ? "-" : "";
    $value = abs($value);
    # Determine order.
    my $o = floor(log($value)/log(10.0));
    my $order = floor(log10($value));
    # Scale value.
    $value /= 10.0**$order;
    # Format value.
    my $format = "%".$sigFigs.".".($sigFigs-1)."f";    
    # Return
    my $result = $isNegative.sprintf($format,$value).($order == 0 ? "" : " \\times 10^{".$order."}");
    $result = "\$".$result."\$"
	unless ( exists($options{'mathMode'}) && $options{'mathMode'} == 1 );
    return $result;
}

sub latexFormatErrors {
    my $value       = shift();
    my $lowerOffset = shift();
    my $upperOffset = shift();
    my $sigFigs     = shift();
    my %options      = @_;
    my $order;
    my $isNegative = "";
    # Handle zero.
    if ( $value == 0.0 ) {
	$order = 0;
    } else {
	# Handle negative values.
	$isNegative = $value < 0.0 ? "-" : "";
	$value = abs($value);
	# Determine order.
	$order = floor(log10($value));
    }
    # Scale value.   
    $value       /= 10.0**$order;
    $lowerOffset /= 10.0**$order;
    $upperOffset /= 10.0**$order;
    $lowerOffset  = $lowerOffset < 0.0 ? 0.0 : $lowerOffset;
    $upperOffset  = $upperOffset < 0.0 ? 0.0 : $upperOffset;
    my $lowerOrder = $lowerOffset == 0.0 ? 0 : floor(log10($lowerOffset));
    my $upperOrder = $upperOffset == 0.0 ? 0 : floor(log10($upperOffset));
    # Format value.
    my $orderLowest  = $lowerOrder < $upperOrder ? $lowerOrder : $upperOrder;
    my $sigFigsValue = $orderLowest < 0 ? $sigFigs-$orderLowest : $sigFigs;    
    my $format       = "%".($sigFigsValue+1).".".($sigFigsValue-1)."f";
    # Return
    my $result = ($order == 0 ? "" : "(").$isNegative.sprintf($format,$value)."^{+".sprintf($format,$upperOffset)."}_{-".sprintf($format,$lowerOffset)."}".($order == 0 ? "" : ") \\times 10^{".$order."}");
    $result = "\$".$result."\$"
	unless ( exists($options{'mathMode'}) && $options{'mathMode'} == 1 );
    return $result;
}
