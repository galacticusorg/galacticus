#!/usr/bin/env perl
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use XML::Simple;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::Local;
use Galacticus::Launch::PBS;
use Galacticus::Launch::MonolithicPBS;
use Galacticus::Launch::Slurm;
use Galacticus::Constraints::Parameters;
use List::ExtraUtils;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Visualize the posterior distribution from MCMC chains using a triangle arrangement.
# Andrew Benson (12-June-2012)

# Ensure large objects are allowed.
$PDL::BIGPDL = 1;
    
# Get file name to process.
die("Usage: posteriorCornerPlot.pl chainFileRoot mcmcParameterFile [options]")
    unless ( scalar(@ARGV) > 1 );
my $chainFileRoot      = $ARGV[0];
my $mcmcConfigFileName = $ARGV[1];
my %options =
    (
     workDirectory  => ".",
     useUnconverged => "no",
     chain          => "all",
     outputFile     => "corner.pdf",
     plotSize       => "8.5in,11in",
     xBorderLower   => 0.10,
     xBorderUpper   => 0.02,
     yBorderLower   => 0.05,
     yBorderUpper   => 0.02,
     lineWeight1D   => "4,2",
     lineWeight2D   => "8,4",
     stepsMaximum   => 1000,
     countGrid      => 50,
     sigFigs        => 2,
     cleanUp        => "no"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Create work directory.
system("mkdir -p ".$options{'workDirectory'});

# Parse property options.
if ( exists($options{'property'}) ) {
    my @properties;
    foreach my $property ( &List::ExtraUtils::as_array($options{'property'}) ) {
	my @propertySpecifiers = split(/(?<!:):(?!:)/,$property);
	my %propertyData = 
	    (
	     name    => $propertySpecifiers[0],
	     scaling => "linear"
	    );
	for(my $i=1;$i<scalar(@propertySpecifiers);++$i) {
	    $propertyData{'scaling'} = "log"
		if ( $propertySpecifiers[$i] eq "logarithmic" );
	    if ( $propertySpecifiers[$i] =~ m/^(label|format|increment)=(.*)/ ) {
		$propertyData{$1} = $2;
	    }
	}
	push(@properties,\%propertyData);
    }
    @{$options{'property'}} = @properties;
}

# Parse lineWeight options.
my @lineWeight1D = split(/,/,$options{'lineWeight1D'});
my @lineWeight2D = split(/,/,$options{'lineWeight2D'});

# Open the parameter file and parse the available parameter names.
my $xml        = new XML::Simple;
my $mcmcConfig = $xml->XMLin($mcmcConfigFileName);
my @propertyNamesAvailable;
foreach my $parameter ( @{$mcmcConfig->{'posteriorSampleSimulation'}->{'modelParameter'}} ) {
    push(@propertyNamesAvailable,$parameter->{'name'}->{'value'})
	 if ( $parameter->{'value'} eq "active" );
}
my @propertyNames;
if ( exists($options{'property'}) ) {
    @propertyNames = map {$_->{'name'}} @{$options{'property'}};
    die("posteriorCornerPlot.pl: at least 2 propertyNames must be specified")
	if ( scalar(@propertyNames) < 2 );
    foreach my $property ( @propertyNames ) {
	die("Property '".$property."' is not available") unless ( grep {$_ eq $property} @propertyNamesAvailable );
    }
} else {
    @propertyNames = @propertyNamesAvailable;
}

# Construct property list, determining the index of each parameter of interest in the state matrix.
my @properties;
my $indexParameter = -1;
foreach my $parameter ( @{$mcmcConfig->{'posteriorSampleSimulation'}->{'modelParameter'}} ) {
    ++$indexParameter;
    if ( grep {$parameter->{'name'}->{'value'} eq $_} @propertyNames ) {
	foreach my $property ( @{$options{'property'}} ) {
	    if ( $property->{'name'} eq $parameter->{'name'}->{'value'} ) {
		foreach ( 'label', 'scaling', 'format', 'increment' ) {
		    $parameter->{$_} = $property->{$_}
		        if ( exists($property->{$_}) );
		}
		$parameter->{'index'} = $indexParameter;
	    }
	}
	push(@properties,$parameter);
    }
}

# Retrieve posterior samples.
my $posteriorSamples = &Galacticus::Constraints::Parameters::parameterMatrix($mcmcConfig,$options{'chain'},\%options);

# Map properties to logarithmic form as needed.
for(my $i=0;$i<scalar(@properties);++$i) {
    $posteriorSamples->(($i),:) .= $posteriorSamples->(($i),:)->log()
	if ( $properties[$i]->{'scaling'} eq "log" );
}

# Create a list of temporary files.
my @tempFiles;

# Perform KDE to find posterior distributions.
my $posteriorX1D         = pdl zeros(scalar(@properties)                    ,$options{'countGrid'}                      );
my $posteriorDensities1D = pdl zeros(scalar(@properties)                    ,$options{'countGrid'}                      );
my $posteriorX2D         = pdl zeros(scalar(@properties)                    ,$options{'countGrid'}                      );
my $posteriorY2D         = pdl zeros(scalar(@properties)                    ,$options{'countGrid'}                      );
my $posteriorDensities2D = pdl zeros(scalar(@properties),scalar(@properties),$options{'countGrid'},$options{'countGrid'});
for(my $i=0;$i<scalar(@properties);++$i) {
    # 1D posteriors.
    my $samples1DFileName = $options{'workDirectory'}."/samples1D_".$i.".txt";
    my $kde1DFileName     = $options{'workDirectory'}."/samples1D_".$i.".kde";
    push(@tempFiles,$samples1DFileName,$kde1DFileName);
    open(my $samples1D,">",$samples1DFileName);
    for(my $k=0;$k<$posteriorSamples->dim(1);++$k) {
	print $samples1D $posteriorSamples->(($properties[$i]->{'index'}),($k))."\n";
    }
    close($samples1D);
    system("python2 constraints/visualization/kernelDensityEstimation.py ".$samples1DFileName." 0 --ngrid=".$options{'countGrid'}." --output=".$kde1DFileName);
    open(my $kde1D,"<",$kde1DFileName);
    for(my $m=0;$m<$options{'countGrid'};++$m) {
	my $line = <$kde1D>;
	my @columns = split(" ",$line);
	$posteriorX1D        ->(($i),($m)) .= $columns[0];
	$posteriorDensities1D->(($i),($m)) .= $columns[1];
    }
    close($kde1D);
    # 2D posteriors.
    for(my $j=$i+1;$j<scalar(@properties);++$j) {
	my $samples2DFileName = $options{'workDirectory'}."/samples2D_".$i."_".$j.".txt";
	my $kde2DFileName     = $options{'workDirectory'}."/samples2D_".$i."_".$j.".kde";
	push(@tempFiles,$samples2DFileName,$kde2DFileName);
	open(my $samples2D,">",$samples2DFileName);
	for(my $k=0;$k<$posteriorSamples->dim(1);++$k) {
	    print $samples2D $posteriorSamples->(($properties[$i]->{'index'}),($k))." ".$posteriorSamples->(($properties[$j]->{'index'}),($k))."\n";
	}
	close($samples2D);	
	system("python2 constraints/visualization/kernelDensityEstimation.py ".$samples2DFileName." 0 1 --ngrid=".$options{'countGrid'}." --output=".$kde2DFileName);
	open(my $kde2D,"<",$kde2DFileName);
	for(my $m=0;$m<$options{'countGrid'};++$m) {
	    for(my $n=0;$n<$options{'countGrid'};++$n) {
		my $line = <$kde2D>;
		my @columns = split(" ",$line);
		$posteriorX2D        ->(($i)     ,($m)     ) .= $columns[0];
		$posteriorY2D        ->(     ($j)     ,($n)) .= $columns[1];
		$posteriorDensities2D->(($i),($j),($m),($n)) .= $columns[2];
	    }
	}
	close($kde2D);
    }
}

# Create the corner plot.
my $gnuPlot;
(my $plotFileTeX = $options{'outputFile'}) =~ s/\.pdf$/.tex/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size ".$options{'plotSize'}."\n";
print $gnuPlot "set output '".$plotFileTeX."'\n";
print $gnuPlot "set multiplot\n";
for(my $i=0;$i<scalar(@properties);++$i) {
    for(my $j=$i;$j<scalar(@properties);++$j) {
	# Compute position for this panel.
	my $lmargin = sprintf("%4.2f",(+$i                    +0)*(1.0-$options{'xBorderLower'}-$options{'xBorderUpper'})/scalar(@properties)+$options{'xBorderLower'});
	my $rmargin = sprintf("%4.2f",(+$i                    +1)*(1.0-$options{'xBorderLower'}-$options{'xBorderUpper'})/scalar(@properties)+$options{'xBorderLower'});
	my $bmargin = sprintf("%4.2f",(-$j+scalar(@properties)-1)*(1.0-$options{'yBorderLower'}-$options{'yBorderUpper'})/scalar(@properties)+$options{'yBorderLower'});
	my $tmargin = sprintf("%4.2f",(-$j+scalar(@properties)+0)*(1.0-$options{'yBorderLower'}-$options{'yBorderUpper'})/scalar(@properties)+$options{'yBorderLower'});
	print $gnuPlot "set lmargin screen ".$lmargin."\n";
	print $gnuPlot "set rmargin screen ".$rmargin."\n";
	print $gnuPlot "set bmargin screen ".$bmargin."\n";
	print $gnuPlot "set tmargin screen ".$tmargin."\n";
	print $gnuPlot "set border; set xtics; set ytics; unset label\n";
	# Handle 1D and 2D posterior panels separately.
	if ( $i == $j ) {
	    # 1D panel.
	    print $gnuPlot "unset ylabel\n";     # Never label the y-axis for 1D panels.	  
	    print $gnuPlot "set format y ''\n";
	    print $gnuPlot "unset logscale y\n"; # y-axis is always linear for 1D panels.
	    print $gnuPlot "unset mytics\n";
	    print $gnuPlot "set ytics autofreq\n";
	    if ( $i == scalar(@properties)-1 ) {
		if ( $properties[$i]->{'scaling'} eq "log" ) {
		    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
		} else {
		    print $gnuPlot "set format x '".(exists($properties[$i]->{'format'}) ? $properties[$i]->{'format'} : "\%g")."'\n";
		}
		print $gnuPlot "set xlabel '".$properties[$i]->{'label'}."'\n";
	    } else {
	        print $gnuPlot "set format x ''\n";
		print $gnuPlot "unset xlabel\n";
	    }
	    if ( $properties[$i]->{'scaling'} eq "log" ) {
		print $gnuPlot "set logscale x\n";
		print $gnuPlot "set mxtics 10\n";
	    } else {
		print $gnuPlot "unset logscale x\n";
		print $gnuPlot "unset mxtics\n";
	    }
	    print $gnuPlot "set xrange [".$posteriorX1D->(($i),:)->minimum().":".$posteriorX1D->(($i),:)->maximum()."]\n";
	    print $gnuPlot "set yrange [0.0:".$posteriorDensities1D->(($i),:)->maximum()."]\n";
	    print $gnuPlot "set pointsize 1.0\n";
	    if ( exists($properties[$i]->{'increment'}) ) {
		my $xticsStart = int($posteriorX2D->(($i),:)->minimum()/$properties[$i]->{'increment'}    )*$properties[$i]->{'increment'};
		my $xticsEnd   = int($posteriorX2D->(($i),:)->maximum()/$properties[$i]->{'increment'}+1.0)*$properties[$i]->{'increment'};
		print $gnuPlot "set xtics ".$xticsStart.",".$properties[$i]->{'increment'}.",".$xticsEnd."\n";
	    } else {
		print $gnuPlot "set xtics autofreq\n";
	    }
	    # Determine contour levels for confidence regions.
	    my $probabilityDifferential                       = $posteriorDensities1D->(($i),:)->flat()->qsort();
	    my $probabilityTotal                              = $probabilityDifferential->sum  ();
	    my $probabilityCumulative                         = $probabilityDifferential->cumusumover()/$probabilityTotal;
	    my $confidence                                    = pdl ( 0.682689492137 );
	    my $excludedFraction                              = 1.0-$confidence;
	    (my $levels, my $error)                           = interpolate($excludedFraction,$probabilityCumulative,$probabilityDifferential);
	    my $posteriorDensities1DConfidence                = $posteriorDensities1D->(($i),:)->copy();
	    (my $confident, my $notConfident)                 = which_both($posteriorDensities1DConfidence >= $levels->((0)));
	    $posteriorDensities1DConfidence->($notConfident) .= 0.0;
	    my $zeros                                         = pdl zeros(nelem($posteriorDensities1DConfidence));
	    my $indexMode                                     = $posteriorDensities1D->(($i),:)->maximum_ind();
	    my $mode                                          = $posteriorX1D->(($i),($indexMode));
	    my $confidenceUpper                               = +$posteriorX1D->(($i),$confident)->maximum()-$mode;
	    my $confidenceLower                               = -$posteriorX1D->(($i),$confident)->minimum()+$mode;
	    print $gnuPlot
		"set title offset 0.0,-0.5'".$properties[$i]->{'label'}." \$ = ".
		&latexFormatErrors($mode,$confidenceLower,$confidenceUpper,$options{'sigFigs'},mathMode => 1).
		" \$'\n";
	    my $plot;
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
	       \$plot                                                         ,
		$posteriorX1D                  ->(($i),:)                     ,
		$posteriorDensities1DConfidence                               ,
		y2          => $zeros,                                        ,
		style       => "filledCurve"                                  ,
		weight      => \@lineWeight1D                                 ,
		linePattern => 1                                              ,
		color       => ["#38A8FF","#38A8FF"]
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                        ,
		$posteriorX1D        ->(($i),:)                               ,
		$posteriorDensities1D->(($i),:)                               ,
		style        => "line"                                        ,
		weight       => \@lineWeight1D                                ,
		linePattern  => 1                                             ,
		color        => $GnuPlot::PrettyPlots::colorPairs{'slateGray'}
		);
	    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot, multiPlot => 1);
	} else {
	    # 2D panel.
	    if ( $properties[$i]->{'scaling'} eq "log" ) {
		print $gnuPlot "set logscale x\n";
		print $gnuPlot "set mxtics 10\n";
	    } else {
		print $gnuPlot "unset logscale x\n";
		print $gnuPlot "unset mxtics\n";
	    }
	    if ( $properties[$j]->{'scaling'} eq "log" ) {
		print $gnuPlot "set logscale y\n";
		print $gnuPlot "set mytics 10\n";
	    } else {
		print $gnuPlot "unset logscale y\n";
		print $gnuPlot "unset mytics\n";
	    }
	    if ( exists($properties[$i]->{'increment'}) ) {
		my $xticsStart = int($posteriorX2D->(($i),:)->minimum()/$properties[$i]->{'increment'}    )*$properties[$i]->{'increment'};
		my $xticsEnd   = int($posteriorX2D->(($i),:)->maximum()/$properties[$i]->{'increment'}+1.0)*$properties[$i]->{'increment'};
		print $gnuPlot "set xtics ".$xticsStart.",".$properties[$i]->{'increment'}.",".$xticsEnd."\n";
	    } else {
		print $gnuPlot "set xtics autofreq\n";
	    }
	    if ( exists($properties[$j]->{'increment'}) ) {
		my $yticsStart = int($posteriorY2D->(($j),:)->minimum()/$properties[$j]->{'increment'}    )*$properties[$j]->{'increment'};
		my $yticsEnd   = int($posteriorY2D->(($j),:)->maximum()/$properties[$j]->{'increment'}+1.0)*$properties[$j]->{'increment'};
		print $gnuPlot "set ytics ".$yticsStart.",".$properties[$i]->{'increment'}.",".$yticsEnd."\n";
	    } else {
		print $gnuPlot "set ytics autofreq\n";
	    }
	    print $gnuPlot "unset title\n";
	    print $gnuPlot "set xrange [".$posteriorX2D->(($i),:)->minimum().":".$posteriorX2D->(($i),:)->maximum()."]\n";
	    print $gnuPlot "set yrange [".$posteriorY2D->(($j),:)->minimum().":".$posteriorY2D->(($j),:)->maximum()."]\n";
	    print $gnuPlot "set pointsize 1.0\n";
	    # Determine contour levels for confidence regions.
	    my $probabilityDifferential = $posteriorDensities2D->(($i),($j),:,:)->flat()->qsort();
	    my $probabilityTotal        = $probabilityDifferential->sum  ();
	    my $probabilityCumulative   = $probabilityDifferential->cumusumover()/$probabilityTotal;
	    my $confidence              = pdl ( 0.682689492137,  0.954499736104 );
	    my @confidenceColor         = ( "#384BFF", "#38A8FF" );
	    my $excludedFraction        = 1.0-$confidence;
	    (my $levels, my $error)     = interpolate($excludedFraction,$probabilityCumulative,$probabilityDifferential);
	    # Generate confidence contours.
	    for(my $level=0;$level<nelem($levels);++$level) {
		open(ppHndl,"|gnuplot");
		print ppHndl "set table 'contour.dat'\n";
		print ppHndl "unset surface\n";
		print ppHndl "set contour base; set cntrparam levels discrete ".$levels->(($level))."\n";
		print ppHndl "splot '-'\n";
		for(my $iGrid=0;$iGrid<$options{'countGrid'};++$iGrid) {
		    for(my $jGrid=0;$jGrid<$options{'countGrid'};++$jGrid) {
			print ppHndl $posteriorX2D->(($i),($iGrid))." ".$posteriorY2D->(($j),($jGrid))." ".$posteriorDensities2D->(($i),($j),($iGrid),($jGrid))."\n";
		    }
		    print ppHndl "\n"
			unless ( $iGrid == $options{'countGrid'}-1 );
		}
		print ppHndl "e\n";
		print ppHndl "unset table\n";
		close(ppHndl);
		system("awk \"NF<2{printf\\\"\\n\\\"}{print}\" < contour.dat > contour".$level.".dat");
	    }
	    # Plot heat map - do not include axes here.
	    print $gnuPlot "unset xlabel\n";
	    print $gnuPlot "unset ylabel\n";
	    print $gnuPlot "set format x ''\n";
	    print $gnuPlot "set format y ''\n";
	    print $gnuPlot "set pm3d map\n";
	    print $gnuPlot "set pm3d explicit\n";
	    print $gnuPlot "set pm3d corners2color c1\n";
	    print $gnuPlot "unset colorbox\n";
	    print $gnuPlot "set palette defined (0 'light-blue', 1 'dark-blue')\n";
	    print $gnuPlot "splot '-' with pm3d notitle\n";
	    my $x = $posteriorX2D->(($i),:);
	    my $y = $posteriorY2D->(($j),:);
	    $x .= exp($x)
		if ( $properties[$i]->{'scaling'} eq "log" );
	    $y .= exp($y)
		if ( $properties[$j]->{'scaling'} eq "log" );
	    for(my $iGrid=0;$iGrid<$options{'countGrid'};++$iGrid) {
		for(my $jGrid=0;$jGrid<$options{'countGrid'};++$jGrid) {
		    print $gnuPlot $x->(($iGrid))." ".$y->(($jGrid))." ".$posteriorDensities2D->(($i),($j),($iGrid),($jGrid))."\n";
		}
		print $gnuPlot "\n"
		    unless ( $iGrid == $options{'countGrid'}-1 );
	    }	    
	    print $gnuPlot "e\n";
	    # Set axes for contour plots.
	    if ( $j == scalar(@properties)-1 ) {
		if ( $properties[$i]->{'scaling'} eq "log" ) {
		    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
		} else {
		    print $gnuPlot "set format x '".(exists($properties[$i]->{'format'}) ? $properties[$i]->{'format'} : "\%g")."'\n";
		}
		print $gnuPlot "set xlabel '".$properties[$i]->{'label'}."'\n";
	    } else {
	        print $gnuPlot "set format x ''\n";
		print $gnuPlot "unset xlabel\n";
	    }
	    if ( $i == 0 ) {
		if ( $properties[$j]->{'scaling'} eq "log" ) {
		    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
		} else {
		    print $gnuPlot "set format y '".(exists($properties[$j]->{'format'}) ? $properties[$j]->{'format'} : "\%g")."'\n";
		}
		print $gnuPlot "set ylabel '".$properties[$j]->{'label'}."'\n";
	    } else {
	        print $gnuPlot "set format y ''\n";
		print $gnuPlot "unset ylabel\n";
	    }
	    # Add contours.
	    print $gnuPlot "plot";
	    for(my $level=0;$level<nelem($levels);++$level) {
		print $gnuPlot ","
		    unless ( $level == 0 );
		print $gnuPlot "  'contour".$level.".dat' with line lw ".$lineWeight2D[0]." lc rgbcolor \"#FFFFFF\" notitle";
		print $gnuPlot ", 'contour".$level.".dat' with line lw ".$lineWeight2D[1]." lc rgbcolor \"".$confidenceColor[$level]."\" notitle";
	    }
	    print $gnuPlot "\n";
	}
    }
}
print $gnuPlot "unset multiplot\n";
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);

# Clean up temporary files.
unlink(@tempFiles)
    if ( $options{'cleanUp'} eq "yes" );

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
