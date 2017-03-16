#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::MatrixOps;
use PDL::IO::HDF5;
use Data::Dumper;
use UNIVERSAL;
use Galacticus::Constraints::Parameters;
use Galacticus::Constraints::Covariances;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use LaTeX::Format;

# Run calculations of posterior predictive checks on a constrained Galacticus model.
# Andrew Benson (19-March-2013)

# Get arguments.
die("Usage: posteriorPredictiveChecks.pl <configFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make           => "yes",
     plotTitle      => "yes",
     showPValue     => "yes",
     sampleFrom     =>    -1,
     sampleCount    =>    -1,
     outliers       =>    "",
     overSampleRate =>     1
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Validate the config file.
die("posteriorPredictiveChecks.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("posteriorPredictiveChecks.pl: compilation must be specified in config file"   ) unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("posteriorPredictiveChecks.pl: baseParameters must be specified in config file") unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Determine the work directory.
my $workDirectory  = $config->{'likelihood'}->{'workDirectory'};
$arguments{'sampleDirectory'} = $workDirectory."/posteriorSample"
    unless ( exists($arguments{'sampleDirectory'}) );

# Get a hash of the parameter values.
my $compilationFile;
if ( exists($arguments{'compilationOverride'}) ) {
    $compilationFile = $arguments{'compilationOverride'};
} else {
    $compilationFile = $config->{'likelihood'}->{'compilation'};
}
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($compilationFile,$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Generate a sample of models.
print "Generating posterior model sample...\n";
(my $sampleCount, my $sampleDirectory) = &Galacticus::Constraints::Parameters::Sample_Models($config,\%arguments);
print "...done\n";

# Begin building LaTeX table of p-values.
system("mkdir -p ".$workDirectory."/posteriorPredictiveChecks");
open(my $pValueTable,">".$workDirectory."/posteriorPredictiveChecks/ppcPValues.tex");
print $pValueTable "\\begin{tabular}{ll}\n";
print $pValueTable "\\hline\n";
print $pValueTable "{\\bf Constraint} & \\boldmath{\$p_{\\rm B}\$} \\\\\n";
print $pValueTable "\\hline\n";
# Read in the results for all models and constraints.
my $xml = new XML::Simple;
# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $xml = new XML::Simple;
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    print "Posterior predictive check for constraint '".$constraintDefinition->{'label'}."'...\n"; 
    # Iterate over models.
    my @modelX;
    my @modelY;
    my $dataY;
    my $dataCovariance;
    for (my $i=0;$i<$sampleCount;++$i) {
	# Specify output directory.
	my $modelDirectory = $sampleDirectory."/".$i."/";
	# Parse the results file.
	my $resultFile         = $modelDirectory."/".$constraintDefinition->{'label'}.".hdf5";
	my $results            = new PDL::IO::HDF5($resultFile);
	my $x                  = $results->dataset('x'             )->get();
	my $y                  = $results->dataset('y'             )->get();
	$dataY                 = $results->dataset('yData'         )->get();
	my $dataCovarianceFlat = $results->dataset('covarianceData')->get();
	my $covarianceFlat     = $results->dataset('covariance'    )->get();
	my $dataSize           = nelem($y);
	$dataCovariance        = reshape($dataCovarianceFlat,$dataSize,$dataSize);
	# Compute the square root of the covariance matrix using eigendecomposition. This is more robust than Cholesky
	# decomposition in this case as it allows us to zero out any negative eigenvalues.
	(my $eigenVectors, my $eigenValues)                = eigens_sym($dataCovariance        );
	my $eigenValuesMatrix                              = zeroes    ($dataCovariance        );
	my $nonNegative                                    = which     ($eigenValues     >= 0.0);
	$eigenValuesMatrix->diagonal(0,1)->($nonNegative) .= $eigenValues->($nonNegative)->sqrt();
	my $dataCovarianceSquareRoot                       = $eigenVectors x $eigenValuesMatrix;
	# Build a realization of the data from the model and the data covariance.
	for(my $i=0;$i<$arguments{'overSampleRate'};++$i) {
	    my $randomDeviates     = pdl grandom($dataSize);
	    my $perturbations      = $dataCovarianceSquareRoot x transpose($randomDeviates);
	    my $yPerturbed         = $y+$perturbations->((0),:);
	    # Store the results.
	    push(@modelY,$yPerturbed);
	    push(@modelX,$x         );
	}
    }
    # Evaluate the mean result and covariance.
    my $meanY      = pdl zeroes(nelem($modelY[0]));
    my $covariance = pdl zeroes(nelem($modelY[0]),nelem($modelY[0]));
    for (my $i=0;$i<scalar(@modelY);++$i) {
	$meanY += $modelY[$i];
    }
    $meanY /= scalar(@modelY);
    for (my $i=0;$i<scalar(@modelY);++$i) {
	my $offset = $modelY[$i]-$meanY;
	$covariance += outer($offset,$offset);
    }
    $covariance /= scalar(@modelY)-1;
    my $inverseCovariance = minv($covariance);
    # Compute the test statistic.
    my $testStatistic = pdl [];
    for (my $i=0;$i<scalar(@modelY);++$i) {
	my $difference = $modelY[$i]-$meanY;
	$testStatistic = $testStatistic->append(sclr($difference x $inverseCovariance x transpose($difference)));
    }
    # Compute test statistic for the data.
    my $difference        = $dataY-$meanY;
    my $testStatisticData = sclr($difference x $inverseCovariance x transpose($difference));
    # Find the p-value.
    my $exceeders = which($testStatistic > $testStatisticData);
    my $pValue    = nelem($exceeders)/nelem($testStatistic);
	my $pValueLabel;
    if ( $pValue == 0.0 ) {
	my $pValueLimit = 1.0/nelem($testStatistic);
	$pValueLabel = "<".($pValueLimit < 1.0e-2 ? &LaTeX::Format::Number($pValueLimit,2, mathMode => 0) : sprintf("%4.2f",$pValueLimit));
    } else { 
	$pValueLabel = $pValue < 1.0e-2 ? &LaTeX::Format::Number($pValue,2,mathMode => 0) : sprintf("%4.2f",$pValue);
    }
    print "   → p_B = ".$pValue."\n";
    print $pValueTable $constraintDefinition->{'name'}." & \$".$pValueLabel."\$ \\\\\n";
    open(oHndl,">".$workDirectory."/posteriorPredictiveChecks/".$constraintDefinition->{'label'}."_testStatistic_pValue.txt");
    print oHndl "Bayesian p-value                           : ".$pValue           ."\n";
    print oHndl "Test statistic, data                       : ".$testStatisticData."\n";
    print oHndl "Test statistic realization number and value:\n";
    for(my $i=0;$i<nelem($testStatistic);++$i) {
	print oHndl "\t".$i."\t".$testStatistic->(($i))."\n";
    }
    close(oHndl);
    # Create a plot showing the test statistic distribution.
    my ($gnuPlot, $plot);
    system("mkdir -p ".$workDirectory."/posteriorPredictiveChecks");
    (my $safeLabel = $constraintDefinition->{'label'}) =~ s/\./_/g;
    my $plotFileName = $workDirectory."/posteriorPredictiveChecks/".$safeLabel."_testStatistic.pdf";
    $plotFileName =~ s/_pdf$/.pdf/;
    (my $plotFileTeX = $plotFileName) =~ s/\.pdf$/.tex/;
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.90\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.2,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    my $xAll     = $testStatistic->append($testStatisticData);
    my $xMinimum = minimum($xAll);
    my $xMaximum = maximum($xAll);
    my $decadesSpanned = log10($xMaximum/$xMinimum);
    if ( $decadesSpanned > 2.0 ) {
    	print $gnuPlot "set logscale x\n";
    	print $gnuPlot "set mxtics 10\n";
    	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    }
    print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
    print $gnuPlot "set yrange [-0.05:1.05]\n";
    print $gnuPlot "set title '".$constraintDefinition->{'name'}."'\n"
    	if ( $arguments{'plotTitle'} eq "yes" );
    print $gnuPlot "set xlabel 'Test statistic; \$\\mathcal{T}\$'\n";
    print $gnuPlot "set ylabel 'Cumulative probability; \$ P ( < \\mathcal{T})\$'\n";
    print $gnuPlot "set arrow from first ".$testStatisticData.", 0.0 to first ".$testStatisticData.", 0.4 filled linewidth 5 linecolor rgbcolor \"#3CB371\"\n";
    if ( $arguments{'showPValue'} eq "yes" ) {
    	print $gnuPlot "set label '\$ p_\\mathrm{B}".($pValue == 0.0 ? "" : "=").$pValueLabel."\$' at graph 0.05, 0.8\n";
    }
    my $testStatisticSorted    = $testStatistic->qsort();
    my $cumulativeProbability  = pdl sequence(nelem($testStatistic));
    $cumulativeProbability    /= nelem($testStatistic);
    &GnuPlot::PrettyPlots::Prepare_Dataset
    	(
    	 \$plot,
    	 $testStatisticSorted,
    	 $cumulativeProbability,
    	 style       => "line",
    	 weight      => [5,3],
    	 color       => $GnuPlot::PrettyPlots::colorPairs{'redYellow'}
    	);
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);    
}
print $pValueTable "\\hline\n";
print $pValueTable "\\end{tabular}\n";
close($pValueTable);

exit;
