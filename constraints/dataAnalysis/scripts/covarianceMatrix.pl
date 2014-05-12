#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
use LaTeX::Encode;
use Clone qw(clone);
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Find the maximum likelihood estimate of the covariance matrix for a mass function.
# Andrew Benson (05-July-2012)

# Get the config file.
die("Usage: covarianceMatrix.pl <configFile> [<redshiftIndex>]")
    unless ( scalar(@ARGV) == 1 || scalar(@ARGV) == 2 );
my $configFile    = $ARGV[0];
my $redshiftIndex;
$redshiftIndex    = $ARGV[1]
    if ( scalar(@ARGV) == 2 );
my $xml        = new XML::Simple;
my $config     = $xml->XMLin($configFile);

# Validate the config file.
die("covarianceMatrix.pl: config file must specify baseDirectory" )
    unless ( exists($config->{'baseDirectory' }) );
die("covarianceMatrix.pl: config file must specify parameterFile" )
    unless ( exists($config->{'parameterFile' }) );
die("covarianceMatrix.pl: config file must specify constraintFile")
    unless ( exists($config->{'constraintFile'}) );
die("covarianceMatrix.pl: config file must specify mcmcConfigFile")
    unless ( exists($config->{'mcmcConfigFile'}) );

# Set any defaults.
my $pbsLabel = "covariance";
$pbsLabel = $config->{'pbsLabel'}
   if ( exists($config->{'pbsLabel'}) );
my $sourceLabel = "data";
$sourceLabel = latex_encode($config->{'sourceLabel'})
    if ( exists($config->{'sourceLabel'}) );
my $massType = "Stellar";
$massType = $config->{'massType'}
    if ( exists($config->{'massType'}) );
my $massVariable = "M_\\star";
$massVariable = $config->{'massVariable'}
    if ( exists($config->{'massVariable'}) );

# Specify the work directory.
my $baseDirectory = $config->{'baseDirectory'};

# Extract leaf and path for parameter and constraint files.
my $parameterFile  = $config->{'parameterFile' };
my $constraintFile = $config->{'constraintFile'};
my $mcmcConfigFile = $config->{'mcmcConfigFile'};
(my $parameterFileLeaf  = $parameterFile ) =~ s/^.*\/([^\/]+)$/$1/;
(my $parameterFilePath  = $parameterFile ) =~ s/^(.*\/).*$/$1/;
(my $constraintFileLeaf = $constraintFile) =~ s/^.*\/([^\/]+)$/$1/;
(my $constraintFilePath = $constraintFile) =~ s/^(.*\/).*$/$1/;
(my $mcmcConfigFileLeaf = $mcmcConfigFile) =~ s/^.*\/([^\/]+)$/$1/;
(my $mcmcConfigFilePath = $mcmcConfigFile) =~ s/^(.*\/).*$/$1/;

# Ensure that the covariance matrix, mass function, and constraint codes are built.
system("make Mass_Function_Covariance.exe Conditional_Mass_Function.exe Constrain_Galacticus.exe");
die("covarianceMatrix.pl: unable to build executables")
    unless ( $? == 0 );

# List that will be used to store the best fit values from one stage to the next.
my @bestFit;

# Get parameter names
my $xmlNames       = new XML::Simple(KeyAttr => []);
my $mcmcConfig     = $xmlNames->XMLin($mcmcConfigFile);
my @parameterNames = map {$_->{'name'}} @{$mcmcConfig->{'parameters'}->{'parameter'}};

# Perform a specified number of stages of the algorithm.
my $stageCount = 4;
$stageCount = $config->{'stageCount'}
   if ( exists($config->{'stageCount'}) );
for(my $stage=0;$stage<=$stageCount;++$stage) {
    # Create an output folder for this stage.
    my $stageDirectory = $baseDirectory;
    $stageDirectory   .= "redshift".$redshiftIndex."_"
	if ( defined($redshiftIndex) );
    $stageDirectory   .= "stage".$stage;
    system("mkdir -p ".$stageDirectory);
    # Generate the covariance matrix file.
    my $covarianceMatrixFile = $stageDirectory."/".$constraintFileLeaf;
    unless ( -e $covarianceMatrixFile ) {
	# Parse the basic covariance matrix parameter file.
	my $xml        = new XML::Simple;
	my $parameters = $xml->XMLin($parameterFile);
	# Modify parameters as required.
	$parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'} = $covarianceMatrixFile;
	my $include = "true";
	$include = "false"
	    if ( $stage == 0 );
	$parameters->{'parameter'}->{$_}->{'value'} = $include
	    foreach ( "massFunctionCovarianceIncludeHalo", "massFunctionCovarianceIncludeLSS" );
	# Modify redshift ranges.
	if ( defined($redshiftIndex) ) {
	    $parameters->{'parameter'}->{$parameters->{'parameter'}->{'surveyGeometryMethod'}->{'value'}.'RedshiftBin'}->{'value'} = $redshiftIndex;
	    $parameters->{'parameter'}->{'conditionalMassFunctionRedshiftMinimum'}->{'value'} = ${$config->{'redshift'}}[$redshiftIndex]->{'minimum'};
	    $parameters->{'parameter'}->{'conditionalMassFunctionRedshiftMaximum'}->{'value'} = ${$config->{'redshift'}}[$redshiftIndex]->{'maximum'};
	    $parameters->{'parameter'}->{'massFunctionCovarianceRedshiftMinimum' }->{'value'} = ${$config->{'redshift'}}[$redshiftIndex]->{'minimum'};
	    $parameters->{'parameter'}->{'massFunctionCovarianceRedshiftMaximum' }->{'value'} = ${$config->{'redshift'}}[$redshiftIndex]->{'maximum'};
	}
	# For stages after the first, insert the best fit parameters from the previous stage.
	if ( $stage > 0 ) {
	    for(my $i=0;$i<scalar(@parameterNames);++$i) {
		$parameters->{'parameter'}->{"conditionalMassFunctionBehroozi".ucfirst($parameterNames[$i])}->{'value'} = $bestFit[$i];
	    }
	}
	# Write out the new parameter file.
	my $covarianceMatrixParameterFile = $stageDirectory."/".$parameterFileLeaf;
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
	open(oHndl,">".$covarianceMatrixParameterFile);
	print oHndl $xmlOutput->XMLout($parameters);
	close(oHndl);
	# Generate the covariance matrix.
	my $command = "mpirun -np 1 -hostfile \$PBS_NODEFILE ".$galacticusPath."constraints/dataAnalysis/scripts/generateCovarianceMatrix.pl ".$stageDirectory."/".$parameterFileLeaf." ".$configFile;
	$command .= " ".$redshiftIndex
	    if ( defined($redshiftIndex) );
	my $matrixJob = 
	{
	    batchFile  => $stageDirectory."/generateMatrix.pbs",
	    queue      => "batch",
	    nodes      => "nodes=1:ppn=12",
	    wallTime   => "200:00:00",
	    outputFile => $stageDirectory."/generateMatrix.log",
	    name       => $pbsLabel."Stage".$stage."MatrixGen",
	    commands   => $command
	};
	&Submit_To_PBS($matrixJob);
	# Since we have a new covariance matrix file, remove the MCMC chains to force it to be regenerated.
	unlink($stageDirectory."/chains_*.log");
	# Also remove any plot of the correlation matrix.
	unlink($stageDirectory."/correlationMatrix.pdf");
	# For stages past the 0th, determine how well converged the matrix is.
	if ( $stage > 0 ) {
	    my $oldStageDirectory = $baseDirectory;
	    $oldStageDirectory   .= "redshift".$redshiftIndex."_"
		if ( defined($redshiftIndex) );
	    $oldStageDirectory   .= "stage".($stage-1);
	    my $oldCovarianceMatrixFile = $oldStageDirectory."/".$constraintFileLeaf;
	    my $new = new PDL::IO::HDF5(">".$covarianceMatrixFile);
	    my $old = new PDL::IO::HDF5( $oldCovarianceMatrixFile);
	    my $newMatrix = $new->dataset("covariance")->get();
	    my $oldMatrix = $old->dataset("covariance")->get();
	    my $error     = abs($newMatrix-$oldMatrix)/0.5/($newMatrix+$oldMatrix);
	    my $errorMaximum = max($error);
	    $new->dataset("covarianceError")->set($error);
	    $new->attrSet(covarianceErrorMaximum => $errorMaximum);
	}
    }
    # Generate a plot of the correlation matrix.
    unless ( -e $stageDirectory."/correlationMatrix.pdf" ) {
	# Read the covariance matrix.
	my $hdfFile      = new PDL::IO::HDF5($covarianceMatrixFile);
	my $correlation  = $hdfFile->dataset("correlation")->get();
	my $mass         = $hdfFile->dataset("mass"       )->get();
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	my $rangeLow  = $mass->(( 0));
	my $rangeHigh = $mass->((-1));
	$plotFileEPS = $stageDirectory."/correlationMatrix.eps";
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
	print $gnuPlot "set xlabel '\$".$massVariable."\$ [\$M_\\odot\$]'\n";
	print $gnuPlot "set ylabel '\$".$massVariable."\$ [\$M_\\odot\$]'\n";
	print $gnuPlot "set xrange [".$rangeLow.":".$rangeHigh."]\n";
	print $gnuPlot "set yrange [".$rangeLow.":".$rangeHigh."]\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set logscale cb\n";
	print $gnuPlot "set format cb '\$10^{\%L}\$'\n";
	print $gnuPlot "set pm3d map\n";
	print $gnuPlot "set pm3d explicit\n";
	print $gnuPlot "set pm3d corners2color c1\n";
	print $gnuPlot "set palette rgbformulae 21,22,23 negative\n";
	print $gnuPlot "splot '-' with pm3d notitle\n";
	for(my $i=0;$i<$correlation->dim(0);++$i) {
	    for(my $j=0;$j<$correlation->dim(1);++$j) {
 		print $gnuPlot $mass->(($i))." ".$mass->(($j))." ".abs($correlation->(($i),($j)))."\n";
	    }
	    print $gnuPlot "\n"
		unless ( $i == $correlation->dim(0)-1 );
	}
	print $gnuPlot "e\n";
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS);
    }

    # Launch MCMC simulation to find parameters which give a good fit to this mass function.
    unless ( -e $stageDirectory."/chains_0000.log" ) {
    	# Clean up any old files.
    	system("rm ".$stageDirectory."/chains_*.log");
	# Parse the basic covariance matrix parameter file.
	my $xmlP       = new XML::Simple;
	my $parameters = $xmlP->XMLin($parameterFile);
    	# Make stage-specific copies of the MCMC configuration file.
	my $xml        = new XML::Simple(KeyAttr => []);
	my $mcmcConfig = $xml->XMLin($mcmcConfigFile);
	@parameterNames = map {$_->{'name'}} @{$mcmcConfig->{'parameters'}->{'parameter'}};
	$mcmcConfig->{'simulation'}->{'logFileRoot'         } = $stageDirectory."/chains";
	# Locate the likelihood descriptor associated with the mass function likelihood function.
	my $likelihood;
	if ( $mcmcConfig->{'likelihood'}->{'type'} eq "gaussianRegression" ) {
	    $likelihood = $mcmcConfig->{'likelihood'}->{'simulatorLikelihood'};
	} else {
	    $likelihood = $mcmcConfig->{'likelihood'}                         ;
	}
	$likelihood->{'massFunctionFileName'} = $covarianceMatrixFile;
	if ( defined($redshiftIndex) ) {
	    $likelihood->{$parameters->{'parameter'}->{'surveyGeometryMethod'}->{'value'}.'RedshiftBin'} = $redshiftIndex;
	    $likelihood->{'redshiftMinimum'} = ${$config->{'redshift'}}[$redshiftIndex]->{'minimum'};
	    $likelihood->{'redshiftMaximum'} = ${$config->{'redshift'}}[$redshiftIndex]->{'maximum'};
	}
	my $xmlOutput = new XML::Simple (NoAttr => 1, KeyAttr => [], RootName=>"simulationConfig");
     	open(oHndl,">".$stageDirectory."/mcmcConfig.xml");
     	print oHndl $xmlOutput->XMLout($mcmcConfig);
     	close(oHndl);	
	# Determine number of threads to use.
	my $threadCount = $config->{'nodeCount'}*$config->{'threadsPerNode'};
    	# Submit the job.
    	my $cmfJob = 
    	{
    	    batchFile  => $stageDirectory."/constrainCMF.pbs",
    	    queue      => "batch",
    	    nodes      => "nodes=".$config->{'nodeCount'}.":ppn=".$config->{'threadsPerNode'},
    	    wallTime   => "300:00:00",
    	    outputFile => $stageDirectory."/constrainCMF.log",
    	    name       => $pbsLabel."Stage".$stage."CMF",
    	    commands   => "mpirun -np ".$threadCount." -hostfile \$PBS_NODEFILE Constrain_Galacticus.exe ".$stageDirectory."/".$parameterFileLeaf." ".$stageDirectory."/mcmcConfig.xml"
    	};
    	&Submit_To_PBS($cmfJob);
    	# Since the posterior has been updated, remove the best fit mass function file.
    	unlink($stageDirectory."/massFunctionBestFit.hdf5");
    	# Since the posterior has been updated, remove the triangle plot files.
    	system("rm -f ".$stageDirectory."/triangle*");
    }

    # Generate a triangle plot.
    unless ( -e $stageDirectory."/triangle.tex" ) {
    	# Construct the command to generate the triangle plot.
    	my $command;
    	$command .= "constraints/visualization/mcmcVisualizeTriangle.pl";
    	$command .= " ".$stageDirectory."/chains ";
	$command .= " ".$mcmcConfigFile;
	$command .= " --workDirectory ".$stageDirectory;
    	$command .= " --scale 0.1";
    	$command .= " --ngood 100000";
    	$command .= " --ngrid 100";
    	$command .= " --output ".$stageDirectory."/triangle";
    	$command .= " --property 'alphaSatellite:linear:xLabel=\$\\alpha\$:zLabel=\${\\rm d}p/{\\rm d}\\alpha\$'";
    	$command .= " --property 'log10M1:linear:xLabel=\$\\log_{10}(M_1/M_\\odot)\$:zLabel=\${\\rm d}p/{\\rm d}\\log_{10} M_1\$'";
    	$command .= " --property 'log10Mstar0:linear:xLabel=\$\\log_{10}(M_{\\star,0}/M_\\odot)\$:zLabel=\${\\rm d}p/{\\rm d}\\log_{10} M_{\\star,0}\$'";
    	$command .= " --property 'beta:linear:xLabel=\$\\beta\$:zLabel=\${\\rm d}p/{\\rm d}\\beta\$'";
    	$command .= " --property 'delta:linear:xLabel=\$\\delta\$:zLabel=\${\\rm d}p/{\\rm d}\\delta\$'";
    	$command .= " --property 'gamma:linear:xLabel=\$\\gamma\$:zLabel=\${\\rm d}p/{\\rm d}\\gamma\$'";
    	$command .= " --property 'sigmaLogMstar:logarithmic:xLabel=\$\\sigma_{\\log ".$massVariable."}\$:zLabel=\${\\rm d}p/{\\rm d}\\log \\sigma_{\\log ".$massVariable."}\$'";
    	$command .= " --property 'BCut:linear:xLabel=\$B_{\\rm cut}\$:zLabel=\${\\rm d}p/{\\rm d}B_{\\rm cut}\$'";
    	$command .= " --property 'BSatellite:linear:xLabel=\$B_{\\rm sat}\$:zLabel=\${\\rm d}p/{\\rm d}B_{\\rm sat}\$'";
    	$command .= " --property 'betaCut:linear:xLabel=\$\\beta_{\\rm cut}\$:zLabel=\${\\rm d}p/{\\rm d}\\beta_{\\rm cut}\$'";
    	$command .= " --property 'betaSatellite:linear:xLabel=\$\\beta_{\\rm sat}\$:zLabel=\${\\rm d}p/{\\rm d}\\beta_{\\rm sat}\$'";
    	# Submit the job.
    	my $triangleJob = 
    	{
    	    batchFile  => $stageDirectory."/triangle.pbs",
    	    queue      => "batch",
    	    nodes      => "nodes=1:ppn=1",
    	    wallTime   => "20:00:00",
    	    outputFile => $stageDirectory."/triangle.log",
    	    name       => $pbsLabel."Stage".$stage."Triangle",
    	    commands   => "mpirun -np 1 -hostfile \$PBS_NODEFILE ".$command
    	};
    	&Submit_To_PBS($triangleJob);
    }

    # Find the maximum likelihood parameters.
    unless ( -e $stageDirectory."/massFunctionBestFit.xml" ) {
	my $maximumLikelihood = -1.0e30;
	opendir(my $stageHndl,$stageDirectory);
	while ( my $file = readdir($stageHndl) ) {
	    if( $file =~ m/chains_\d+\.log/ ) {
		open(iHndl,$stageDirectory."/".$file);
		while ( my $line = <iHndl> ) {
		    $line =~ s/^\s*//;
		    $line =~ s/\s*$//;
		    my @columns = split(/\s+/,$line);
		    if ( $columns[4] > $maximumLikelihood ) {
			$maximumLikelihood = $columns[4];
			@bestFit = @columns[5..$#columns];
		    }
		}
		close(iHndl);
	    }
	}
	closedir($stageHndl);
    	# Write the results to file.
    	my $parameters;
	my $xml        = new XML::Simple(KeyAttr => []);
	my $mcmcConfig = $xml->XMLin($mcmcConfigFile);
    	for(my $i=0;$i<scalar(@parameterNames);++$i) {
	    my $value = $bestFit[$i];
	    $value = exp($value)
		if ( ${$mcmcConfig->{'parameters'}->{'parameter'}}[$i]->{'mapping'}->{'type'} eq "logarithmic" );
    	    push(
    		@{$parameters->{'parameter'}},
    		{
    		    name  => $parameterNames[$i],
    		    value => $bestFit       [$i]
    		}
    		);
    	}
    	# Output the best fit parameters.
    	my $xmlOutput = new XML::Simple(NoAttr=>1, RootName=>"parameters");
    	open(oHndl,">".$stageDirectory."/massFunctionBestFit.xml");
    	print oHndl $xmlOutput->XMLout($parameters);
    	close(oHndl);
    } else {
    	# Read in the best-fit parameters.
    	my $xmlIn      = new XML::Simple;
    	my $parameters = $xmlIn->XMLin($stageDirectory."/massFunctionBestFit.xml", KeyAttr => []);
    	# Extract to arrays.
    	foreach my $parameter ( @{$parameters->{'parameter'}} ) {
    	    push(@bestFit ,$parameter->{'value'});
    	}
    }

    # Generate the maximum likelihood mass function.
    unless ( -e $stageDirectory."/massFunctionBestFit.hdf5" ) {
	# The maximum likelihood mass function has changed, so remove the covariance matrix
	# for the next stage, forcing it to be updated.
	my $stageNext = $stage+1;
	(my $covarianceMatrixFileNext = $stageDirectory."/".$constraintFileLeaf) =~ s/stage$stage/stage$stageNext/;
	system("rm -f ".$covarianceMatrixFileNext)
	    if ( -e $covarianceMatrixFileNext );	
	# Generate the maximum likelihood mass function.
	my $xml        = new XML::Simple;
	my $parameters = $xml->XMLin($parameterFile);
	for(my $i=0;$i<scalar(@parameterNames);++$i) {
	    my $name = "conditionalMassFunctionBehroozi".ucfirst($parameterNames[$i]);
	    $parameters->{'parameter'}->{$name}->{'value'} = $bestFit[$i];
	}
	# Modify redshift ranges.
	if ( defined($redshiftIndex) ) {
	    $parameters->{'parameter'}->{'conditionalMassFunctionRedshiftMinimum'}->{'value'} = ${$config->{'redshift'}}[$redshiftIndex]->{'minimum'};
	    $parameters->{'parameter'}->{'conditionalMassFunctionRedshiftMaximum'}->{'value'} = ${$config->{'redshift'}}[$redshiftIndex]->{'maximum'};
	    $parameters->{'parameter'}->{$parameters->{'parameter'}->{'surveyGeometryMethod'}->{'value'}.'RedshiftBin'}->{'value'} = $redshiftIndex;
	}
	# Set file name.
	$parameters->{'parameter'}->{'conditionalMassFunctionOutputFileName'}->{'value'} = $stageDirectory."/massFunctionBestFit.hdf5";
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
	open(oHndl,">".$stageDirectory."/massFunctionGenerate.xml");
	print oHndl $xmlOutput->XMLout($parameters);
	close(oHndl);
	system($galacticusPath."Conditional_Mass_Function.exe ".$stageDirectory."/massFunctionGenerate.xml");
	# Read the best fit mass function data.
	my $bestFit      = new PDL::IO::HDF5($stageDirectory."/massFunctionBestFit.hdf5");
	my $mass         = $bestFit->dataset('mass'        )->get();
	my $massFunction = $bestFit->dataset('massFunction')->get();
	# Read the observational data.
	my $observed             = new PDL::IO::HDF5($stageDirectory."/".$constraintFileLeaf);
	my $massObserved         = $observed->dataset('mass'                )->get();
	my $massFunctionObserved = $observed->dataset('massFunctionObserved')->get();
	my $covariance           = $observed->dataset('covariance'          )->get();
	my $errorObserved        = sqrt($covariance->diagonal(0,1));
	# Create a plot of this.
	my $plot;
	my $gnuPlot;
	my $plotFile = $stageDirectory."/massFunctionBestFit.pdf";
	(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
	print $gnuPlot "set title 'Best fit mass function for Stage ".$stage."'\n";
	print $gnuPlot "set xlabel '".ucfirst($massType)." mass; \$".$massVariable."\,[M_\\odot]\$'\n";
	print $gnuPlot "set ylabel 'Mass function; \$\\phi(M)\$ [Mpc\$^{-3} \\log(".$massVariable.")\$]'\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.275,0.16\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	my $xMinimum = 0.5*minimum($massObserved        );
	my $xMaximum = 2.0*maximum($massObserved        );
	my $yMinimum = 0.5*minimum($massFunctionObserved);
	my $yMaximum = 2.0*maximum($massFunctionObserved);
	print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	print $gnuPlot "set pointsize 2.0\n";
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $massObserved,
	    $massFunctionObserved,
	    errorUp    => $errorObserved,
	    errorDown  => $errorObserved,
	    style      => "point",
	    weight     => [5,3],
	    symbol     => [6,7],
	    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
	    title      => $sourceLabel
	    );
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $mass,
	    $massFunction,
	    style      => "point",
	    weight     => [5,3],
	    symbol     => [6,7],
	    pointSize  => 0.5,
	    color      => $PrettyPlots::colorPairs{'redYellow'},
	    title      => "Maximum likelihood fit"
	    );
	&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS);
    }
}

exit;

sub Submit_To_PBS {
    # Submit a job to the PBS queue and wait for it to complete.
    my $jobDescriptor = shift;
    # Assert that job commands must be present.
    die("Submit_To_PBS: no job commands were supplied")
	unless ( exists($jobDescriptor->{'commands'}) );
    # Assert that batch file must be present.
    die("Submit_To_PBS: no batchFile was supplied")
	unless ( exists($jobDescriptor->{'batchFile'}) );
    # Generate the bacth file.
    open(oHndl,">".$jobDescriptor->{'batchFile'});
    print oHndl "#!/bin/bash\n";
    print oHndl "#PBS -q ".$jobDescriptor->{'queue'}."\n"
	if ( exists($jobDescriptor->{'queue'}) );
    print oHndl "#PBS -l ".$jobDescriptor->{'nodes'}."\n"
	if ( exists($jobDescriptor->{'nodes'}) );
    print oHndl "#PBS -j oe\n";
    print oHndl "#PBS -o ".$jobDescriptor->{'outputFile'}."\n"
	if ( exists($jobDescriptor->{'outputFile'}) );
    print oHndl "#PBS -l walltime=".$jobDescriptor->{'wallTime'}."\n"
	if ( exists($jobDescriptor->{'wallTime'}) );
    print oHndl "#PBS -N ".$jobDescriptor->{'name'}."\n"
	if ( exists($jobDescriptor->{'name'}) );
    print oHndl "#PBS -V\n";
    print oHndl "cd \$PBS_O_WORKDIR\n";
    print oHndl "export LD_LIBRARY_PATH=\$HOME/Galacticus/OpenMPI/lib:\$HOME/Galacticus/Tools/lib:\$HOME/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
    print oHndl "export PATH=\$HOME/Galacticus/OpenMPI/bin:\$HOME/Galacticus/Tools/bin:\$HOME/perl5/bin:\$PATH\n";
    print oHndl "export GFORTRAN_ERROR_DUMPCORE=YES\n";
    print oHndl "export PERL_LOCAL_LIB_ROOT=\"\$HOME/perl5\"\n";
    print oHndl "export PERL_MB_OPT=\"--install_base \$HOME/perl5\"\n";
    print oHndl "export PERL_MM_OPT=\"INSTALL_BASE=\$HOME/perl5\"\n";
    print oHndl "export PERL5LIB=\"\$HOME/perl5/lib/perl5/x86_64-linux-thread-multi:\$HOME/perl5/lib/perl5\"\n";
    print oHndl "export PYTHONPATH=\"\$HOME/Galacticus/Tools/lib/python:\$HOME/Galacticus/Tools/lib/python2.7:/share/apps/atipa/acms/lib\"\n";
    print oHndl "ulimit -t unlimited\n";
    print oHndl "ulimit -c unlimited\n";
    print oHndl $jobDescriptor->{'commands'};
    close(oHndl);
    open(pHndl,"qsub ".$jobDescriptor->{'batchFile'}." |");
    my $jobID = "";
    while ( my $line = <pHndl> ) {
	if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
    }
    close(pHndl);
    print "PBS job ";
    print "\"".$jobDescriptor->{'name'}."\" "
	if ( exists($jobDescriptor->{'name'}) );
    print "submitted as ".$jobID."\n";
    # Wait for the job to finish.
    print "Waiting for PBS job to finish....\n";
    my $jobIsRunning = 1;
    while ( $jobIsRunning == 1 ) {
	$jobIsRunning = 0;
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {
		$jobIsRunning = 1
		    if ( $1 eq $jobID );
	    }
	}
	close(pHndl);
	sleep 10;
    }

}
