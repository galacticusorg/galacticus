#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
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
die("Usage: covarianceMatrix.pl <configFile>")
    unless ( scalar(@ARGV) == 1 );
my $configFile = $ARGV[0];
my $xml        = new XML::Simple;
my $config     = $xml->XMLin($configFile);

# Validate the config file.
die("covarianceMatrix.pl: config file must specify baseDirectory" )
    unless ( exists($config->{'baseDirectory' }) );
die("covarianceMatrix.pl: config file must specify parameterFile" )
    unless ( exists($config->{'parameterFile' }) );
die("covarianceMatrix.pl: config file must specify constraintFile")
    unless ( exists($config->{'constraintFile'}) );
die("covarianceMatrix.pl: config file must specify bieFile"       )
    unless ( exists($config->{'bieFile'       }) );

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
my $bieFile        = $config->{'bieFile'       };
(my $parameterFileLeaf  = $parameterFile ) =~ s/^.*\/([^\/]+)$/$1/;
(my $parameterFilePath  = $parameterFile ) =~ s/^(.*\/).*$/$1/;
(my $constraintFileLeaf = $constraintFile) =~ s/^.*\/([^\/]+)$/$1/;
(my $constraintFilePath = $constraintFile) =~ s/^(.*\/).*$/$1/;
(my $bieFileLeaf        = $bieFile       ) =~ s/^.*\/([^\/]+)$/$1/;
(my $bieFilePath        = $bieFile       ) =~ s/^(.*\/).*$/$1/;

# Ensure that the covariance matrix and mass function codes are built.
system("make Mass_Function_Covariance.exe Conditional_Mass_Function.exe");

# Lists that will be used to store the best fit parameter names and values from one stage to the next.
my @parameterNames;
my @bestFit;

# Perform a specified number of stages of the algorithm.
my $stageCount = 4;
$stageCount = $config->{'stageCount'}
   if ( exists($config->{'stageCount'}) );
for(my $stage=0;$stage<=$stageCount;++$stage) {
    # Create an output folder for this stage.
    system("mkdir -p ".$baseDirectory."stage".$stage);
    # Generate the covariance matrix file.
    my $covarianceMatrixFile = $baseDirectory."stage".$stage."/".$constraintFileLeaf;
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
	# For stages after the first, insert the best fit parameters from the previous stage.
	if ( $stage > 0 ) {
	    for(my $i=0;$i<scalar(@parameterNames);++$i) {
		$parameters->{'parameter'}->{$parameterNames[$i]}->{'value'} = $bestFit[$i];
	    }
	}
	# Write out the new parameter file.
	my $covarianceMatrixParameterFile = $baseDirectory."stage".$stage."/".$parameterFileLeaf;
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
	open(oHndl,">".$covarianceMatrixParameterFile);
	print oHndl $xmlOutput->XMLout($parameters);
	close(oHndl);
	# Generate the covariance matrix.
	my $matrixJob = 
	{
	    batchFile  => $baseDirectory."stage".$stage."/generateMatrix.pbs",
	    queue      => "batch",
	    nodes      => "nodes=1:ppn=12",
	    wallTime   => "200:00:00",
	    outputFile => $baseDirectory."stage".$stage."/generateMatrix.log",
	    name       => $pbsLabel."Stage".$stage."MatrixGen",
	    commands   => "mpirun -np 1 -hostfile \$PBS_NODEFILE ".$galacticusPath."constraints/dataAnalysis/scripts/generateCovarianceMatrix.pl ".$parameterFilePath."stage".$stage."/".$parameterFileLeaf." ".$configFile
	};
	&Submit_To_PBS($matrixJob);
	# Since we have a new covariance matrix file, remove the BIE statelog to force it to be regenerated.
	unlink($baseDirectory."stage".$stage."/bie.statelog");
	# Also remove any plot of the correlation matrix.
	unlink($baseDirectory."stage".$stage."/correlationMatrix.pdf");
	# For stages past the 0th, determine how well converged the matrix is.
	if ( $stage > 0 ) {
	    my $oldCovarianceMatrixFile = $constraintFilePath."stage".($stage-1)."/".$constraintFileLeaf;
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
    unless ( -e $baseDirectory."stage".$stage."/correlationMatrix.pdf" ) {
	# Read the covariance matrix.
	my $hdfFile      = new PDL::IO::HDF5($covarianceMatrixFile);
	my $correlation  = $hdfFile->dataset("correlation")->get();
	my $mass         = $hdfFile->dataset("mass"       )->get();

	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	my $rangeLow  = $mass->(( 0));
	my $rangeHigh = $mass->((-1));
	$plotFileEPS = $baseDirectory."stage".$stage."/correlationMatrix.eps";
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

    # Launch BIE to find parameters which give a good fit to this mass function.
    unless ( -e $baseDirectory."stage".$stage."/bie.statelog" ) {
    	# Clean up any old files.
    	system("rm ".$baseDirectory."stage".$stage."/bie*");
    	# Make stage-specific copies of the bie and control scripts.
    	&Make_Stage_Specific
    	    (
    	     $stage,
    	     $bieFile,
    	     $baseDirectory."stage".$stage."/covarianceMatrix.bie"
    	    );
	
        my $configCopy = clone($config);
     	$configCopy->{'bieFile'       } = $baseDirectory."stage".$stage."/covarianceMatrix.bie";
     	$configCopy->{'constraintFile'} = $baseDirectory."stage".$stage."/".$constraintFileLeaf;
     	$configCopy->{'parameterFile' } = $baseDirectory."stage".$stage."/".$parameterFileLeaf;
     	$configCopy->{'workDirectory' } = $baseDirectory."stage".$stage."/";
     	my $xml = new XML::Simple(RootName => 'constrain');
     	open(oHndl,">".$baseDirectory."stage".$stage."/covarianceMatrixControl.xml");
     	print oHndl $xml->XMLout($configCopy);
     	close(oHndl);	
    	# Submit the job.
    	my $cmfJob = 
    	{
    	    batchFile  => $baseDirectory."stage".$stage."/constrainCMF.pbs",
    	    queue      => "batch",
    	    nodes      => "nodes=4:ppn=12",
    	    wallTime   => "300:00:00",
    	    outputFile => $baseDirectory."stage".$stage."/constrainCMF.log",
    	    name       => $pbsLabel."Stage".$stage."CMF",
    	    commands   => "rm -rf pdir/covarianceMatrixStage".$stage."\nmpirun -np 48 -hostfile \$PBS_NODEFILE bie -f ".$baseDirectory."stage".$stage."/covarianceMatrix.bie"
    	};
    	&Submit_To_PBS($cmfJob);
    	# Since the posterior has been updated, remove the best fit mass function file.
    	unlink($baseDirectory."stage".$stage."/massFunctionBestFit.hdf5");
    	# Since the posterior has been updated, remove the triangle plot files.
    	system("rm -f ".$baseDirectory."stage".$stage."/triangle*");
    }

    # Generate a triangle plot.
    unless ( -e $baseDirectory."stage".$stage."/triangle.tex" ) {
    	# Construct the command to generate the triangle plot.
    	my $command;
    	$command .= "constraints/visualization/bieStatelogVisualizeTriangle.pl";
    	$command .= " ".$baseDirectory."stage".$stage."/bie.statelog ";
    	$command .= " --scale 0.1";
    	$command .= " --ngood 500000";
    	$command .= " --ngrid 100";
    	$command .= " --output ".$baseDirectory."stage".$stage."/triangle";
    	$command .= " --property 'conditionalMassFunctionBehrooziAlphaSatellite:linear:xLabel=\$\\alpha\$:zLabel=\${\\rm d}p/{\\rm d}\\alpha\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziLog10M1:linear:xLabel=\$\\log_{10}(M_1/M_\\odot)\$:zLabel=\${\\rm d}p/{\\rm d}\\log_{10} M_1\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziLog10Mstar0:linear:xLabel=\$\\log_{10}(M_{\\star,0}/M_\\odot)\$:zLabel=\${\\rm d}p/{\rm d}\\log_{10} M_{\\star,0}\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziBeta:linear:xLabel=\$\\beta\$:zLabel=\${\\rm d}p/{\\rm d}\\beta\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziDelta:linear:xLabel=\$\\delta\$:zLabel=\${\\rm d}p/{\\rm d}\\delta\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziGamma:linear:xLabel=\$\\gamma\$:zLabel=\${\\rm d}p/{\\rm d}\\gamma\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziSigmaLogMstar:linear:xLabel=\$\\sigma_{\\log ".$massVariable."}\$:zLabel=\${\\rm d}p/{\\rm d}\\sigma_{\\log ".$massVariable."}\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziBCut:linear:xLabel=\$B_{\\rm cut}\$:zLabel=\${\\rm d}p/{\\rm d}B_{\\rm cut}\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziBSatellite:linear:xLabel=\$B_{\\rm sat}\$:zLabel=\${\\rm d}p/{\\rm d}B_{\\rm sat}\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziBetaCut:linear:xLabel=\$\\beta_{\\rm cut}\$:zLabel=\${\\rm d}p/{\\rm d}\\beta_{\\rm cut}\$'";
    	$command .= " --property 'conditionalMassFunctionBehrooziBetaSatellite:linear:xLabel=\$\\beta_{\\rm sat}\$:zLabel=\${\\rm d}p/{\\rm d}\\beta_{\\rm sat}\$'";
    	# Submit the job.
    	my $triangleJob = 
    	{
    	    batchFile  => $baseDirectory."stage".$stage."/triangle.pbs",
    	    queue      => "batch",
    	    nodes      => "nodes=1:ppn=1",
    	    wallTime   => "20:00:00",
    	    outputFile => $baseDirectory."stage".$stage."/triangle.log",
    	    name       => $pbsLabel."Stage".$stage."Triangle",
    	    commands   => "mpirun -np 1 -hostfile \$PBS_NODEFILE ".$command
    	};
    	&Submit_To_PBS($triangleJob);
    }

    # Find the maximum likelihood parameters.
    unless ( -e $baseDirectory."stage".$stage."/massFunctionBestFit.xml" ) {
    	open(iHndl,$baseDirectory."stage".$stage."/bie.statelog");
    	my $maximumLikelihood = -1.0e30;
    	while ( my $line = <iHndl> ) {
    	    $line =~ s/^\s+//;
    	    if ( $line =~ m/^\"/ ) {
    		$line =~ s/\"/ /g;
    		$line =~ s/^\s+//;
    		my @columns = split(/\s+/,$line);
    		@parameterNames = @columns[5..$#columns];
    	    } else {
    		my @columns = split(/\s+/,$line);
    		if ( $columns[2] > $maximumLikelihood ) {
    		    $maximumLikelihood = $columns[2];
    		    @bestFit = @columns[5..$#columns];
    		}
    	    }
    	}
    	close(iHndl);
    	# Write the results to file.
    	my $parameters;
    	for(my $i=0;$i<scalar(@parameterNames);++$i) {
    	    push(
    		@{$parameters->{'parameter'}},
    		{
    		    name  => $parameterNames[$i],
    		    value => $bestFit[$i]
    		}
    		);
    	}
    	# Output the best fit parameters.
    	my $xmlOutput = new XML::Simple(NoAttr=>1, RootName=>"parameters");
    	open(oHndl,">".$baseDirectory."stage".$stage."/massFunctionBestFit.xml");
    	print oHndl $xmlOutput->XMLout($parameters);
    	close(oHndl);
    } else {
    	# Read in the best-fit parameters.
    	my $xmlIn      = new XML::Simple;
    	my $parameters = $xmlIn->XMLin($baseDirectory."stage".$stage."/massFunctionBestFit.xml", KeyAttr => []);
    	# Extract to arrays.
    	foreach my $parameter ( @{$parameters->{'parameter'}} ) {
    	    push(@parameterNames,$parameter->{'name' });
    	    push(@bestFit       ,$parameter->{'value'});
    	}
    }

    # Generate the maximum likelihood mass function.
    unless ( -e $baseDirectory."stage".$stage."/massFunctionBestFit.hdf5" ) {
	# The maximum likelihood mass function has changed, so remove the covariance matrix
	# for the next stage, forcing it to be updated.
	my $stageNext = $stage+1;
	my $covarianceMatrixFileNext = $baseDirectory."stage".$stageNext."/".$constraintFileLeaf;
	system("rm -f ".$covarianceMatrixFileNext)
	    if ( -e $covarianceMatrixFileNext );	
	# Generate the maximum likelihood mass function.
	system($galacticusPath."constraints/dataAnalysis/scripts/computeLikelihood.pl ".$baseDirectory."stage".$stage."/covarianceMatrixControl.xml 0 ".join(" ",@bestFit));
	system("mv ".$baseDirectory."stage".$stage."/massFunction_0000.hdf5 ".$baseDirectory."stage".$stage."/massFunctionBestFit.hdf5");
	# Read the best fit mass function data.
	my $bestFit      = new PDL::IO::HDF5($baseDirectory."stage".$stage."/massFunctionBestFit.hdf5");
	my $mass         = $bestFit->dataset('mass'        )->get();
	my $massFunction = $bestFit->dataset('massFunction')->get();
	# Read the observational data.
	my $observed             = new PDL::IO::HDF5($baseDirectory."stage".$stage."/".$constraintFileLeaf);
	my $massObserved         = $observed->dataset('mass'                )->get();
	my $massFunctionObserved = $observed->dataset('massFunctionObserved')->get();
	my $covariance           = $observed->dataset('covariance'          )->get();
	my $errorObserved        = sqrt($covariance->diagonal(0,1));
	# Create a plot of this.
	my $plot;
	my $gnuPlot;
	my $plotFile = $baseDirectory."stage".$stage."/massFunctionBestFit.pdf";
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
    print oHndl "export LD_LIBRARY_PATH=\$HOME/Galacticus/Tools/lib:\$HOME/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
    print oHndl "export PATH=\$HOME/Galacticus/Tools/bin:\$HOME/perl5/bin:\$PATH\n";
    print oHndl "export GFORTRAN_ERROR_DUMPCORE=YES\n";
    print oHndl "export PERL_LOCAL_LIB_ROOT=\"\$HOME/perl5\"\n";
    print oHndl "export PERL_MB_OPT=\"--install_base \$HOME/perl5\"\n";
    print oHndl "export PERL_MM_OPT=\"INSTALL_BASE=\$HOME/perl5\"\n";
    print oHndl "export PERL5LIB=\"\$HOME/perl5/lib/perl5/x86_64-linux-thread-multi:\$HOME/perl5/lib/perl5\"\n";
    print oHndl "export PYTHONPATH=\"\$HOME/Galacticus/Tools/py-lib\"\n";
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
	sleep 60;
    }

}

sub Make_Stage_Specific {
    # Make a copy of a file that is specific to the given stage.
    my $stage = shift;
    my $from  = shift;
    my $to    = shift;
    # Copy the file, inserting the stage number where necessary.
    open(iHndl,$from);
    open(oHndl,">".$to);
    while ( my $line = <iHndl> ) {
	$line =~ s/\%STAGE/$stage/g;
	print oHndl $line;
    }
    close(iHndl);
    close(oHndl);
}
