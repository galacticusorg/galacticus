#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
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
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Find the maximum likelihood estimate of the covariance matrix for the Li & White (2009) SDSS stellar mass function.
# Andrew Benson (05-July-2012)

# Ensure that the covariance matrix and mass function codesare built.
system("make Mass_Function_Covariance.exe Conditional_Stellar_Mass_Function.exe");

# Specify the work directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/";

# Lists that will be used to store the best fit parameter names and values from one stage to the next.
my @parameterNames;
my @bestFit;

# Perform a specified number of stages of the algorithm.
my $stageCount = 4;
for(my $stage=0;$stage<=$stageCount;++$stage) {
    # Create an output folder for this stage.
    system("mkdir -p ".$workDirectory."stage".$stage);
    # Generate the covariance matrix file.
    my $covarianceMatrixFile = $workDirectory."stage".$stage."/covarianceMatrix.hdf5";
    unless ( -e $covarianceMatrixFile ) {
	# Parse the basic covariance matrix parameter file.
	my $xml        = new XML::Simple;
	my $parameters = $xml->XMLin($workDirectory."covarianceMatrix.xml");
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
	my $covarianceMatrixParameterFile = $workDirectory."stage".$stage."/covarianceMatrix.xml";
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
	open(oHndl,">".$covarianceMatrixParameterFile);
	print oHndl $xmlOutput->XMLout($parameters);
	close(oHndl);
	# Generate the covariance matrix.
	my $matrixJob = 
	{
	    batchFile  => $workDirectory."stage".$stage."/generateMatrix.pbs",
	    queue      => "batch",
	    nodes      => "nodes=1:ppn=12",
	    wallTime   => "200:00:00",
	    outputFile => $workDirectory."stage".$stage."/generateMatrix.log",
	    name       => "stage".$stage."MatrixGen",
	    commands   => "mpirun -np 1 -hostfile \$PBS_NODEFILE ".$workDirectory."generateCovarianceMatrix.pl ".$workDirectory."stage".$stage."/covarianceMatrix.xml"
	};
	&Submit_To_PBS($matrixJob);
	# Since we have a new covariance matrix file, remove the BIE statelog to force it to be regenerated.
	unlink($workDirectory."stage".$stage."/bie.statelog");
	# Also remove any plot of the correlation matrix.
	unlink($workDirectory."stage".$stage."/correlationMatrix.pdf");
	# For stages past the 0th, determine how well converged the matrix is.
	if ( $stage > 0 ) {
	    my $oldCovarianceMatrixFile = $workDirectory."stage".($stage-1)."/covarianceMatrix.hdf5";
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
    unless ( -e $workDirectory."stage".$stage."/correlationMatrix.pdf" ) {
	# Read the covariance matrix.
	my $hdfFile      = new PDL::IO::HDF5($covarianceMatrixFile);
	my $correlation  = $hdfFile->dataset("correlation")->get();
	my $mass         = $hdfFile->dataset("mass"       )->get();

	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	my $rangeLow  = $mass->(( 0));
	my $rangeHigh = $mass->((-1));
	$plotFileEPS = $workDirectory."stage".$stage."/correlationMatrix.eps";
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
	print $gnuPlot "set xlabel '\$M_\\star\$ [\$M_\\odot\$]'\n";
	print $gnuPlot "set ylabel '\$M_\\star\$ [\$M_\\odot\$]'\n";
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
 		print $gnuPlot $mass->(($i))." ".$mass->(($j))." ".$correlation->(($i),($j))."\n";
	    }
	    print $gnuPlot "\n"
		unless ( $i == $correlation->dim(0)-1 );
	}
	print $gnuPlot "e\n";
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS);
    }

    # Launch BIE to find parameters which give a good fit to this mass function.
    unless ( -e $workDirectory."stage".$stage."/bie.statelog" ) {
	# Clean up any old files.
	system("rm ".$workDirectory."stage".$stage."/bie*");
	# Make stage-specific copies of the bie and control scripts.
	&Make_Stage_Specific
	    (
	     $stage,
	     $workDirectory."covarianceMatrix.bie",
	     $workDirectory."stage".$stage."/covarianceMatrix.bie"
	    );
	&Make_Stage_Specific
	    (
	     $stage,
	     $workDirectory."covarianceMatrixControl.xml",
	     $workDirectory."stage".$stage."/covarianceMatrixControl.xml"
	    );
	# Submit the job.
	my $csmfJob = 
	{
	    batchFile  => $workDirectory."stage".$stage."/constrainCSMF.pbs",
	    queue      => "batch",
	    nodes      => "nodes=4:ppn=12",
	    wallTime   => "300:00:00",
	    outputFile => $workDirectory."stage".$stage."/constrainCSMF.log",
	    name       => "stage".$stage."CSMF",
	    commands   => "rm -rf pdir/covarianceMatrixStage".$stage."\nmpirun -np 48 -hostfile \$PBS_NODEFILE bie -f ".$workDirectory."stage".$stage."/covarianceMatrix.bie"
	};
	&Submit_To_PBS($csmfJob);
	# Since the posterior has been updated, remove the best fit mass function file.
	unlink($workDirectory."stage".$stage."/massFunctionBestFit.hdf5");
	# Since the posterior has been updated, remove the triangle plot files.
	system("rm -f ".$workDirectory."stage".$stage."/triangle*");
    }

    # Generate a triangle plot.
    unless ( -e $workDirectory."stage".$stage."/triangle.tex" ) {
	# Construct the command to generate the triangle plot.
	my $command;
	$command .= "constraints/visualization/bieStatelogVisualizeTriangle.pl";
	$command .= " ".$workDirectory."stage".$stage."/bie.statelog ";
	$command .= " --scale 0.1";
	$command .= " --ngood 500000";
	$command .= " --ngrid 100";
	$command .= " --output ".$workDirectory."stage".$stage."/triangle";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziAlphaSatellite:linear:xLabel=\$\\alpha\$:zLabel=\${\\rm d}p/{\\rm d}\\alpha\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziLog10M1:linear:xLabel=\$\\log_{10}(M_1/M_\\odot)\$:zLabel=\${\\rm d}p/{\\rm d}\\log_{10} M_1\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziLog10Mstar0:linear:xLabel=\$\\log_{10}(M_{\\star,0}/M_\\odot)\$:zLabel=\${\\rm d}p/{\rm d}\\log_{10} M_{\\star,0}\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziBeta:linear:xLabel=\$\\beta\$:zLabel=\${\\rm d}p/{\\rm d}\\beta\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziDelta:linear:xLabel=\$\\delta\$:zLabel=\${\\rm d}p/{\\rm d}\\delta\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziGamma:linear:xLabel=\$\\gamma\$:zLabel=\${\\rm d}p/{\\rm d}\\gamma\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziSigmaLogMstar:linear:xLabel=\$\\sigma_{\\log M_\\star}\$:zLabel=\${\\rm d}p/{\\rm d}\\sigma_{\\log M_\\star}\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziBCut:linear:xLabel=\$B_{\\rm cut}\$:zLabel=\${\\rm d}p/{\\rm d}B_{\\rm cut}\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziBSatellite:linear:xLabel=\$B_{\\rm sat}\$:zLabel=\${\\rm d}p/{\\rm d}B_{\\rm sat}\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziBetaCut:linear:xLabel=\$\\beta_{\\rm cut}\$:zLabel=\${\\rm d}p/{\\rm d}\\beta_{\\rm cut}\$'";
	$command .= " --property 'conditionalStellarMassFunctionBehrooziBetaSatellite:linear:xLabel=\$\\beta_{\\rm sat}\$:zLabel=\${\\rm d}p/{\\rm d}\\beta_{\\rm sat}\$'";
	# Submit the job.
	my $triangleJob = 
	{
	    batchFile  => $workDirectory."stage".$stage."/triangle.pbs",
	    queue      => "batch",
	    nodes      => "nodes=1:ppn=1",
	    wallTime   => "20:00:00",
	    outputFile => $workDirectory."stage".$stage."/triangle.log",
	    name       => "stage".$stage."Triangle",
	    commands   => "mpirun -np 1 -hostfile \$PBS_NODEFILE ".$command
	};
	&Submit_To_PBS($triangleJob);
    }

    # Find the maximum likelihood parameters.
    unless ( -e $workDirectory."stage".$stage."/massFunctionBestFit.xml" ) {
	open(iHndl,$workDirectory."stage".$stage."/bie.statelog");
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
	open(oHndl,">".$workDirectory."stage".$stage."/massFunctionBestFit.xml");
	print oHndl $xmlOutput->XMLout($parameters);
	close(oHndl);
    } else {
	# Read in the best-fit parameters.
	my $xmlIn      = new XML::Simple;
	my $parameters = $xmlIn->XMLin($workDirectory."stage".$stage."/massFunctionBestFit.xml", KeyAttr => []);
	# Extract to arrays.
	foreach my $parameter ( @{$parameters->{'parameter'}} ) {
	    push(@parameterNames,$parameter->{'name' });
	    push(@bestFit       ,$parameter->{'value'});
	}
    }

    # Generate the maximum likelihood mass function.
    unless ( -e $workDirectory."stage".$stage."/massFunctionBestFit.hdf5" ) {
	# The maximum likelihood mass function has changed, so remove the covariance matrix
	# for the next stage, forcing it to be updated.
	my $stageNext = $stage+1;
	my $covarianceMatrixFileNext = $workDirectory."stage".$stageNext."/covarianceMatrix.hdf5";
	system("rm -f ".$covarianceMatrixFileNext)
	    if ( -e $covarianceMatrixFileNext );	
	# Generate the maximum likelihood mass function.
	system($workDirectory."computeLikelihood.pl ".$workDirectory."stage".$stage."/covarianceMatrixControl.xml 0 ".join(" ",@bestFit));
	system("mv ".$workDirectory."stage".$stage."/massFunction_0000.hdf5 ".$workDirectory."stage".$stage."/massFunctionBestFit.hdf5");
	# Read the best fit mass function data.
	my $bestFit      = new PDL::IO::HDF5($workDirectory."stage".$stage."/massFunctionBestFit.hdf5");
	my $stellarMass  = $bestFit->dataset('stellarMass' )->get();
	my $massFunction = $bestFit->dataset('massFunction')->get();
	# Read the observational data.
	my $observed     = new PDL::IO::HDF5($workDirectory."stage".$stage."/covarianceMatrix.hdf5");
	my $stellarMassObserved  = $observed->dataset('mass'                )->get();
	my $massFunctionObserved = $observed->dataset('massFunctionObserved')->get();
	my $covariance           = $observed->dataset('covariance'          )->get();
	my $errorObserved        = sqrt($covariance->diagonal(0,1));
	# Create a plot of this.
	my $plot;
	my $gnuPlot;
	my $plotFile = $workDirectory."stage".$stage."/massFunctionBestFit.pdf";
	(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
	print $gnuPlot "set title 'Best fit mass function for Stage ".$stage."'\n";
	print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\,[M_\\odot]\$'\n";
	print $gnuPlot "set ylabel 'Mass function; \$\\phi(M)\$ [Mpc\$^{-3} \\log(M_\\star)\$]'\n";
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
	print $gnuPlot "set xrange [1.0e8:1.0e13]\n";
	print $gnuPlot "set yrange [1.0e-9:1.0e-1]\n";
	print $gnuPlot "set pointsize 2.0\n";
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $stellarMassObserved,
	    $massFunctionObserved,
	    errorUp    => $errorObserved,
	    errorDown  => $errorObserved,
	    style      => "point",
	    weight     => [5,3],
	    symbol     => [6,7],
	    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
	    title      => "Li \\\\& White (2009)"
	    );
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $stellarMass,
	    $massFunction,
	    style      => "point",
	    weight     => [5,3],
	    symbol     => [6,7],
	    pointSize  => 1.0,
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
