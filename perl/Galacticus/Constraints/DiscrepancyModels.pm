# Contains a Perl module which provides tools for running models for discrepancy analysis.
# Andrew Benson (15-March-2015)

package DiscrepancyModels;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath  = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use Clone qw(clone);
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::DiscrepancySystematics;
require Galacticus::Launch::PBS;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

sub RunModels {
    # Get arguments.
    my $discrepancyName =   shift() ;
    my $description     =   shift() ;
    my $configFile      =   shift() ;
    my %options         = %{shift()};
    my $models          =   shift() ;
    
    # Parse the Galacticus config file.
    my $xml    = new XML::Simple;
    my $config = $xml->XMLin($configFile, KeyAttr => 0);

    # Validate the config file.
    die("DiscrepancyModels::RunModels(): workDirectory must be specified in config file" )
	unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
    die("DiscrepancyModels::RunModels(): compilation must be specified in config file"   )
	unless ( exists($config->{'likelihood'}->{'compilation'   }) );
    die("DiscrepancyModels::RunModels(): baseParameters must be specified in config file")
	unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

    # Determine the scratch and work directories.
    my $workDirectory    =        $config->{'likelihood'}->{'workDirectory'   };
    my $scratchDirectory = exists($config->{'likelihood'}->{'scratchDirectory'})
	                   ?
                              	  $config->{'likelihood'}->{'scratchDirectory'}
                           :
                   	          $config->{'likelihood'}->{'workDirectory'   };
    
    # Create the work and scratch directories.
    system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'});

    # Determine base parameters to use.
    my $baseParameters = exists($options{'baseParameters'}) 
	                 ?
	                        $options{'baseParameters'}
                         :
	                        $config->{'likelihood'}->{'baseParameters'};

    # Ensure that Galacticus is built.
    if ( $options{'make'} eq "yes" ) {
	system("make Galacticus.exe");
	die("DiscrepancyModels::RunModels(): failed to build Galacticus.exe")
	    unless ( $? == 0 );
    }
    
    # Get a hash of the parameter values.
    (my $constraintsRef, my $parameters) = 
	&Parameters::Compilation
	(
	 $config->{'likelihood'}->{'compilation'},
	 $baseParameters
	);
    my @constraints = @{$constraintsRef};
    
    # Switch off thread locking.
    $parameters->{'treeEvolveThreadLock'}->{'value'} = "false";

    # Initialize a stack for PBS models.
    my @pbsStack;

    # Iterate over models.
    foreach my $model ( "default", "alternate" ) {
	next
	    if ( exists($options{'runOnly'}) && $options{'runOnly'} ne $model );
	# Specify the output name.
	my $modelDirectory = $workDirectory."/modelDiscrepancy/".$discrepancyName."/".$models->{$model}->{'label'};
	system("mkdir -p ".$modelDirectory);
	my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
	push(
	    @{$models->{$model}->{'parameters'}},
	    {
		name  => "galacticusOutputFileName",
		value => $galacticusFileName
	    }
	    );
	my $newParameters = clone($parameters);
	foreach my $parameter ( @{$models->{$model}->{'parameters'}} ) {
	    my $thisParameter = $newParameters;
	    foreach my $thisParameterName ( split(/\-\>/,$parameter->{'name'}) ) {
		my $arrayParameter = 0;
		if ( $thisParameterName =~ m/^\@(.*)/ ) {
		    $thisParameterName = $1;
		    $arrayParameter    =  1;
		}
		if ( $arrayParameter ) {
		    ${$thisParameter->{$thisParameterName}}[++$#{$thisParameter->{$thisParameterName}}]->{'value'} = undef();
		} else {
		    $thisParameter->{$thisParameterName}->{'value'} = undef()
			unless ( exists($thisParameter->{$thisParameterName}) );
		}
		if ( UNIVERSAL::isa($thisParameter->{$thisParameterName},"ARRAY") ) {
		    $thisParameter = ${$thisParameter->{$thisParameterName}}[$#{$thisParameter->{$thisParameterName}}];
		} else {
		    $thisParameter = $thisParameter->{$thisParameterName};
		}
	    }
	    $thisParameter->{'value'} = $parameter->{'value'};
	}
	# Adjust the number of trees to run if specified, unless the model specifically requests that trees be read from file.
	if ( exists($options{'treesPerDecade'}) && ! grep {$_->{'name'} eq "mergerTreeConstructMethod" && $_->{'value'} eq "read"} @{$models->{$model}->{'parameters'}} ) {
	    # Must also specify that trees are to be built in this case.
	    $newParameters->{'mergerTreeConstructMethod'    }->{'value'} = "build";
	    $newParameters->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $options{'treesPerDecade'};
	}
	# Run the model.
	unless ( -e $galacticusFileName ) {
	    # Generate the parameter file.
	    my $parameterFileName = $modelDirectory."/parameters.xml";
	    &Parameters::Output($newParameters,$parameterFileName);
	    # Create PBS job.
	    my $command = "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	    foreach my $constraint ( @constraints ) {
		# Parse the definition file.
		my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
		# Insert code to run the analysis code.
		my $analysisCode = $constraintDefinition->{'analysis'};
		$command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5\n";
	    }
	    my %job =
		(
		 launchFile => $modelDirectory."/launch.pbs",
		 label      => $discrepancyName.ucfirst($models->{$model}->{'label'}),
		 logFile    => $modelDirectory."/launch.log",
		 command    => $command
		);
	    foreach ( 'ppn', 'walltime', 'mem' ) {
		if ( exists($models->{$model}->{$_}) ) {
		    $job{$_} = $models->{$model}->{$_};
		} elsif ( exists($options{$_}) ) {
		    $job{$_} = $options{$_};
		}
	    }
	    # Queue the calculation.
	    push(
		@pbsStack,
		\%job
		);   
	}
    }
    # Send jobs to PBS.
    &PBS::SubmitJobs(\%options,@pbsStack)
	if ( scalar(@pbsStack) > 0 );
    # Iterate over constraints.
    unless ( exists($options{'analyze'}) && $options{'analyze'} eq "no" ) {
	foreach my $constraint ( @constraints ) {
	    # Parse the definition file.
	    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
	    print "Computing discrepancy for constraint: ".$constraintDefinition->{'label'}."\n";    
	    # Locate the model results.
	    my $alternateResultFileName = 
		$workDirectory                   .
		"/modelDiscrepancy/"             .
		$discrepancyName                 .
		"/"                              .
		$models->{'alternate'}->{'label'}.
		"/"                              .
		$constraintDefinition->{'label'} .
		".hdf5";
	    my $defaultResultFileName     = 
		$workDirectory                  .
		"/modelDiscrepancy/"            .
		$discrepancyName                .
		"/"                             .
		$models->{'default'}->{'label'} .
		"/"                             .
		$constraintDefinition->{'label'}.
		".hdf5";
	    # Read the results.
	    my $alternateResult = new PDL::IO::HDF5($alternateResultFileName);
	    my $defaultResult   = new PDL::IO::HDF5($defaultResultFileName  );
	    # Extract the results.
	    my $defaultX            = $defaultResult  ->dataset('x'         )->get();
	    my $defaultY            = $defaultResult  ->dataset('y'         )->get();
	    my $defaultCovariance   = $alternateResult->dataset('covariance')->get();
	    my $alternateY          = $alternateResult->dataset('y'         )->get();
	    my $alternateCovariance = $alternateResult->dataset('covariance')->get();
	    # Apply any systematics models.
	    my %systematicResults;
	    foreach my $argument ( keys(%options) ) {
		if ( $argument =~ m/^systematic(.*)/ ) {
		    my $systematicModel = $1;
		    if 
			( 
			  exists($DiscrepancySystematics::models{$systematicModel}) 
			  && 
			  $options{$argument} eq "yes"
			  &&
			  exists($constraintDefinition->{$argument})
			  &&
			  $constraintDefinition->{$argument} eq "yes"
			) {
			    %{$systematicResults{$systematicModel}} =
				&{$DiscrepancySystematics::models{$systematicModel}}(
				\%options             ,
				$constraintDefinition ,
				$defaultX             ,
				$defaultY             ,
				$defaultCovariance    ,
				$alternateY         ,
				$alternateCovariance
			    );
		    }
		}
	    }
	    # Find the multiplicative discrepancy between these two models.
	    (my $nonZero, my $zero)                      = which_both($defaultY > 0.0);
	    my $modelDiscrepancyMultiplicative           = $alternateY->copy();
	    $modelDiscrepancyMultiplicative->($nonZero) /= $defaultY->($nonZero);
	    $modelDiscrepancyMultiplicative->(   $zero) .= 1.0
		if ( nelem($zero) > 0 );
	    # Compute the covariance.
	    my $modelDiscrepancyCovarianceMultiplicative = 
		+outer($defaultY-$alternateY,$defaultY-$alternateY)
		*outer(      1.0/$alternateY,      1.0/$alternateY);
	    # Output the model discrepancy to file.
	    my $outputFile = 
		new PDL::IO::HDF5(
		    ">".
		    $workDirectory.
		    "/modelDiscrepancy/".$discrepancyName."/discrepancy".
		    ucfirst($constraintDefinition->{'label'}).
		    ".hdf5"
		);
	    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
	    $outputFile->attrSet(
		description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to ".$description."."
		);
	    # Add results of systematics models.
	    my $systematicGroup = $outputFile->group("systematicModels");
	    foreach my $model ( keys(%systematicResults) ) {
		my $modelGroup = $systematicGroup->group($model);
		my %modelResults = %{$systematicResults{$model}};
		foreach my $parameter ( keys(%modelResults) ) {
		    $modelGroup->attrSet($parameter => $modelResults{$parameter});
		}
	    }
	    # Create a plot if requested.
	    if ( $options{'plot'} eq "yes" ) {
		# Determine suitable x and y ranges.
		my $yCombined = $alternateY->append($defaultY);
		my $yNonZero  = which($yCombined > 0.0);
		my $x    = $defaultX;
		my $xMin = 0.67*min($x     );
		my $xMax = 1.50*max($x     );
		my $yMin = 0.67*min($yCombined->($yNonZero));
		my $yMax = 1.50*max($yCombined->($yNonZero));
		$yMin = 1.0e-10*$yMax
		    if ( $yMin < 1.0e-10*$yMax );
		# Declare standards for GnuPlot;
		my ($gnuPlot, $plotFileEPS, $plot);
		# Open a pipe to GnuPlot.
		(my $label = ucfirst($constraintDefinition->{'label'})) =~ s/\./_/;
		my $plotFile = $workDirectory."/modelDiscrepancy/".$discrepancyName."/discrepancy".$label.".pdf";
		($plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
		open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
		print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
		print $gnuPlot "set output '".$plotFileEPS."'\n";
		print $gnuPlot "set lmargin screen 0.15\n";
		print $gnuPlot "set rmargin screen 0.95\n";
		print $gnuPlot "set bmargin screen 0.15\n";
		print $gnuPlot "set tmargin screen 0.95\n";
		print $gnuPlot "set key spacing 1.2\n";
		print $gnuPlot "set key at screen 0.4,0.2\n";
		print $gnuPlot "set key left\n";
		print $gnuPlot "set key bottom\n";
		print $gnuPlot "set logscale xy\n";
		print $gnuPlot "set mxtics 10\n";
		print $gnuPlot "set mytics 10\n";
		print $gnuPlot "set format x '\$10^{\%L}\$'\n";
		print $gnuPlot "set format y '\$10^{\%L}\$'\n";
		print $gnuPlot "set xrange [".$xMin.":".$xMax."]\n";
		print $gnuPlot "set yrange [".$yMin.":".$yMax."]\n";
		print $gnuPlot "set xlabel '\$x\$'\n";
		print $gnuPlot "set ylabel '\$y\$'\n";	
		print $gnuPlot "set title '".$constraintDefinition->{'label'}."'\n";
		&PrettyPlots::Prepare_Dataset(\$plot,
					      $x,$alternateY,
					      errorUp   => sqrt($alternateCovariance->diagonal(0,1)),
					      errorDown => sqrt($alternateCovariance->diagonal(0,1)),
					      style     => "point",
					      symbol    => [6,7], 
					      weight    => [5,3],
					      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
					      title     => $models->{'alternate'}->{'label'}
		    );
		&PrettyPlots::Prepare_Dataset(\$plot,
					      $x,$defaultY,
					      errorUp   => sqrt($defaultCovariance->diagonal(0,1)),
					      errorDown => sqrt($defaultCovariance->diagonal(0,1)),
					      style     => "point",
					      symbol    => [6,7], 
					      weight    => [5,3],
					      color     => $PrettyPlots::colorPairs{'mediumSeaGreen'},
					      title     => $models->{'default'}->{'label'}
		    );
		&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
		close($gnuPlot);
		&LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
	    }
	}
    }   
}

sub Apply_Discrepancies {
    my $discrepancyFileName = shift();
    my $discrepancyPath     = shift();
    my $value               = shift();
    my $error               = shift();
    my $covariance          = shift();
    # Get any options.
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );
    # Scan the path for discrepancy files.
    opendir(discrepDir,$discrepancyPath);
    while ( my $discrepancy = readdir(discrepDir) ) {
	my $discrepancyFileName = $discrepancyPath."/".$discrepancy."/".$discrepancyFileName;
	if ( -e $discrepancyFileName ) {
	    my $discrepancyFile = new PDL::IO::HDF5($discrepancyFileName);
	    my @datasets = $discrepancyFile->datasets();
	    foreach my $dataset ( @datasets ) {
		if ( $dataset eq "multiplicative" ) {
		    # Read the multiplicative discrepancy
		    my $multiplier  = $discrepancyFile->dataset('multiplicative')->get();
		    $value         *= $multiplier;
		    $error         *= $multiplier;
		    $covariance    .= $covariance*outer($multiplier,$multiplier);
		}		    
		if ( $dataset eq "multiplicativeCovariance" ) {
		    # Adjust the model accordingly.
		    my $covarianceMultiplier  = $discrepancyFile->dataset('multiplicativeCovariance')->get();
		    if ( exists($options{'limitMultiplicativeCovariance'}) ) {
			# First ensure that no diagonal terms in the covariance matrix exceed the allowed limit. If they do,
			# truncate to that limit while preserving the correlation structure of the matrix.
			my $variance     = $covarianceMultiplier->diagonal(0,1);
			my $highVariance = which($variance > $options{'limitMultiplicativeCovariance'});
			my $correlation  = $covarianceMultiplier/outer($variance->sqrt(),$variance->sqrt());
			$variance->($highVariance) .= $options{'limitMultiplicativeCovariance'};
			$covarianceMultiplier .= $correlation*outer($variance->sqrt(),$variance->sqrt());
			# Check for any off-diagonal terms that exceed the limit and limit them.
			my $low  = which($covarianceMultiplier->flat() < -$options{'limitMultiplicativeCovariance'});
			my $high = which($covarianceMultiplier->flat() > +$options{'limitMultiplicativeCovariance'});
			$covarianceMultiplier->flat()->($low ) .= -$options{'limitMultiplicativeCovariance'};
			$covarianceMultiplier->flat()->($high) .= +$options{'limitMultiplicativeCovariance'};
		    }
		    $covariance              += $covarianceMultiplier*outer($value,$value);
		}		    
		if ( $dataset eq "additive" ) {
		    # Read the additive discrepancy
		    my $addition  = $discrepancyFile->dataset('additive')->get();
		    # Adjust the model accordingly.
		    $value       += $addition;
		}
		if ( $dataset eq "additiveCovariance" ) {
		    # Read the covariance of the discrepancy.
		    my $covariance  = $discrepancyFile->dataset('additiveCovariance')->get();
		    # Adjust the model discrepancy covariance accordingly.
		    $covariance    += $covariance;
		}
	    }
	}
    }
}


1;
