#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::Stats::Basic;
use Stats::Histograms;
use Stats::Percentiles;
use Galacticus::HDF5;
use Galacticus::StellarMass;
use Galacticus::GasMass;
use Galacticus::Options;
use Galacticus::Launch::PBS;

# Run Galacticus models for integration testing.
# Andrew Benson (05-December-2014)

# Get arguments.
my %options = (
    calibrate           => "no" ,
    calibrateCount      => 100  ,
    calibratePercentile => 1.0  ,
    instance            => "1:1"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Check for an instance number for this launch.
my $instance;
my $instanceCount;
if ( $options{"instance"} =~ m/(\d+):(\d+)/ ) {
    $instance      = $1;
    $instanceCount = $2;
    print " -> launching instance ".$instance." of ".$instanceCount."\n";
} else {
    die("'instance' argument syntax error");
}

# Specify test statistics to compute.
my @tests =
    (
     {
	 name  => "halo mass function z=0.0",
	 label => "massFunctionHalo",
	 test  => \&massFunctionHalo
     },
     {
	 name  => "stellar mass function z=0.0",
	 label => "massFunctionStellar",
	 test  => \&massFunctionStellar
     },
     {
	 name  => "ISM gas mass function z=0.0",
	 label => "massFunctionISM",
	 test  => \&massFunctionISM
     },
     {
	 name  => "median disk sizes z=0.0",
	 label => "medianSizes",
	 test  => \&medianSizes
     },
    );

# Get Git revision.
my $gitRevision;
open(my $gitHndl,"git rev-parse HEAD|");
$gitRevision = <$gitHndl>;
close($gitHndl);
$gitRevision = "Unknown"
    unless ( defined($gitRevision) );
chomp($gitRevision);

# Choose a random seed.
my $randomSeed = int(rand(1000));

# Get PBS configuration and determine number of threads.
my $pbsConfig = &Galacticus::Options::Config("pbs");
my $ppn       = exists($pbsConfig->{'ppn'}) ? $pbsConfig->{'ppn'} : 1;

# Find integration models to run.
my $modelCount = 0;
opendir(my $testSuite,".");
while ( my $fileName = readdir($testSuite) ) {
    # Skip non-model integration files.
    next
	unless ( $fileName =~ m/^test\-model\-integration\-([a-zA-Z0-9]+)\.xml$/ );
    my $modelName = $1;    
    # Determine whether to run this model.
    ++$modelCount;
    next
	unless ( $modelCount % $instanceCount == $instance-1 );
    print "Running model '".$modelName."'\n";
    # Iterate over realizations.
    my @pbsJobs;
    for(my $i=0;$i<($options{'calibrate'} eq "yes" ? $options{'calibrateCount'} : 1);++$i) {
    	# Run the model.
    	print "--> Generating model '".$modelName."' ".($options{'calibrate'} eq "yes" ? "[realization: ".$i."]" : "")." for model integration testing...\n";
    	print "   --> ".($options{'calibrate'} eq "yes" ? "Queueing" : "Running")." model...\n";
    	system("mkdir -p outputs/test-model-integration/".$modelName.($options{'calibrate'} eq "yes" ? $i : ""));
    	my $xml = new XML::Simple();
    	my $parameters = $xml->XMLin($fileName);
    	if ( $options{'calibrate'} eq "yes" ) {
    	    $parameters->{'galacticusOutputFileName'}          ->{'value'} = "testSuite/outputs/test-model-integration/".$modelName.$i."/galacticus.hdf5";
    	    $parameters->{'randomNumberGenerator'   }->{'seed'}->{'value'} = $randomSeed+$i;
    	} else {
    	    $parameters->{'galacticusOutputFileName'}          ->{'value'} = "testSuite/outputs/test-model-integration/".$modelName   ."/galacticus.hdf5";
    	    $parameters->{'randomNumberGenerator'   }->{'seed'}->{'value'} = $randomSeed+$i;
    	}
    	open(my $parameterFile,">outputs/test-model-integration/".$modelName.($options{'calibrate'} eq "yes" ? $i : "")."/parameters.xml");
    	print $parameterFile $xml->XMLout($parameters, RootName => "parameters");
    	close($parameterFile);
    	if ( $options{'calibrate'} eq "yes" ) {
    	    my %job =
    		(
    		 launchFile   => "outputs/test-model-integration/".$modelName.$i."/launch.pbs",
    		 label        => $modelName.$i,
    		 logFile      => "outputs/test-model-integration/".$modelName.$i."/launch.log",
    		 command      => "cd ..; ./Galacticus.exe testSuite/outputs/test-model-integration/".$modelName.$i."/parameters.xml",
    		 ppn          => $ppn
    		);
    	    push(@pbsJobs,\%job);
    	} else {
    	    system("cd ..; ./Galacticus.exe testSuite/outputs/test-model-integration/".$modelName."/parameters.xml");
    	}
    	print "   <-- ...done\n";
    	print "<-- ...done\n";
    }
    &Galacticus::Launch::PBS::SubmitJobs(\%options,@pbsJobs)
    	if ( $options{'calibrate'} eq "yes" );
    # If calibrating, aggregate statistics and compute means and variances.
    if ( $options{'calibrate'} eq "yes" ) {
	for(my $i=0;$i<$options{'calibrateCount'};++$i) {
	    # Generate test statistics.
	    print "--> Testing model '".$modelName."' [realization: ".$i."] for model integration testing...\n";
	    foreach my $test ( @tests ) {
		print "   --> Test '".$test->{'name'}."'...\n";
		(my $failureCount, my $testCount) = &{$test->{'test'}}($modelName,"outputs/test-model-integration/".$modelName.$i,$test->{'label'},$gitRevision);
	    }
	    print "<-- ...done\n";
	}
	foreach my $test ( @tests ) {
	    print "--> Calibrating test '".$test->{'name'}."'...\n";
	    my $calibration;
	    for(my $i=0;$i<$options{'calibrateCount'};++$i) {
		(my $values) = rcols("outputs/test-model-integration/".$modelName.$i."/".$test->{'label'}."_r".$gitRevision.".txt",1);
		$calibration = pdl zeroes($options{'calibrateCount'},nelem($values))
		    unless ( defined($calibration) );
		$calibration->(($i),:) .= $values;
	    }
	    my $percentiles       = pdl [ $options{'calibratePercentile'}, 50.0, 100.0-$options{'calibratePercentile'} ];
	    my $calibrationMedian = pdl zeros($calibration->dim(1));
	    my $calibrationLow    = pdl zeros($calibration->dim(1));
	    my $calibrationHigh   = pdl zeros($calibration->dim(1));
	    for(my $i=0;$i<nelem($calibrationLow);++$i) {
		my $rank       = $calibration->(:,($i))->qsorti();
		my $percentile = pdl 100.0*(sequence($options{'calibrateCount'})+1.0)/$options{'calibrateCount'};
		(my $interpolants) = interpolate($percentiles,$percentile,$calibration->($rank,($i)));
		$calibrationLow   ->(($i)) .= $interpolants->((0));
		$calibrationMedian->(($i)) .= $interpolants->((1));
		$calibrationHigh  ->(($i)) .= $interpolants->((2));
	    }
	    open(my $calibrationFile,">data/model-integration/".$test->{'label'}."_".$modelName.".txt");
	    print $calibrationFile "# Model integration test: ".$test->{'name'}." [Git revision: ".$gitRevision."]\n";
	    for(my $i=0;$i<nelem($calibrationLow);++$i) {
		print $calibrationFile $i."\t".$calibrationLow->(($i))."\t".$calibrationMedian->(($i))."\t".$calibrationHigh->(($i))."\n";
	    }
	    close($calibrationFile);
	    print "<-- ...done\n";
	}
    } else {
	# Test models.
	my $failureTotal;
	my $testTotal;
	# Generate test statistics.
	print "--> Testing model '".$modelName."' for model integration testing...\n";
	foreach my $test ( @tests ) {
	    print "   --> Test '".$test->{'name'}."'...\n";
	    (my $failureCount, my $testCount) = &{$test->{'test'}}($modelName,"outputs/test-model-integration/".$modelName,$test->{'label'},$gitRevision);
	    $failureTotal += $failureCount;
	    $testTotal    += $testCount;
	    print "   --> failure rate: ".$failureCount."/".$testCount."\n";
	}
	print "<-- ...done\n";
	# Compute probability of failure count.
	my $probabilityFailure = 2.0*$options{'calibratePercentile'}/100.0;
	my $probability        = 0.0;
	for(my $k=0;$k<=$failureTotal;++$k) {
	    $probability += &binomialCoefficient($testTotal,$k)*$probabilityFailure**$k*(1.0-$probabilityFailure)**($testTotal-$k);
	}
	print "Probability of this number or fewer failures = ".$probability."\n";
	my $probabilityExcess = 1.0-$probability;
	my $status;
	if      ( $probabilityExcess < 0.001 ) {
	    $status = "FAILED";
	} elsif ( $probabilityExcess < 0.010 ) {
	    $status = "WARNING";
	} else {
	    $status = "success";
	}
	print "Status: ".$status."\n";
    }
}
closedir($testSuite);    
exit 0;

sub binomialCoefficient {
    # Return the binomial coefficient (k,n).
    my $r=1;$r*=$_/($_[0]-$_+1)for(1+pop..$_[0]);$r
}

sub massFunctionHalo {
    my $modelName          = shift();
    my $modelDirectoryName = shift();
    my $label              = shift();
    my $gitRevision        = shift();
    # Create mass function bins.
    my $massHaloLogarithmicBins = pdl sequence(10)/2.0+9.0;
    # Create data structure to read the results.
    my $model;
    $model->{'file' } = $modelDirectoryName."/galacticus.hdf5";
    $model->{'store'} = 0;
    $model->{'tree' } = "all";
    &Galacticus::HDF5::Get_Parameters($model    );
    &Galacticus::HDF5::Get_Times     ($model    );
    &Galacticus::HDF5::Select_Output ($model,0.0);
    &Galacticus::HDF5::Count_Trees   ($model    );
    &Galacticus::HDF5::Get_Dataset($model,['mergerTreeWeight','basicMass']);
    my $massHaloLogarithmic = log10($model->{'dataSets'}->{'basicMass'       });
    my $weight              =       $model->{'dataSets'}->{'mergerTreeWeight'} ;
    (my $massFunction, my $massFunctionError) = &Stats::Histograms::Histogram($massHaloLogarithmicBins,$massHaloLogarithmic,$weight,differential => 1);
    $massFunction      /= log(10.0);
    $massFunctionError /= log(10.0);
    # Output the mass function.
    open(my $outputFile,">".$modelDirectoryName."/".$label."_r".$gitRevision.".txt");
    print $outputFile "# Model integration test: Halo mass function at z=0.0 [git revision: ".$gitRevision."]\n";
    for(my $i=0;$i<nelem($massHaloLogarithmicBins);++$i) {
	print $outputFile $i."\t".$massFunction->(($i))."\n";
    }
    close($outputFile);
    # Read the reference dataset.
    my $referenceDataFileName = "data/model-integration/".$label."_".$modelName.".txt";
    if ( -e $referenceDataFileName ) {
	(my $referenceLow, my $referenceHigh) = rcols($referenceDataFileName,1,3);
	my $failures = which(($massFunction < $referenceLow*(1.0-1.0e-6)) | ($massFunction > $referenceHigh*(1.0+1.0e-6)));
	return ( nelem($failures), nelem($massFunction) );
    } else {
	return ( 0               , nelem($massFunction) );
    }
}

sub massFunctionStellar {
    my $modelName          = shift();
    my $modelDirectoryName = shift();
    my $label              = shift();
    my $gitRevision        = shift();
    # Create mass function bins.
    my $massStellarLogarithmicBins = pdl sequence(10)/2.0+8.0;
    # Create data structure to read the results.
    my $model;
    $model->{'file' } = $modelDirectoryName."/galacticus.hdf5";
    $model->{'store'} = 0;
    $model->{'tree' } = "all";
    &Galacticus::HDF5::Get_Parameters($model    );
    &Galacticus::HDF5::Get_Times     ($model    );
    &Galacticus::HDF5::Select_Output ($model,0.0);
    &Galacticus::HDF5::Count_Trees   ($model    );
    &Galacticus::HDF5::Get_Dataset($model,['mergerTreeWeight','massStellar']);
    my $massStellarLogarithmic = log10($model->{'dataSets'}->{'massStellar'});
    my $weight                 = $model->{'dataSets'}->{'mergerTreeWeight'};
    (my $massFunction, my $massFunctionError) = &Stats::Histograms::Histogram($massStellarLogarithmicBins,$massStellarLogarithmic,$weight,differential => 1);
    $massFunction      /= log(10.0);
    $massFunctionError /= log(10.0);
    # Output the mass function.
    open(my $outputFile,">".$modelDirectoryName."/".$label."_r".$gitRevision.".txt");
    print $outputFile "# Model integration test: Stellar mass function at z=0.0 [git revision: ".$gitRevision."]\n";
    for(my $i=0;$i<nelem($massStellarLogarithmicBins);++$i) {
	print $outputFile $i."\t".$massFunction->(($i))."\n";
    }
    close($outputFile);
    # Read the reference dataset.
    my $referenceDataFileName = "data/model-integration/".$label."_".$modelName.".txt";
    if ( -e $referenceDataFileName ) {
	(my $referenceLow, my $referenceHigh) = rcols($referenceDataFileName,1,3);
	my $failures = which(($massFunction < $referenceLow*(1.0-1.0e-6)) | ($massFunction > $referenceHigh*(1.0+1.0e-6)));
	return ( nelem($failures), nelem($massFunction) );
    } else {
	return ( 0               , nelem($massFunction) );
    }
}

sub massFunctionISM {
    my $modelName          = shift();
    my $modelDirectoryName = shift();
    my $label              = shift();
    my $gitRevision        = shift();
    # Create mass function bins.
    my $massColdGasLogarithmicBins = pdl sequence(10)/2.0+6.0;
    # Create data structure to read the results.
    my $model;
    $model->{'file' } = $modelDirectoryName."/galacticus.hdf5";
    $model->{'store'} = 0;
    $model->{'tree' } = "all";
    &Galacticus::HDF5::Get_Parameters($model    );
    &Galacticus::HDF5::Get_Times     ($model    );
    &Galacticus::HDF5::Select_Output ($model,0.0);
    &Galacticus::HDF5::Count_Trees   ($model    );
    &Galacticus::HDF5::Get_Dataset($model,['mergerTreeWeight','massColdGas']);
    my $massColdGasLogarithmic = log10($model->{'dataSets'}->{'massColdGas'});
    my $weight                 = $model->{'dataSets'}->{'mergerTreeWeight'};
    (my $massFunction, my $massFunctionError) = &Stats::Histograms::Histogram($massColdGasLogarithmicBins,$massColdGasLogarithmic,$weight,differential => 1);
    $massFunction      /= log(10.0);
    $massFunctionError /= log(10.0);
    # Output the mass function.
    open(my $outputFile,">".$modelDirectoryName."/".$label."_r".$gitRevision.".txt");
    print $outputFile "# Model integration test: ISM mass function at z=0.0 [git revision: ".$gitRevision."]\n";
    for(my $i=0;$i<nelem($massColdGasLogarithmicBins);++$i) {
	print $outputFile $i."\t".$massFunction->(($i))."\n";
    }
    close($outputFile);
    # Read the reference dataset.
    my $referenceDataFileName = "data/model-integration/".$label."_".$modelName.".txt";
    if ( -e $referenceDataFileName ) {
	(my $referenceLow, my $referenceHigh) = rcols($referenceDataFileName,1,3);
	my $failures = which(($massFunction < $referenceLow*(1.0-1.0e-6)) | ($massFunction > $referenceHigh*(1.0+1.0e-6)));
	return ( nelem($failures), nelem($massFunction) );
    } else {
	return ( 0               , nelem($massFunction) );
    }
}

sub medianSizes {
    my $modelName          = shift();
    my $modelDirectoryName = shift();
    my $label              = shift();
    my $gitRevision        = shift();
    # Create mass function bins.
    my $massStellarLogarithmicBins = pdl sequence(10)/2.0+8.0;
    # Create data structure to read the results.
    my $model;
    $model->{'file' } = $modelDirectoryName."/galacticus.hdf5";
    $model->{'store'} = 0;
    $model->{'tree' } = "all";
    &Galacticus::HDF5::Get_Parameters($model    );
    &Galacticus::HDF5::Get_Times     ($model    );
    &Galacticus::HDF5::Select_Output ($model,0.0);
    &Galacticus::HDF5::Count_Trees   ($model    );
    &Galacticus::HDF5::Get_Dataset($model,['mergerTreeWeight','massStellar','diskRadius']);
    my $size                   =       $model->{'dataSets'}->{'diskRadius'      } ;
    my $massStellarLogarithmic = log10($model->{'dataSets'}->{'massStellar'     });
    my $weight                 =       $model->{'dataSets'}->{'mergerTreeWeight'};
    my $percentiles            = pdl [50.0];
    my $quantiles              = &Stats::Percentiles::BinnedPercentiles($massStellarLogarithmicBins,$massStellarLogarithmic,$size,$weight,$percentiles);
    # Output the size distribution.
    open(my $outputFile,">".$modelDirectoryName."/".$label."_r".$gitRevision.".txt");
    print $outputFile "# Model integration test: Median disk sizes at z=0.0 [git revision: ".$gitRevision."]\n";
    for(my $i=0;$i<nelem($massStellarLogarithmicBins);++$i) {
    	print $outputFile $i."\t".$quantiles->(($i),(0))."\n";
    }
    close($outputFile);
    # Read the reference dataset.
    my $referenceDataFileName = "data/model-integration/".$label."_".$modelName.".txt";
    if ( -e $referenceDataFileName ) {
	(my $referenceLow, my $referenceHigh) = rcols($referenceDataFileName,1,3);
	my $failures = which(($quantiles->(:,(0)) < $referenceLow*(1.0-1.0e-6)) | ($quantiles->(:,(0)) > $referenceHigh*(1.0+1.0e-6)));
	return ( nelem($failures), $quantiles->dim(0) );
    } else {
	return ( 0               , $quantiles->dim(0) );
    }
}
