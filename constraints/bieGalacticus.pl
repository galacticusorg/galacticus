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
use UNIVERSAL;
use Data::Dumper;
require Galacticus::Constraints::Parameters;

# Generate an input file for the Bayesian Inference Engine to run Galacticus.
# Andrew Benson (04-October-2011)

# Get arguments.
die("Usage: bieGalacticus.pl <configFile> [options]") unless ( scalar(@ARGV) > 0 );
my $configFile = $ARGV[0];

# Create a hash of named arguments.
my $iArg = 0;
my %arguments = (

    make => "no",
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Validate the config file.
die("bieGalacticus.pl: workDirectory must be specified in config file") unless ( exists($config->{'workDirectory'}) );
die("bieGalacticus.pl: compilation must be specified in config file"  ) unless ( exists($config->{'compilation'  }) );
die("bieGalacticus.pl: parameters must be specified in config file"   ) unless ( exists($config->{'parameters'   }) );
die("bieGalacticus.pl: nodes must be specified in config file"        ) unless ( exists($config->{'nodes'        }) );

# Determine number of threads to run.
my $threadsPerNode = 1;
$threadsPerNode = $config->{'threadsPerNode'}
    if ( exists($config->{'threadsPerNode'}) );
my $nodeCount = 1;
$nodeCount = $config->{'nodes'}
    if ( exists($config->{'nodes'}) );
my $biePerNode = $threadsPerNode;
$biePerNode = $config->{'bieThreads'}
   if ( exists($config->{'bieThreads'}) );
my $bieThreadsTotal = $biePerNode*$nodeCount;

# Determine the scratch directory.
my $scratchDirectory = $config->{'workDirectory'}."/mcmc";
$scratchDirectory = $config->{'scratchDirectory'} if ( exists($config->{'scratchDirectory'}) );

# Determine the BIE executable to use.
my $bie = "bie";
$bie = $config->{'bie'} if ( exists($config->{'bie'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'workDirectory'}."/mcmc");
system("cp -f ".$configFile." ".$config->{'workDirectory'}."/mcmc/config.xml");

# Parse parameters from the config file.
my @parameters;
if ( UNIVERSAL::isa($config->{'parameters'}->{'parameter'},"ARRAY") ) {
    @parameters = @{$config->{'parameters'}->{'parameter'}};
} else {
    push(@parameters,$config->{'parameters'}->{'parameter'});
}
my $parameterCount = 0;
for(my $i=0;$i<scalar(@parameters);++$i) {
    ++$parameterCount if ( exists($parameters[$i]->{'prior'}) );
}

# Detect which algorithm is to be used.
my $algorithm = "metropolisHastings";
$algorithm = $config->{'method'} if ( exists($config->{'method'}) );

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("constrainGalacticus.pl: failed to build Galacticus.exe") unless ( $? == 0 );
}

# Open the BIE input file.
open(oHndl,">".$config->{'workDirectory'}."/mcmc/bieInput.txt");

# Initialize the persistence system.
print oHndl "pnewsession ".$config->{'name'}."\n";
print oHndl "cktoggle\n";
my $checkpointInterval = 1;
$checkpointInterval = $config->{'checkpoint'}->{'interval'} if ( exists($config->{'checkpoint'}->{'interval'}) );
print oHndl "ckinterval ".$checkpointInterval."\n";

# Set names for output files.
print oHndl "set nametag = \"".$config->{'workDirectory'}."/mcmc/galacticusBIE\"\n";
print oHndl "set outfile = \"".$config->{'workDirectory'}."/mcmc/galacticusBIE.statelog\"\n";

# Select a likelihood function.
print oHndl "set ndim = ".$parameterCount."\n";
print oHndl "set fct = new GalacticusLikelihoodFunction(ndim,\"constraints/bieGalacticusWrapper.pl\",\"".$configFile."\",\"\%\%SCRATCHDIRECTORY\%\%/glcLikelihood.dat\")\n";

# Create a StateInfo and set labels.
print oHndl "set si = new StateInfo(ndim)\n";
print oHndl "set l = new clivectors(ndim)\n";
my $j = -1;
for(my $i=0;$i<scalar(@parameters);++$i) {
    if ( exists($parameters[$i]->{'prior'}) ) {
	++$j;
	print oHndl "l->setval(".$j.",\"".$parameters[$i]->{'name'}."\")\n";
    }
}
print oHndl "si->labelAll(l)\n";

# Construct the "eps" distributions (from which random perturbations to the proposal are drawn).
my $eps = 0.01;
$eps    = $config->{'monteCarlo'}->{'epsWidth'}
    if ( exists($config->{'monteCarlo'}->{'epsWidth'}) );
print oHndl "set eps = new clivectordist(".$parameterCount.")\n"
    if ( $algorithm eq "differentialEvolution" );

# Define prior distributions on parameters.
$j = -1;
for(my $i=0;$i<scalar(@parameters);++$i) {
    if ( exists($parameters[$i]->{'prior'}) ) {
	++$j;
	print oHndl "set dis".$j." = new ";
	# Map limits to logarithmic space if necessary.
	if ( $parameters[$i]->{'prior'}->{'distribution'} eq "uniform" ) {
	    if ( exists($parameters[$i]->{'prior'}->{'mapping'}) ) {
		if ( $parameters[$i]->{'prior'}->{'mapping'} eq "logarithmic" ) {
		    $parameters[$i]->{'prior'}->{'lowerLimit'} = log($parameters[$i]->{'prior'}->{'lowerLimit'});
		    $parameters[$i]->{'prior'}->{'upperLimit'} = log($parameters[$i]->{'prior'}->{'upperLimit'});
		}
	    }
	    print oHndl "UniformDist(".$parameters[$i]->{'prior'}->{'lowerLimit'}.", ".$parameters[$i]->{'prior'}->{'upperLimit'}.")\n";

	    # Define the eps distribution for this parameter if necessary.
	    if ( $algorithm eq "differentialEvolution" ) {
		my $parameterRange = $parameters[$i]->{'prior'}->{'upperLimit'}-$parameters[$i]->{'prior'}->{'lowerLimit'};
		my $width = $eps*$parameterRange;
		my $pWidth = sprintf("%.10f",$width);
		print oHndl "set cdf".$j." = new CauchyDist(".$pWidth.")\n";
		print oHndl "eps->setval(".$j.", cdf".$j.")\n";
	    }
	} elsif ( $parameters[$i]->{'prior'}->{'distribution'} eq "normal" ) {
	    my $variance = $parameters[$i]->{'prior'}->{'variance'};
	    $variance *= $config->{'temperature'}
		if ( exists($config->{'temperature'}) );
	    print oHndl "NormalDist(".$parameters[$i]->{'prior'}->{'mean'}.", ".$variance;
	    print oHndl ",".$parameters[$i]->{'prior'}->{'lowerLimit'}
	      if ( exists($parameters[$i]->{'prior'}->{'lowerLimit'}) );
	    print oHndl ",".$parameters[$i]->{'prior'}->{'upperLimit'}
	      if ( exists($parameters[$i]->{'prior'}->{'upperLimit'}) );       
	    print oHndl ")\n";

	    # Define the eps distribution for this parameter if necessary.
	    if ( $algorithm eq "differentialEvolution" ) {

		my $parameterRange;
		$parameterRange = $parameters[$i]->{'prior'}->{'upperLimit'}-$parameters[$i]->{'prior'}->{'lowerLimit'}
		    if ( exists($parameters[$i]->{'prior'}->{'upperLimit'}) && exists($parameters[$i]->{'prior'}->{'lowerLimit'}));
		if ( exists($parameters[$i]->{'prior'}->{'variance'}) ) {
		    my $dispersion = sqrt($parameters[$i]->{'prior'}->{'variance'});
		    $parameterRange = $dispersion
			if ( ! defined($parameterRange) || $parameterRange > $dispersion );
		}	       
		my $width = $eps*$parameterRange;
		my $pWidth = sprintf("%.10f",$width);
		print oHndl "set cdf".$j." = new CauchyDist(".$pWidth.")\n";
		print oHndl "eps->setval(".$j.", cdf".$j.")\n";
	    }
	} else {
	    die("bieGalacticus.pl: unknown prior distribution");
	}
    }
}

# Construct the prior distribution vector.
print oHndl "set pvec = new clivectordist(ndim)\n";
for(my $i=0;$i<$parameterCount;++$i) {
    print oHndl "pvec->setval(".$i.", dis".$i.")\n";
}
my $priorObjectName = "prior";
$priorObjectName = "oldPrior" if ( exists($config->{'postPrior'}) );
print oHndl "set ".$priorObjectName." = new Prior(si, pvec)\n";

# If we are to use a previously computed MCMC to estimate the posterior distribution of the parameters
# (presumably computed for some other constraint) as a prior, then set that up now.
if ( exists($config->{'postPrior'}) ) {
    my $burnCount = 0;
    $burnCount = $config->{'postPrior'}->{'burn'} if ( exists($config->{'postPrior'}->{'burn'}) );
    print oHndl "set postsi = new StateInfo(ndim)\n";
    print oHndl "set nburn  = ".$burnCount."\n";
    print oHndl "set post   = new EnsembleKD(postsi,0,nburn,\"".$config->{'postPrior'}->{'stateLog'}."\",0)\n";
    my $restartPrior = "no";
    $restartPrior = $config->{'postPrior'}->{'restart'}
       if ( exists($config->{'postPrior'}->{'restart'}) );
    if ( $restartPrior eq "yes" ) {
	print oHndl "set prior  = new RestartPrior(oldPrior,post)\n";
    } else {
	print oHndl "set prior  = new PostPrior(oldPrior,post)\n";
    }
}

# Create a Metropolis-Hastings helper class to define transition probability.
my $proposalWidth = 0.01;
$proposalWidth = $config->{'monteCarlo'}->{'proposalWidth'} if ( exists($config->{'monteCarlo'}->{'proposalWidth'}) );
print oHndl "set mvec = new clivectord(ndim, ".$proposalWidth.")\n";
print oHndl "set mhwidth = new MHWidthOne(si, mvec)\n";

# Create the Monte Carlo algorithm.
my $nSteps = 200;
$nSteps = $config->{'monteCarlo'}->{'steps'} if ( exists($config->{'monteCarlo'}->{'steps'}) );
print oHndl "set nsteps = ".$nSteps."\n";
my $statState = "si";
$statState = "postsi" if ( exists($config->{'postPrior'}) );
print oHndl "set sstat = new EnsembleStat(".$statState.")\n";
print oHndl "set convrg = new GelmanRubinConverge(0, sstat, \"convergence\")\n";
my $nStepsSkip = 0;
$nStepsSkip = $config->{'monteCarlo'}->{'stepsSkip'} if ( exists($config->{'monteCarlo'}->{'stepsSkip'}) );
print oHndl "convrg->setNskip(".$nStepsSkip.")\n";
my $nGood = 0;
$nGood = $config->{'monteCarlo'}->{'goodSteps'} if ( exists($config->{'monteCarlo'}->{'goodSteps'}) );
print oHndl "convrg->setNgood(".$nGood.")\n";
my $alpha = 0.05;
$alpha = $config->{'monteCarlo'}->{'grAlpha'} if ( exists($config->{'monteCarlo'}->{'grAlpha'}) );
print oHndl "convrg->setAlpha(".$alpha.")\n";
my $nOutlier = 500;
$nOutlier = $config->{'monteCarlo'}->{'grNOutlier'} if ( exists($config->{'monteCarlo'}->{'grNOutlier'}) );
print oHndl "convrg->setNoutlier(".$nOutlier.")\n";
my $maxOutlier = 6;
$maxOutlier = $config->{'monteCarlo'}->{'grMaxOut'} if ( exists($config->{'monteCarlo'}->{'grMaxOut'}) );
print oHndl "convrg->setMaxout(".$maxOutlier.")\n";
my $RhatMax = 1.2;
$RhatMax = $config->{'monteCarlo'}->{'grRhatMax'} if ( exists($config->{'monteCarlo'}->{'grRhatMax'}) );
print oHndl "convrg->setRhatMax(".$RhatMax.")\n";
print oHndl "set mca = new MetropolisHastings()\n";
print oHndl "set like = new LikelihoodComputationSerial()\n";
my $maxTemperature = 1.0;
$maxTemperature = $config->{'monteCarlo'}->{'maxTemperature'} if ( exists($config->{'monteCarlo'}->{'maxTemperature'}) );
if ( $maxTemperature > 1.0 ) {
    print oHndl "set Tmax = ".$maxTemperature."\n";
    my $minmc = 6;
    $minmc = $config->{'monteCarlo'}->{'minTemperatureStates'} if ( exists($config->{'monteCarlo'}->{'minTemperatureStates'}) );
    print oHndl "set minmc = ".$minmc."\n";
}
if ( $algorithm eq "metropolisHastings" ) {
    print oHndl "set sim = new ParallelChains(si, minmc, Tmax,  mhwidth, convrg, prior, like, mca)\n";
    print oHndl "sim->SetAlgorithm(1)\n";
    print oHndl "sim->NewNumber(".$bieThreadsTotal.")\n";
} elsif ( $algorithm eq "differentialEvolution" ) {
    if ( $maxTemperature <= 1.0 ) {
	print oHndl "set sim = new DifferentialEvolution(si, ".$bieThreadsTotal.", eps, convrg, prior, like, mca)\n";
    } else {
	print oHndl "set sim = new TemperedDifferentialEvolution(si, ".$bieThreadsTotal.", minmc, Tmax, eps, convrg, prior, like, mca)\n";
	my $tempFreq = 10;
	$tempFreq = $config->{'monteCarlo'}->{'temperingFrequency'}
	if ( exists($config->{'monteCarlo'}->{'temperingFrequency'}) );
	print oHndl "sim->SetTempFreq(".$tempFreq.")\n";
	my $tempSteps = 10;
	$tempSteps = $config->{'monteCarlo'}->{'temperingSteps'}
	if ( exists($config->{'monteCarlo'}->{'temperingSteps'}) );
	print oHndl "sim->SetEquilSteps(".$tempSteps.")\n";
	my $Tpow = 0.5;
	$Tpow = $config->{'monteCarlo'}->{'tPower'}
	if ( exists($config->{'monteCarlo'}->{'tPower'}) );
	print oHndl "sim->SetTpow(".$Tpow.")\n";
    }
    print oHndl "sim->SetLinearMapping(1)\n";
    print oHndl "sim->SetJumpFreq(10)\n";
    my $gamma = 0.2;
    $gamma = $config->{'monteCarlo'}->{'gamma'}
      if ( exists($config->{'monteCarlo'}->{'gamma'}) );
    print oHndl "set gamma = ".$gamma."\n";
    print oHndl "sim->NewGamma(gamma)\n";
}
print oHndl "sim->SetControl(1)\n";
print oHndl "sim->EnableLogging()\n";

# Register our likelihood function.
print oHndl "sim->SetUserLikelihood(fct)\n";

# Run the simulation.
print oHndl "set run = new RunOneSimulation(nsteps, 0, sstat,  prior, sim)\n";
print oHndl "run->Run()\n";

# Compute ensemble statistics.
print oHndl "sstat->ComputeDistribution()\n";

# Close the BIE input file.
close(oHndl);

# Remove old persistence store.
system("rm -rf pdir/".$config->{'name'});

# Create a PBS submit script.
open(oHndl,">".$config->{'workDirectory'}."/mcmc/bieLaunch.sh");
print oHndl "#!/bin/bash\n";
my $queue = "";
$queue = "#PBS -q ".$config->{'queue'}."\n"
    if ( exists($config->{'queue'}) );
print oHndl $queue;
print oHndl "#PBS -l nodes=".$nodeCount.":ppn=".$threadsPerNode."\n";
print oHndl "#PBS -o ".$config->{'workDirectory'}."/mcmc/bie.out\n";
print oHndl "#PBS -e ".$config->{'workDirectory'}."/mcmc/bie.err\n";
print oHndl "#PBS -l walltime=".$config->{'walltimeLimit'}."\n"
    if ( exists($config->{'walltimeLimit'}) );
print oHndl "#PBS -l mem=".$config->{'memoryLimit'}."\n"
    if ( exists($config->{'memoryLimit'}) );
print oHndl "#PBS -N galacticusBIE\n";
print oHndl "#PBS -V\n";
my $pwd = `pwd`;
print oHndl "cd ".$pwd;
# Set any custom environment.
if ( exists($config->{'environment'}) ) {
    my @environment;
    if ( UNIVERSAL::isa($config->{'environment'},"ARRAY") ) {
	@environment = @{$config->{'environment'}};
    } else {
	push(@environment,$config->{'environment'});
    }
    foreach ( @environment ) {
	print oHndl "export ".$_."\n";
    }
}
# Execute any custom commands.
if ( exists($config->{'execute'}) ) {
    my @execute;
    if ( UNIVERSAL::isa($config->{'execute'},"ARRAY") ) {
	@execute = @{$config->{'execute'}};
    } else {
	push(@execute,$config->{'execute'});
    }
    foreach ( @execute ) {
	print oHndl $_."\n";
    }
}
print oHndl "rm -rf pdir/".$config->{'name'}."\n";
print oHndl "rm -rf pdir/".$config->{'name'}."_restart*\n";
print oHndl "scratch=\$(echo ".$scratchDirectory." | sed -r s/'\\/'/'\\\\\\/'/g)\n";
print oHndl "sed -i~ -r s/\"\%\%SCRATCHDIRECTORY\%\%\"/\"\$scratch\"/g ".$config->{'workDirectory'}."/mcmc/bieInput.txt\n";
print oHndl "mpirun --mca btl ^openib --mca btl_tcp_if_include eth0 --mca mpi_yield_when_idle 1 -npernode ".$biePerNode." -hostfile \$PBS_NODEFILE ".$bie." -f ".$config->{'workDirectory'}."/mcmc/bieInput.txt\n";
print oHndl "mv EnsembleStat.log.0 ".$config->{'workDirectory'}."/mcmc/\n";
print oHndl "echo done > ".$config->{'workDirectory'}."/mcmc/done\n";
print oHndl "rm -rf ".$scratchDirectory."\n";
close(oHndl);

# Submit the job to PBS.
my $launchCommand = "qsub ".$config->{'workDirectory'}."/mcmc/bieLaunch.sh";
my $jobID         = `$launchCommand`;
chomp($jobID);
my $jobNumber = "";
if ( $jobID =~ m/^(\d+)/ ) {
    $jobNumber = $1;
}

# Create a BIE restart script.
open(oHndl,">".$config->{'workDirectory'}."/mcmc/bieInput_restart1.txt");
print oHndl "prestore ".$config->{'name'}."\n";
print oHndl "pnewsession ".$config->{'name'}."_restart1\n";
print oHndl "cktoggle\n";
$checkpointInterval = 1;
$checkpointInterval = $config->{'checkpoint'}->{'interval'} if ( exists($config->{'checkpoint'}->{'interval'}) );
print oHndl "ckinterval ".$checkpointInterval."\n";
print oHndl "run->Run()\n";
print oHndl "sstat->ComputeDistribution()\n";
close(oHndl);
# If autorestart is requested, submit a recursive job which will restart the simulation.
my $autoRestart = "false";
$autoRestart = $config->{'checkpoint'}->{'autoRestart'} if ( exists($config->{'checkpoint'}->{'autoRestart'}) );
if ( $autoRestart eq "true" ) {
    # Create a PBS submit script.
    open(oHndl,">".$config->{'workDirectory'}."/mcmc/bieLaunch_restart1.sh");
    print oHndl "#!/bin/bash\n";
    print oHndl $queue;
    print oHndl "#PBS -l nodes=".$config->{'nodes'}.":ppn=".$threadsPerNode."\n";
    print oHndl "#PBS -o ".$config->{'workDirectory'}."/mcmc/bie_restart1.out\n";
    print oHndl "#PBS -e ".$config->{'workDirectory'}."/mcmc/bie_restart1.err\n";
    print oHndl "#PBS -l walltime=".$config->{'walltimeLimit'}."\n"
	if ( exists($config->{'walltimeLimit'}) );
    print oHndl "#PBS -l mem=".$config->{'memoryLimit'}."\n"
	if ( exists($config->{'memoryLimit'}) );
    print oHndl "#PBS -N galacticusBIE\n";
    print oHndl "#PBS -V\n";
    print oHndl "#PBS -W depend=afternotok:".$jobNumber."\n";
    my $pwd = `pwd`;
    print oHndl "cd ".$pwd;
    # Set any custom environment.
    if ( exists($config->{'environment'}) ) {
	my @environment;
	if ( UNIVERSAL::isa($config->{'environment'},"ARRAY") ) {
	    @environment = @{$config->{'environment'}};
	} else {
	    push(@environment,$config->{'environment'});
	}
	foreach ( @environment ) {
	    print oHndl "export ".$_."\n";
	}
    }
    # Execute any custom commands.
    if ( exists($config->{'execute'}) ) {
	my @execute;
	if ( UNIVERSAL::isa($config->{'execute'},"ARRAY") ) {
	    @execute = @{$config->{'execute'}};
	} else {
	    push(@execute,$config->{'execute'});
	}
	foreach ( @execute ) {
	    print oHndl $_."\n";
	}
    }
    print oHndl "if [ -e ".$config->{'workDirectory'}."/mcmc/done ]; then\n";
    print oHndl " echo run is finished, will not restart\n";
    print oHndl "else\n";
    print oHndl " perl -pe 'my \$replace = sub { \"restart\".(\$_[0]+1) };s {restart(\\d)} {\$replace->(\$1)}ge' ".$config->{'workDirectory'}."/mcmc/bieLaunch_restart1.sh > ".$config->{'workDirectory'}."/mcmc/bieLaunch_restart2.sh\n";
    print oHndl " perl -pe 'my \$replace = sub { \"restart\".(\$_[0]+1) };s {restart(\\d)} {\$replace->(\$1)}ge' ".$config->{'workDirectory'}."/mcmc/bieInput_restart1.txt > ".$config->{'workDirectory'}."/mcmc/bieInput_restart2.txt\n";
    print oHndl " sed -i~ -r s/\"afternotok:[0-9]+\"/\"afternotok:\$PBS_JOBID\"/ ".$config->{'workDirectory'}."/mcmc/bieLaunch_restart2.sh\n";
    print oHndl "# qsub ".$config->{'workDirectory'}."/mcmc/bieLaunch_restart2.sh\n";
    print oHndl " mpirun --mca btl ^openib --mca btl_tcp_if_include eth0 -npernode ".$biePerNode." -hostfile \$PBS_NODEFILE ".$bie." -f ".$config->{'workDirectory'}."/mcmc/bieInput_restart1.txt\n";
    print oHndl " mv EnsembleStat.log.0 ".$config->{'workDirectory'}."/mcmc/\n";
    print oHndl " echo done > ".$config->{'workDirectory'}."/mcmc/done\n";
    print oHndl " rm -rf ".$scratchDirectory."\n";
    print oHndl "fi\n";
    close(oHndl);
    # Submit the job.
    system("qsub ".$config->{'workDirectory'}."/mcmc/bieLaunch_restart1.sh");
} else {
    # Automatric restart is not required. Create a script that can be used to restart manually.
    open(oHndl,">".$config->{'workDirectory'}."/mcmc/bieLaunch_restart1.sh");
    print oHndl "#!/bin/bash\n";
    print oHndl $queue;
    print oHndl "#PBS -l nodes=".$config->{'nodes'}.":ppn=".$threadsPerNode."\n";
    print oHndl "#PBS -o ".$config->{'workDirectory'}."/mcmc/bie_restart1.out\n";
    print oHndl "#PBS -e ".$config->{'workDirectory'}."/mcmc/bie_restart1.err\n";
    print oHndl "#PBS -l walltime=".$config->{'walltimeLimit'}."\n"
	if ( exists($config->{'walltimeLimit'}) );
    print oHndl "#PBS -l mem=".$config->{'memoryLimit'}."\n"
	if ( exists($config->{'memoryLimit'}) );
    print oHndl "#PBS -N galacticusBIE\n";
    print oHndl "#PBS -V\n";
    my $pwd = `pwd`;
    print oHndl "cd ".$pwd;
    # Set any custom environment.
    if ( exists($config->{'environment'}) ) {
	my @environment;
	if ( UNIVERSAL::isa($config->{'environment'},"ARRAY") ) {
	    @environment = @{$config->{'environment'}};
	} else {
	    push(@environment,$config->{'environment'});
	}
	foreach ( @environment ) {
	    print oHndl "export ".$_."\n";
    }
    }
    # Execute any custom commands.
    if ( exists($config->{'execute'}) ) {
	my @execute;
	if ( UNIVERSAL::isa($config->{'execute'},"ARRAY") ) {
	    @execute = @{$config->{'execute'}};
	} else {
	    push(@execute,$config->{'execute'});
	}
    foreach ( @execute ) {
	print oHndl $_."\n";
    }
    }
    print oHndl "mpirun --mca btl ^openib --mca btl_tcp_if_include eth0 --mca mpi_yield_when_idle 1 -npernode ".$biePerNode." -hostfile \$PBS_NODEFILE ".$bie." -f ".$config->{'workDirectory'}."/mcmc/bieInput_restart1.txt\n";
    print oHndl "mv EnsembleStat.log.0 ".$config->{'workDirectory'}."/mcmc/\n";
    print oHndl "echo done > ".$config->{'workDirectory'}."/mcmc/done\n";
    print oHndl "rm -rf ".$scratchDirectory."\n";
    close(oHndl);
}

exit;
