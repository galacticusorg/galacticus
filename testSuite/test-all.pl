#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use threads;
use Date::Format;
use XML::Simple;
use MIME::Lite;
use Net::SMTP::SSL;
use Data::Dumper;
use File::Slurp qw( slurp );
use File::Find;
use Term::ReadKey;
use System::Redirect;
use Galacticus::Launch::PBS;
use Galacticus::Options;

# Run a suite of tests on the Galacticus code.
# Andrew Benson (19-Aug-2010).

# Read in any configuration options.
my $config = &Galacticus::Options::LoadConfig();

my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     })
    if ( defined($queueManager) );

# Set default options.
my %options =
    (
     'from-scratch'        => "yes",
     'skip-slow'           => "no" ,
     'timing'              => "no" ,
     'build-documentation' => "yes",
     'pbsJobMaximum'       => (defined($queueConfig) && exists($queueConfig->{'jobMaximum'})) ? $queueConfig->{'jobMaximum'} : 100,
     'processesPerNode'    => (defined($queueConfig) && exists($queueConfig->{'ppn'       })) ? $queueConfig->{'ppn'       } :   1,
     'submitSleepDuration' =>                                                                                                    1,
     'waitSleepDuration'   =>                                                                                                   10
    );

# Get any command line options.
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Identify e-mail options for this host.
my $emailConfig;
my $smtpPassword;
if ( exists($config->{'email'}->{'host'}->{$ENV{'HOSTNAME'}}) ) {
    $emailConfig = $config->{'email'}->{'host'}->{$ENV{'HOSTNAME'}};
} elsif ( exists($config->{'email'}->{'host'}->{'default'}) ) {
    $emailConfig = $config->{'email'}->{'host'}->{'default'};
} else {
    $emailConfig->{'method'} = "sendmail";
}
if ( $emailConfig->{'method'} eq "smtp" && exists($emailConfig->{'passwordFrom'}) ) {
    # Get any password now.
    if ( $emailConfig->{'passwordFrom'} eq "input" ) {
	print "Please enter your e-mail SMTP password:\n";
	$smtpPassword = &getPassword;
    }
    elsif ( $emailConfig->{'passwordFrom'} eq "kdewallet" ) {
	my $appName          = "Galacticus";
	my $folderName       = "glc-test-all";
	require Net::DBus;
	my $bus           = Net::DBus->find;
	my $walletService = $bus->get_service("org.kde.kwalletd");
	my $walletObject  = $walletService->get_object("/modules/kwalletd");
	my $walletID      = $walletObject->open("kdewallet",0,$appName);
	if ( $walletObject->hasEntry($walletID,$folderName,"smtpPassword",$appName) == 1 ) {
	    $smtpPassword = $walletObject->readPassword($walletID,$folderName,"smtpPassword",$appName); 
	} else {
	    print "Please enter your e-mail SMTP password:\n";
	    $smtpPassword = &getPassword;
	    $walletObject->writePassword($walletID,$folderName,"smtpPassword",$smtpPassword,$appName); 
	}
    }
}

# Read any existing timing data.
my $timing;
if ( -e "testSuite/timing.xml" ) {
    my $xml = new XML::Simple;
    $timing = $xml->XMLin("testSuite/timing.xml");
}

# Open a log file.
my $logFile = "testSuite/allTests.log";
open(lHndl,">".$logFile);

# Create a directory for test suite outputs.
if ( exists($options{'outputPath'}) ) {
    system("rm -rf "  .$options{'outputPath'}." testSuite/outputs");
    system("mkdir -p ".$options{'outputPath'});
    system("cd testSuite; ln -sf ".$options{'outputPath'}." outputs");
} else {
    system("rm -rf testSuite/outputs");
    system("mkdir -p testSuite/outputs");
}

# Write header to log file.
print lHndl ":-> Running test suite:\n";
print lHndl "    -> Host:\t".$ENV{'HOSTNAME'}."\n";
print lHndl "    -> Time:\t".time2str("%a %b %e %T (%Z) %Y", time)."\n";

# Stack to be used for PBS jobs.
my @jobStack;

# Determine if slow tests are to be skipped.
my $skipSlow = $options{'skip-slow'} eq "yes";

# Define a list of executables to run. Each hash must give the name of the executable, should specify whether or not the
# executable should be run inside of Valgrind (this is useful for detecting errors which lead to misuse of memory but which don't
# necessary cause a crash), and should specify if the executable should be built and run under mpi.
my @executablesToRun = (
    {
	name     => "tests.nodes.exe",                                                    # Tests of Galacticus nodes.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.parameters.exe",                                               # Tests of parameter input.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.files.exe",                                                    # Tests of file functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.IO.HDF5.exe",                                                  # Tests of HDF5 IO routines.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.IO.XML.exe",                                                   # Tests of XML IO routines.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.ODE_solver.exe",                                               # Tests of ODE solver routines.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.random.exe",                                                   # Tests of random number generators.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.random.quasi.exe",                                             # Tests of quasi-random number generators.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.arrays.exe",                                                   # Tests of array functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.meshes.exe",                                                   # Tests of mesh functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.convex_hulls.exe",                                             # Tests of convex hull functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.comparisons.exe",                                              # Tests of comparison functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.geometry.coordinate_systems.exe",                              # Tests of coordinate system functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.hashes.exe",                                                   # Tests of hashing utilities.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.hashes.perfect.exe",                                           # Tests of perfect hashing utilities.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.regular_expressions.exe",                                      # Tests of regular expression utilities.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.hashes.cryptographic.exe",                                     # Tests of cryptographic hashing utilities.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.integration.exe",                                              # Tests of integration functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.integration2.exe",                                             # Tests of integration functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.differentiation.exe",                                          # Tests of differentiation functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.tables.exe",                                                   # Tests of table functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.interpolation.exe",                                            # Tests of interpolation functions.
	valgrind => 0,
	mpi      => 0             
    },
    {
	name     => "tests.interpolation.2D.exe",                                         # Tests of 2D interpolation function.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.make_ranges.exe",                                              # Tests of numerical range building functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.mass_distributions.exe",                                       # Tests of mass distributions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.math_special_functions.exe",                                   # Tests of mathematical special functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.math_distributions.exe",                                       # Tests of mathematical distributions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.math.fast.exe",                                                # Tests of fast mathematical functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.math.linear_algebra.exe",                                      # Tests of linear algebra functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.root_finding.exe",                                             # Tests of root finding functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.multi_dimensional_minimizer.exe",                              # Tests of multidimensional minimization functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.search.exe",                                                   # Tests of searching functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.Roman_numerals.exe",                                           # Tests of Roman numeral conversion functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.sort.exe",                                                     # Tests of sorting functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.sort.topological.exe",                                         # Tests of topological sorting functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.string_utilities.exe",                                         # Tests of string handling utilities.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.vectors.exe",                                                  # Tests of vector functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.multi_counters.exe",                                           # Tests of multi-counters.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.tensors.exe",                                                  # Tests of tensor functions.
	valgrind => 0,
	mpi      => 0
    },
    {	name     => "tests.cosmic_age.exe",                                               # Tests of cosmic age calculations.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.open.exe",                                  # Tests of spherical collapse calculations.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.flat.exe",                                  # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.EdS.exe",                       # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.open.exe",                      # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.lambda.exe",                    # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminusHalf.exe",      # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminusTwoThirds.exe", # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminus0.6.exe",       # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminus0.8.exe",       # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spherical_collapse.baryons_dark_matter.exe",                   # .
	valgrind => 0,
	mpi      => 0,
	isSlow   => 1
    },
    {
	name     => "tests.spherical_collapse.nonlinear.exe",                             # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.warm_dark_matter.exe",                                         # Tests of critical overdensity for collapse in warm dark matter model.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.linear_growth.cosmological_constant.exe",                      # Tests of linear growth factor.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.linear_growth.EdS.exe",                                        # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.linear_growth.open.exe",                                       # .
 	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.linear_growth.dark_energy.exe",                                # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.linear_growth.baryons.EdS.exe",                                # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.halo_mass_function.Tinker.exe",                                # Tests of dark matter halo mass functions.
 	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.halo_mass_function.Reed2007.exe",                                # Tests of dark matter halo mass functions.
 	valgrind => 0,
	mpi      => 0
    },
    {
	name     =>"tests.comoving_distance.exe",                                         # Tests of comoving distance calculations.
 	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.mass_accretion_history.Correa2015.exe",                        # Tests of mass accretion histories.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.Zhao2009_algorithms.dark_energy.exe",                          # Tests of Zhao et al. (2009) algorithms.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.Zhao2009_algorithms.EdS.exe",                                  # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.Zhao2009_algorithms.open.exe",                                 # .
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.NFW96_concentration.dark_energy.exe",                          # Tests of Navarro, Frenk & White (1996) halo concentration algorithm.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.Prada2011_concentration.exe",                                  # Tests of Prada et al. (2011) halo concentration algorithm.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.DiemerKravtsov2014_concentration.exe",                         # Tests of Diemer & Kravtsov (2014) halo concentration algorithm.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.concentration.Correa2015.exe",                                 # Tests of Correa et al. (2015) halo concentration algorithm.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.concentrations.exe",                                           # Tests of various halo concentration algorithms.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.biases.exe",                                                   # Tests of various halo bias algorithms.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.kepler_orbits.exe",                                            # Keplerian orbital parameter conversions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.abundances.exe",                                               # Abundances objects.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.sigma.exe",                                                    # Sigma(M).
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.power_spectrum.exe",                                           # Power spectrum.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.power_spectrum.primordial.exe",                                # Primordial power spectrum.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.transfer_functions.exe",                                       # Transfer functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.black_hole_fundamentals.exe",                                  # Black hole fundamentals.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.bug745815.exe",                                                # Regressions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.tree_branch_destroy.exe",                                      # Tests of merger tree walking.
	valgrind => 1,
	valgrindOptions => "--undef-value-errors=no",
	mpi      => 0
    },
    {
	name     => "tests.gaunt_factors.exe",                                            # Tests of Gaunt factors.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.cooling_functions.exe",                                        # Tests of cooling functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.accretion_disks.exe",                                          # Tests of accretion disks.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.mass_accretion_history.Hearin2021.exe",                        # Tests of dark matter halo mass accretion histories.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.mass_accretion_history.Hearin2021_stochastic.exe",             # Tests of stochastic dark matter halo mass accretion histories.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.dark_matter_profiles.exe",                                     # Tests of dark matter profiles.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.dark_matter_profiles.heated.exe",                              # Tests of heated dark matter profiles.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.dark_matter_halo_radius_enclosing_mass.exe",                   # Tests of dark matter halo radius enclosing mass functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.dark_matter_profiles.tidal_tracks.exe",                        # Tests of dark matter profile tidal tracks.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.dark_matter_profiles.generic.exe",                             # Tests of generic, numerical implementations of dark matter profile functions.
	valgrind => 0,
	mpi      => 0,
	isSlow   => 1
    },
    {
	name     => "tests.dark_matter_profiles.finite_resolution.exe",                   # Tests of finite resolution dark matter profiles.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.spectra.postprocess.Inoue2014.exe",                            # Tests of Inoue (2014) IGM attenutation model.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.MPI.exe",                                                      # Tests of MPI functionality.
	valgrind => 0,
	mpi      => 4
    },
    {
	name     => "tests.radiative_transfer.atomic_matter.state_solver.exe",            # Tests of radiative transfer atomic matter state solving.
	valgrind => 0,
	mpi      => 1
    },
    {
	name     => "tests.recombination_computed.exe",                                   # Tests of computed radiative recombination coefficients.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.recombination_cooling.Hummer.exe",                             # Tests of Hummer radiative recombination coefficients.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.locks.exe",                                                    # Tests of OpenMP locking functionality.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.initial_mass_functions.exe",                                   # Tests of initial mass functions.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.stellar_populations.exe",                                      # Tests of stellar populations.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.stellar_populations.luminosities.exe",                         # Tests of stellar population luminosities.
	valgrind => 0,
	mpi      => 0
    },
    {
	name     => "tests.crash.exe",                                                    # Tests that crashes are detected.
	valgrind => 0,
	mpi      => 0,
	expect   => "crash"
    },
    {
	name     => "tests.fail.exe",                                                     # Tests that failures are detected.
	valgrind => 0,
	mpi      => 0,
	expect   => "fail"
    },
    {
	name     => "tests.excursion_sets.exe",                                           # Tests of excursion set solvers.
	valgrind =>  0,
	ppn      => $options{'processesPerNode'},
	mpi      =>  0,
	isSlow   =>  1
    },
    {
	name     => "tests.merger_tree_branching.exe",                                    # Tests of merger tree branching rate functions.
	valgrind => 0,
	mpi      => 0,
	isSlow   => 1
    },
    {
	name     => "tests.event_hooks.exe",                                              # Tests of event hook infrastructure.
	valgrind => 0,
	mpi      => 0
    }
    );

# Build all executables.
@jobStack = ();
my @executablesNonMPI  = map {($_->{'mpi'} == 0) && (! $skipSlow || ! exists($_->{'isSlow'}) || $_->{'isSlow'} == 0) ? $_->{'name'} : ()} @executablesToRun;
my @executablesMPI     = map {($_->{'mpi'} >  0) && (! $skipSlow || ! exists($_->{'isSlow'}) || $_->{'isSlow'} == 0) ? $_->{'name'} : ()} @executablesToRun;
{
    my $compileCommand;
    $compileCommand .= "rm -rf ./work/build\n"
	if ( $options{'from-scratch'} eq "yes" );
    $compileCommand .= "make -j".$options{'processesPerNode'}." GALACTICUS_BUILD_DOCS=yes "             .join(" ",@executablesNonMPI)." Galacticus.exe\n";
    my %testBuildJob =
	(
	 launchFile   => "testSuite/compileCodeNonMPI.pbs",
	 label        => "testSuite-compileCode-nonMPI"   ,
	 logFile      => "testSuite/compileCodeNonMPI.log",
	 command      =>  $compileCommand                 ,
	 ppn          =>  $options{'processesPerNode'}    ,
	 tracejob     => "yes"                            ,
	 onCompletion => 
	 {
	     function  => \&testCompileFailure,
	     arguments => [ "testSuite/compileCodeNonMPI.log", "Test code (non-MPI) compilation" ]
	 }
	);
    push(@jobStack,\%testBuildJob);
}
{
    my $compileCommand;
    $compileCommand .= "rm -rf ./work/buildMPI\n"
	if ( $options{'from-scratch'} eq "yes" );
    $compileCommand .= "make -j".$options{'processesPerNode'}." GALACTICUS_BUILD_OPTION=MPI "           .join(" ",@executablesMPI   ).               "\n"
	if ( scalar(@executablesMPI) > 0 );
    $compileCommand .= "make -j".$options{'processesPerNode'}." GALACTICUS_BUILD_OPTION=MPI SUFFIX=_MPI".     " "                    ." Galacticus.exe\n";
    my %testBuildJob =
	(
	 launchFile   => "testSuite/compileCodeMPI.pbs"   ,
	 label        => "testSuite-compileCode-MPI"      ,
	 logFile      => "testSuite/compileCodeMPI.log"   ,
	 command      =>  $compileCommand                 ,
	 ppn          =>  $options{'processesPerNode'}    ,
	 tracejob     => "yes"                            ,
	 onCompletion => 
	 {
	     function  => \&testCompileFailure,
	     arguments => [ "testSuite/compileCodeMPI.log", "Test code (MPI) compilation"       ]
	 }
	);
    push(@jobStack,\%testBuildJob);
}
&Galacticus::Launch::PBS::SubmitJobs(\%options,@jobStack);
unlink("testSuite/compileCodeNonMPI.pbs","testSuite/compileCodeMPI.pbs");

# Create a job queue.
my @launchFiles;
@jobStack = ();

# Queue documentation build.
if ( $options{'build-documentation'} eq "yes" ) {
    my %job =
	(
	 launchFile   => "testSuite/docBuild.pbs",
	 label        => "testSuite-docBuild"    ,
	 logFile      => "testSuite/docBuild.log",
	 ppn          => $options{'processesPerNode'},
	 command      => "./scripts/doc/Build_Documentation.sh -p ".$options{'processesPerNode'}." -f no",
	 tracejob     => "yes"                   ,
	 onCompletion => 
	 {
	     function  => \&testDocBuild,
	     arguments => [ "testSuite/docBuild.log" ]
	 }
	);
    push(@jobStack,\%job);
}

# Queue all executables.
foreach my $executable ( @executablesToRun ) {
    next
	if ( $skipSlow && exists($executable->{'isSlow'}) && $executable->{'isSlow'} == 1 );
    # Generate the job.
    if ( exists($executable->{'name'}) ) {
	(my $label = $executable->{'name'}) =~ s/\./_/;
	my $ppn = 1;
	$ppn = $executable->{'ppn'}
	    if ( exists($executable->{'ppn'}) );
	$ppn = $executable->{'mpi'}
	    if ( exists($executable->{'mpi'}) && $executable->{'mpi'} > 0 );
	my $launchFile = "testSuite/".$label.".pbs";
	push(@launchFiles,$launchFile);
	$executable->{'expect'} = "success"
	    unless ( exists($executable->{'expect'}) );
	my %job =
	    (
	     launchFile   => $launchFile               ,
	     label        => "testSuite-".$label       ,
	     logFile      => "testSuite/".$label.".log",
	     ppn          => $ppn                      ,
	     tracejob     => "yes"                     ,
	     onCompletion => 
	     {
		 function  => \&testFailure,
		 arguments => [ "testSuite/".$label.".log", "Test code: ".$executable->{'name'}, $executable->{'expect'} ]
	     }
	    );
	if ( $executable->{'valgrind'} == 1 ) {	    
	    $job{'command'} = "valgrind --error-exitcode=1 ".$executable->{'valgrindOptions'}." ".$executable->{'name'};
	} elsif ( $executable->{'mpi'} > 0 ) {
	    $job{'command'} = "mpirun -np "                 .$executable->{'mpi'            }." ".$executable->{'name'};
	} else {
	    $job{'command'} =                                                                     $executable->{'name'};
	}
	push(@jobStack,\%job);
    }
}

# Perform tests utilizing Galacticus itself, both without and with MPI.
our $mpi;
my  @launchLocal;
foreach $mpi ( "noMPI", "MPI" ) {
    my $suffix = $mpi eq "MPI" ? "_MPI" : "";
    if ( -e "./Galacticus.exe".$suffix ) {
	# Find all test scripts to run.
	my @testDirs = ( "testSuite/" );
	find(\&runTestScript,@testDirs);
    }
}

# Run scripts that can launch themselves using PBS.
my $thread = threads->create(\&launchLocalTests, \@launchLocal, $options{'pbsJobMaximum'}, $options{'submitSleepDuration'}, $options{'waitSleepDuration'});

# Launch all PBS job tests.
&Galacticus::Launch::PBS::SubmitJobs(\%options,@jobStack);
unlink(@launchFiles);

# Wait for local jobs to complete.
$thread->join();
# Scan for any job meta-data from local jobs.
open(my $localJobLog,"testSuite/allTests.tmp");
while ( my $localLine = <$localJobLog> ) {
    # Append local log to the main log.
    print lHndl $localLine;
    # Detect meta-data.
    if ( $localLine =~ m/\-> Job metadata:.*=/ ) {
	chomp($localLine);
	$localLine     =~ s/^.*=\s*//;
	my @data       =  split(/\s+:\s+/,$localLine);
	(my $jobPath   =  $data[0]) =~ s/[\/:]/_/g;
	my $exitStatus =  $data[1];
	my $startAt    =  $data[2];
	my $endAt      =  $data[3];
	$timing->{$jobPath} =
	{
	    startTime => $startAt,
	    endTime   => $endAt
	}
	if ( $exitStatus == 0 );
    }
}
close($localJobLog);
# Remove the local job output.
unlink("testSuite/allTests.tmp");

# Close the log file.
close(lHndl);

# Scan the log file for FAILED.
my $lineNumber = 0;
my @failLines;
open(lHndl,$logFile);
while ( my $line = <lHndl> ) {
    ++$lineNumber;
    if ( $line =~ m/FAILED/ ) {
	push(@failLines,$lineNumber);
    }
    if ( $line =~ m/SKIPPED/ ) {
	push(@failLines,$lineNumber);
    }
}
close(lHndl);
open(lHndl,">>".$logFile);
my $emailSubject = "Galacticus test suite log";
my $exitStatus;
if ( scalar(@failLines) == 0 ) {
    print lHndl "\n\n:-> All tests were successful.\n";
    print       "All tests were successful.\n";
    $emailSubject .= " [success]";
    $exitStatus = 0;
} else {
    print lHndl "\n\n:-> Failures found. See following lines in log file:\n\t".join("\n\t",@failLines)."\n";
    print "Failure(s) found - see ".$logFile." for details.\n";
    $emailSubject .= " [FAILURE]";
    $exitStatus = 1;
}
close(lHndl);

# Output timing data.
if ( defined($timing) && $options{'timing'} eq "yes" ) {
    my $xml = new XML::Simple;
    open(my $timingFile,">testSuite/timing.xml");
    print $timingFile $xml->XMLout($timing,RootName => "timing");
    close($timingFile);
}

# If we have an e-mail address to send the log to, then do so.
if ( defined($config->{'contact'}->{'email'}) ) {
    if ( $config->{'contact'}->{'email'} =~ m/\@/ ) {
	# Get e-mail configuration.
	my $sendMethod = $emailConfig->{'method'};
	# Construct the message.
	my $message  = "Galacticus test suite log is attached.\n";
	my $msg = MIME::Lite->new(
	    From    => '',
	    To      => $config->{'contact'}->{'email'},
	    Subject => $emailSubject,
	    Type    => 'TEXT',
	    Data    => $message
	    );
	system("bzip2 -f ".$logFile);
	$msg->attach(
	    Type     => "application/x-bzip",
	    Path     => $logFile.".bz2",
	    Filename => "allTests.log.bz2"
	    );
	if ( $sendMethod eq "sendmail" ) {
	    $msg->send;
	}
	elsif ( $sendMethod eq "smtp" ) {
	    my $smtp; 
	    $smtp = Net::SMTP::SSL->new($config->{'email'}->{'host'}, Port=>465) or die "Can't connect";
	    $smtp->auth($config->{'email'}->{'user'},$smtpPassword) or die "Can't authenticate:".$smtp->message();
	    $smtp->mail( $config->{'contact'}->{'email'}) or die "Error:".$smtp->message();
	    $smtp->to( $config->{'contact'}->{'email'}) or die "Error:".$smtp->message();
	    $smtp->data() or die "Error:".$smtp->message();
	    $smtp->datasend($msg->as_string) or die "Error:".$smtp->message();
	    $smtp->dataend() or die "Error:".$smtp->message();
	    $smtp->quit() or die "Error:".$smtp->message();
	}
    }
}

exit $exitStatus;

sub runTestScript {
    # Run a test script.
    my $fileName = $_;
    chomp($fileName);

    # Test if this is a script to run.
    if ( $fileName =~ m/^test\-.*\.pl$/ && $fileName ne "test-all.pl" ) {
	if ( ( $mpi eq "MPI" && $fileName =~ m/_MPI\.pl/ ) || ( $mpi eq "noMPI" && $fileName !~ m/_MPI\.pl/ ) ) {
	    system("grep -q -e launch.pl -e 'selfManage: true' ".$fileName);
	    if ( $? == 0 ) {
		# This script will launch its own models.
		push(
		    @launchLocal,
		    $fileName
		    );
	    } else {
		# We need to launch this script.
		(my $label = $fileName) =~ s/\.pl$//;
		push(
		    @jobStack,
		    {
			launchFile   => "testSuite/".$label.".pbs"  ,
			label        => "testSuite-".$label         ,
			logFile      => "testSuite/".$label.".log"  ,
			command      => "cd testSuite; ".$fileName." ".join(" ",map {"--".$_." ".$options{$_}} keys(%options)),
			ppn          => $options{'processesPerNode'},
			tracejob     => "yes"                       ,
			onCompletion => 
			{
			    function  => \&testFailure,
			    arguments => [ "testSuite/".$label.".log", "Test script '".$label."'", "success" ]
			}
		    }
		    );
	    }
	}
    }
}

sub getPassword {
    # Read a password from standard input while echoing asterisks to the screen.
    ReadMode('noecho');
    ReadMode('raw');
    my $password = '';
    while (1) {
	my $c;
	1 until defined($c = ReadKey(-1));
	last if $c eq "\n";
	print "*";
	$password .= $c;
    }
    ReadMode('restore');
    print "\n";
    return $password;
}

sub testFailure {
    # Callback function which checks for failure of jobs run in PBS.
    my $logFile    = shift();
    my $jobMessage = shift();
    my $expect     = shift();
    my $jobID      = shift();
    my $exitStatus = shift();
    shift();
    my $job        = shift();
    my $startTime  = shift();
    my $endTime    = shift();
    my $result;
    # Branch on expected behavior.
    if ( $expect eq "crash" ) {
	if ( $exitStatus == 0 ) {
	    # We expected a crash, but didn't get one. This is a failure.
	    $result = "FAILED: ".$jobMessage." (crash expected)\n";
	} else {
	    # We expected a crash and we got one. This is a success.
	    $result = "SUCCESS: ".$jobMessage."\n";
	}
    } else {
	# Check for crash.
	if ( $exitStatus == 0 ) {
	    # Check for failure message in log file.
	    system("grep -q FAIL ".$logFile);
	    if ( $? == 0 ) {
		# Failure detected.
		if ( $expect eq "fail" ) {
		    # We expected failure, and we did fail. This is a success.
		    $result = "SUCCESS: ".$jobMessage."\n";
		} elsif ( $expect eq "success" ) {
		    # We expected failure, but we did not fail. This is a failure.
		    $result = "FAILED: ".$jobMessage."\n";
		} else {
		    # Unknown expectation.
		    $result = "FAILED: ".$jobMessage." (unknown expectation)\n";
		}
	    } else {
		# No failure detected.
		if ( $expect eq "fail" ) {
		    # We expected failure, but we did not fail. This is a failure.
		    $result = "FAILED: ".$jobMessage." (failure expected)\n";
		} elsif ( $expect eq "success" ) {
		    # We expected failure, but we did not fail. This is a failure.
		    $result = "SUCCESS: ".$jobMessage."\n";
		} else {
		    # Unknown expectation.
		    $result = "FAILED: ".$jobMessage." (unknown expectation)\n";
		}
	    }
	} else {
	    # We did not expect a crash, but got one. This is a failure.
	    $result = "FAILED: ".$jobMessage." (crashed unexpectedly)\n";
	}
    }
    # Add metadata to the timing data file.
    if ( $result =~ m/SUCCESS/ ) {
	$timing->{$job->{'label'}} =
	{
	    startTime => $startTime,
	    endTime   => $endTime
	};
    }
    # Report success or failure.
    print lHndl $result;
    if ( $result =~ m/^FAILED:/ ) {
	# Job failed.
	print lHndl "Job output follows:\n";
	print lHndl slurp($logFile);
    }
}

sub testCompileFailure {
    # Callback function which checks for failure of compile jobs run in PBS.
    my $logFile     = shift();
    my $jobMessage  = shift();
    my $jobID       = shift();
    my $errorStatus = shift();
    shift();
    my $job        = shift();
    my $startTime  = shift();
    my $endTime    = shift();
    # Check for failure message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q FAIL ".$logFile);
	$errorStatus = 1
	    if ( $? == 0 );	
    }
    # Check for compiler error message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q Error: ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = "Compiler errors issued\n".$jobMessage;
	}
    }
    # Check for compiler warning message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q Warning: ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = "Compiler warnings issued\n".$jobMessage;
	}
    }
    # Check for make error message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -v 'is up to date' ".$logFile." | grep -q make: ");
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = "Make errors issued\n".$jobMessage;
	}
    }
    # Report success or failure.
    if ( $errorStatus == 0 ) {
	# Job succeeded.
	print lHndl "SUCCESS: ".$jobMessage."\n";
	unlink($logFile); 
	# Add metadata to the timing data file.
	$timing->{$job->{'label'}} =
	{
	    startTime => $startTime,
	    endTime   => $endTime
	};
    } else {
	# Job failed.
	print lHndl "FAILED: ".$jobMessage."\n";
	print lHndl "Job output follows:\n";
	print lHndl slurp($logFile);
    }
}

sub testDocBuild {
    # Callback function which checks for failure of documentation build jobs run in PBS.
    my $logFile     = shift();
    my $jobID       = shift();
    my $errorStatus = shift();
    shift();
    my $job        = shift();
    my $startTime  = shift();
    my $endTime    = shift();
    # Check for failure message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q FAIL ".$logFile);
	$errorStatus = 1
	    if ( $? == 0 );	
    }
    my $jobMessage = ($errorStatus == 0 ? "SUCCESS" : "FAILED").": documentation build\n";
    # Check for compiler error message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q \"Warning: missing method descriptions in class\" ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = $jobMessage."FAILED [undefined method]\n";
	}
	system("grep -q -e \"pdfTeX warning\" -e \"referenced but does not exist\" ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = $jobMessage."FAILED [undefined hyperreferences]\n";
	}
	system("grep -q \"There were multiply-defined labels\" ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = $jobMessage."FAILED [multiply-defined labels]\n";
	}
	system("grep -q \"There were undefined references\" ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = $jobMessage."FAILED [undefined labels]\n";
	}
    }
    # Add metadata to the timing data file.
    if ( $errorStatus == 0 ) {
	$timing->{$job->{'label'}} =
	{
	    startTime => $startTime,
	    endTime   => $endTime
	};
    }
    # Report success or failure.
    print lHndl $jobMessage;
    if ( $errorStatus == 0 ) {
	# Job succeeded.
	unlink($logFile);
    } else {
	# Job failed.
	print lHndl "Job output follows:\n";
	print lHndl slurp($logFile);
    }
}

sub launchLocalTests {
    # Run scripts that can launch themselves using PBS.
    my @launchLocal         = @{shift()};
    my $pbsJobMaximum       =   shift() ;
    my $submitSleepDuration =   shift() ;
    my $waitSleepDuration   =   shift() ;
    print           ":-> Running test scripts:\n";
    print lHndl "\n\n:-> Running test scripts:\n";
    print       join("\n",map {"\t".$_} @launchLocal)."\n";
    print lHndl join("\n",map {"\t".$_} @launchLocal)."\n";
    open(my $script,">testSuite/outputs/launchLocal.sh");
    print $script "cd testSuite\n";
    foreach my $localScript ( @launchLocal ) {
	print $script $localScript." --pbsJobMaximum ".$pbsJobMaximum." --submitSleepDuration ".$submitSleepDuration." --waitSleepDuration ".$waitSleepDuration." &\n";
    }
    print $script "wait\n";
    print $script "exit\n";
    close($script);
    &System::Redirect::tofile("chmod u=wrx testSuite/outputs/launchLocal.sh; testSuite/outputs/launchLocal.sh","testSuite/allTests.tmp");
}
