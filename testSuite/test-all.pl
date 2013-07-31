#!/usr/bin/env perl
use strict;
use warnings;
use lib "./perl";
use System::Redirect;
use Date::Format;
use XML::Simple;
use MIME::Lite;
use Net::SMTP::SSL;
use Data::Dumper;
use File::Slurp qw( slurp );
use File::Find;
use Switch;
use Term::ReadKey;

# Run a suite of tests on the Galacticus code.
# Andrew Benson (19-Aug-2010).

# Read in any configuration options.
my $config;
if ( -e "galacticusConfig.xml" ) {
    my $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Identify e-mail options for this host.
my $emailConfig;
my $smtpPassword;
if ( exists($config->{'email'}->{'host'}->{$ENV{'HOST'}}) ) {
    $emailConfig = $config->{'email'}->{'host'}->{$ENV{'HOST'}};
} elsif ( exists($config->{'email'}->{'host'}->{'default'}) ) {
    $emailConfig = $config->{'email'}->{'host'}->{'default'};
} else {
    $emailConfig->{'method'} = "sendmail";
}
if ( $emailConfig->{'method'} eq "smtp" && exists($emailConfig->{'passwordFrom'}) ) {
    # Get any password now.
    switch ( $emailConfig->{'passwordFrom'} ) {
	case ( "input" ) {
	    print "Please enter your e-mail SMTP password:\n";
	    $smtpPassword = &getPassword;
	}
	case ( "kdewallet" ) {
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
}

# Open a log file.
my $logFile = "testSuite/allTests.log";
open(lHndl,">".$logFile);

# Clean up previous build.
system("rm -rf work/build/*");

# Create a directory for test suite outputs.
system("rm -rf testSuite/outputs");
system("mkdir -p testSuite/outputs");

# Write header to log file.
print lHndl ":-> Running test suite:\n";
print lHndl "    -> Host:\t".$ENV{'HOST'}."\n";
print lHndl "    -> Time:\t".time2str("%a %b %e %T (%Z) %Y", time)."\n";

# Define a list of executables to run. Each hash must give the name of the executable and should specify whether or not the
# executable should be run inside of Valgrind (this is useful for detecting errors which lead to misuse of memory but which don't
# necessary cause a crash).
my @executablesToRun = (
    {
	name     => "tests.nodes.exe",                                                    # Tests of Galacticus nodes.
	valgrind => 0
    },
    {
	name     => "tests.IO.HDF5.exe",                                                  # Tests of HDF5 IO routines.
	valgrind => 0
    },
    {
	name     => "tests.ODE_solver.exe",                                               # Tests of ODE solver routines.
	valgrind => 0
    },
    {
	name     => "tests.arrays.exe",                                                   # Tests of array functions.
	valgrind => 0
    },
    {
	name     => "tests.meshes.exe",                                                   # Tests of mesh functions.
	valgrind => 0
    },
    {
	name     => "tests.comparisons.exe",                                              # Tests of comparison functions.
	valgrind => 0
    },
    {
	name     => "tests.geometry.coordinate_systems.exe",                              # Tests of coordinate system functions.
	valgrind => 0
    },
    {
	name     => "tests.hashes.exe",                                                   # Tests of hashing utilities.
	valgrind => 0
    },
    {
	name     => "tests.hashes.perfect.exe",                                           # Tests of perfect hashing utilities.
	valgrind => 0
    },
    {
	name     => "tests.hashes.cryptographic.exe",                                     # Tests of cryptographic hashing utilities.
	valgrind => 0
    },
    {
	name     => "tests.integration.exe",                                              # Tests of integration functions.
	valgrind => 0
    },
    {
	name     => "tests.tables.exe",                                                   # Tests of table functions.
	valgrind => 0
    },
    {
	name     => "tests.interpolation.exe",                                            # Tests of interpolation functions.
	valgrind => 0             
    },
    {
	name     => "tests.interpolation.2D.exe",                                         # Tests of 2D interpolation function.
	valgrind => 0
    },
    {
	name     => "tests.make_ranges.exe",                                              # Tests of numerical range building functions.
	valgrind => 0
    },
    {
	name     => "tests.mass_distributions.exe",                                       # Tests of mass distributions.
	valgrind => 0
    },
    {
	name     => "tests.math_special_functions.exe",                                   # Tests of mathematical special functions.
	valgrind => 0
    },
    {
	name     => "tests.root_finding.exe",                                             # Tests of root finding functions.
	valgrind => 0
    },
    {
	name     => "tests.search.exe",                                                   # Tests of searching functions.
	valgrind => 0
    },
    {
	name     => "tests.sort.exe",                                                     # Tests of sorting functions.
	valgrind => 0
    },
    {
	name     => "tests.string_utilities.exe",                                         # Tests of string handling utilities.
	valgrind => 0
    },
    {
	name     => "tests.vectors.exe",                                                  # Tests of vector functions.
	valgrind => 0
    },
    {
	name     => "tests.cosmic_age.cosmological_constant.exe",                         # Tests of cosmic age calculations.
	valgrind => 0
    },
    {
	name     => "tests.cosmic_age.EdS.exe",                                           # .
	valgrind => 0
    },
    {
	name     => "tests.cosmic_age.open.exe",                                          # .
	valgrind => 0
    },
    {
	name     => "tests.cosmic_age.dark_energy.cosmological_constant.exe",             # .
	valgrind => 0
    },
    {
	name     => "tests.cosmic_age.dark_energy.omegaMinusOneThird.exe",                # .
	valgrind => 0
    },
    {
	name     => "tests.cosmic_age.dark_energy.closed.exe",                            # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.open.exe",                                  # Tests of spherical collapse calculations.
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.flat.exe",                                  # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.EdS.exe",                       # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.open.exe",                      # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.lambda.exe",                    # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminusTwoThirds.exe", # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminus0.6.exe",       # .
	valgrind => 0
    },
    {
	name     => "tests.spherical_collapse.dark_energy.constantEoSminus0.8.exe",       # .
	valgrind => 0
    },
    {
	name     => "tests.linear_growth.cosmological_constant.exe",                      # Tests of linear growth factor.
	valgrind => 0
    },
    {
	name     => "tests.linear_growth.EdS.exe",                                        # .
	valgrind => 0
    },
    {
	name     => "tests.linear_growth.open.exe",                                       # .
 	valgrind => 0
    },
    {
	name     => "tests.linear_growth.dark_energy.exe",                                # .
	valgrind => 0
    },
    {
	name     => "tests.halo_mass_function.Tinker.exe",                                # Tests of dark matter halo mass functions.
 	valgrind => 0
    },
    {
	name     =>"tests.comoving_distance.dark_energy.exe",                             # Tests of comoving distance calculations.
 	valgrind => 0
    },
    {
	name     =>"tests.comoving_distance.EdS.exe",                                     # .
	valgrind => 0
    },
    {
	name     => "tests.comoving_distance.open.exe",                                   # .
	valgrind => 0
    },
    {
	name     => "tests.Zhao2009_algorithms.dark_energy.exe",                          # Tests of Zhao et al. (2009) algorithms.
	valgrind => 0
    },
    {
	name     => "tests.Zhao2009_algorithms.EdS.exe",                                  # .
	valgrind => 0
    },
    {
	name     => "tests.Zhao2009_algorithms.open.exe",                                 # .
	valgrind => 0
    },
    {
	name     => "tests.NFW96_concentration.dark_energy.exe",                          # Tests of Navarro, Frenk & White (1996) halo concentration algorithm.
	valgrind => 0
    },
    {
	name     => "tests.Prada2011_concentration.exe",                                  # Tests of Prada et al. (2011) halo concentration algorithm.
	valgrind => 0
    },
    {
	name     => "tests.kepler_orbits.exe",                                            # Keplerian orbital parameter conversions.
	valgrind => 0
    },
    {
	name     => "tests.abundances.exe",                                               # Abundances objects.
	valgrind => 0
    },
    {
	name     => "tests.sigma.exe",                                                    # Sigma(M).
	valgrind => 0
    },
    {
	name     => "tests.power_spectrum.exe",                                           # Power spectrum.
	valgrind => 0
    },
    {
	name     => "tests.black_hole_fundamentals.exe",                                  # Black hole fundamentals.
	valgrind => 0
    },
    {
	name     => "tests.bug745815.exe",                                                # Regresssions.
	valgrind => 0
    },
    {
	name     => "tests.tree_branch_destroy.exe",                                      # Tests of merger tree walking.
	valgrind => 1,
	valgrindOptions => "--undef-value-errors=no"
    }
    );

# Run all executables.
foreach my $executable ( @executablesToRun ) {
    print lHndl "\n\n";
    print lHndl ":-> Running test: ".$executable->{'name'}."\n";
    &SystemRedirect::tofile("make ".$executable->{'name'},"allTestsBuild.tmp");
    my $buildSuccess = $?;
    print lHndl slurp("allTestsBuild.tmp");
    unlink("allTestsBuild.tmp");
    if ( $buildSuccess == 0 ) {
	# Run the test and copy any output to our log file.
	if ( $executable->{'valgrind'} == 1 ) {	    
	    &SystemRedirect::tofile("valgrind --error-exitcode=1 ".$executable->{'valgrindOptions'}." ".$executable->{'name'},"allTests.tmp");
	} else {
	    &SystemRedirect::tofile($executable->{'name'},"allTests.tmp");
	}
	my $runSuccess = $?;
	print lHndl "FAILED: running ".$executable->{'name'}." failed\n" if ( $runSuccess != 0 );
	print lHndl slurp("allTests.tmp");
	unlink("allTests.tmp",$executable->{'name'});
    } else {
	# Build failed, report an error in the log file.
	print lHndl "FAILED: building ".$executable->{'name'}." failed\n";
    }
}

# Build Galacticus itself.
print lHndl "\n\n";
print lHndl ":-> Building Galacticus...\n";
&SystemRedirect::tofile("make Galacticus.exe","allTestsBuild.tmp");
my $buildSuccess = $?;
print lHndl slurp("allTestsBuild.tmp");
unlink("allTestsBuild.tmp");
if ( $buildSuccess == 0 ) {
    # Run all tests.
    my @testDirs = ( "testSuite" );
    find(\&runTestScript,@testDirs);
} else {
    # Build failed, report an error in the log file.
    print lHndl "FAILED: building Galacticus.exe failed\n";
}

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
	
	switch ( $sendMethod ) {
	    case ( "sendmail" ) {
	    $msg->send;
	    }
	    case ( "smtp" ) {
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
}

exit $exitStatus;

sub runTestScript {
    # Run a test script.
    my $fileName = $_;
    chomp($fileName);

    # Test if this is a script to run.
    if ( $fileName =~ m/^test\-.*\.pl$/ && $fileName ne "test-all.pl" ) {
	print lHndl "\n\n:-> Running test script: ".$fileName."\n";
	&SystemRedirect::tofile($fileName,"allTests.tmp");
	print lHndl slurp("allTests.tmp");
	unlink("allTests.tmp");
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

