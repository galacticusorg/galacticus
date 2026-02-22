#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Slurp qw(slurp);
use XML::Simple;
use List::ExtraUtils;
use Data::Dumper;
use Galacticus::Options;
use File::Find;

# Check for memory leaks in various models.
# Andrew Benson (14-April-2023)

# Extract options.
my %options = 
    (
     mpi => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Models to check.
my @models =
    (
     # {
     # 	 label      => "quickTest",
     # 	 parameters => "parameters/quickTest.xml"
     # },
     {
	 label      => "formationHalos",
	 parameters => "testSuite/parameters/memoryLeakFormationHalos.xml"
     },
     {
	 label      => "MCMC",
	 parameters => "testSuite/parameters/mcmcConfig.xml",
	 mpi        =>
	 {
	     processes => 4,
	     threads   => 1
	 }
     }
    );

# Iterate over models.
my $status = "SUCCESS";
foreach my $model ( @models ) {

    # Skip models that do not match MPI option.
    next
	if (
	    (   $options{'mpi'} && ! exists($model->{'mpi'}) )
	    ||
	    ( ! $options{'mpi'} &&   exists($model->{'mpi'}) )
	);
    
    # Run the model.
    print "Running model '".$model->{'label'}."'...\n";
    system("mkdir -p outputs/memoryLeaks/".$model->{'label'}."; cd ..; ".($model->{'mpi'} ? "export OMP_NUM_THREADS=".$model->{'mpi'}->{'threads'}."; mpirun --oversubscribe --allow-run-as-root --n ".$model->{'mpi'}->{'processes'}." " : "")."valgrind --leak-check=full --xml=yes --xml-file=testSuite/outputs/memoryLeaks/".$model->{'label'}."/memory-leaks-%p.xml ./Galacticus.exe ".$model->{'parameters'});
    unless ( $? == 0 ) {
    	print "\tFAILED:  model run:\n";
    } else {
    	print "\tSUCCESS: model run\n";
    }

    # Find the Valgrind outputs.
    my @valgrindOutputs;
    &File::Find::find(sub {push(@valgrindOutputs,$_) if ($_ =~ m/^memory\-leaks\-\d+\.xml$/)},"outputs/memoryLeaks/".$model->{'label'});
    
    # Read the Valgrind outputs.
    my $xml = new XML::Simple;
    foreach my $valgrindOutput ( @valgrindOutputs ) {
	print "\tAnalyzing output file '".$valgrindOutput."'\n";
	my $valgrind = eval { $xml->XMLin("outputs/memoryLeaks/".$model->{'label'}."/".$valgrindOutput) };
	if ( $@ ) {
           print "   malformed XML\n";
           next;
        }

	# Iterate over errors, looking for memory leaks.
	foreach my $error ( &List::ExtraUtils::as_array($valgrind->{'error'}) ) {
	    # Skip errors that are not memory leaks.
	    next
		unless (
		    $error->{'kind'} eq "Leak_PossiblyLost"
		    ||
		    $error->{'kind'} eq "Leak_DefinitelyLost"
		);
	    # Walk the stack frame.
	    my $ignore = 0;
	    my $frameFinal;
	    foreach my $frame ( &List::ExtraUtils::as_array($error->{'stack'}->{'frame'})) {
		$frameFinal = $frame
		    if ( ! defined($frameFinal) && $frame->{'obj'} =~ m/Galacticus\.exe/ );
		# Ignore reported possible leaks in libgomp.
		$ignore = 1
		    if ( $error->{'kind'} eq "Leak_PossiblyLost" && $frame->{'obj'} =~ m/libgomp/ );
		# Ignore reported leaks in libmpi.
		$ignore = 1
		    if ( $frame->{'obj'} =~ m/libmpi/ );
	    }
	    next
		if ( $ignore || ! defined($frameFinal) );
	    print "\t\tMemory leak (".$error->{'kind'}.") for model '".$model->{'label'}."' in: '".(exists($frameFinal->{'file'}) ? $frameFinal->{'file'} : "(UNKNOWN)")."' line ".(exists($frameFinal->{'line'}) ? $frameFinal->{'line'} : "(UNKNOWN)")."\n";
	    $status = "FAILED";
	}

    }

}

print $status.": memory leaks\n";

exit 0;
