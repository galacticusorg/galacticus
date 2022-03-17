#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;

# Generates a set of input files for Galacticus which divide a single model between multiple workers.
# Andrew Benson (16-Oct-2011)

# Get the input parameter file and the number of workers to split between.
die("Usage: Split_Models_For_Workers.pl <parameterFile> <workerCount>") unless ( scalar(@ARGV) == 2 );
my $parameterFile = $ARGV[0];
my $workerCount   = $ARGV[1];

# Read in the input parameter file.
my $xml                 = new XML::Simple;
my $parameters         = $xml->XMLin($parameterFile);  
my $outputFileNameBase = "galacticus.hdf5";
my $stateFileNameBase  = "none";
$outputFileNameBase    = $parameters->{'parameter'}->{'outputFileName'}->{'value'}
             if ( exists($parameters->{'parameter'}->{'outputFileName'}) );
$stateFileNameBase     = $parameters->{'parameter'}->{'stateFileRoot'           }->{'value'}
             if ( exists($parameters->{'parameter'}->{'stateFileRoot'           }) );
for(my $worker=1;$worker<=$workerCount;++$worker) {
    ($parameters->{'parameter'}->{'outputFileName'}->{'value'}  = $outputFileNameBase ) =~ s/\.hdf5/_$worker\.hdf5/;
    $parameters ->{'parameter'}->{'stateFileRoot'           }->{'value'}  = $stateFileNameBase;
    $parameters ->{'parameter'}->{'stateFileRoot'           }->{'value'} .= "_".$worker unless ( $stateFileNameBase eq "none" );
    $parameters ->{'parameter'}->{'treeEvolveWorkerCount'   }->{'value'}  = $workerCount;
    $parameters ->{'parameter'}->{'treeEvolveWorkerNumber'  }->{'value'}  = $worker;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
    (my $outputFile = $parameterFile) =~ s/\.xml/_$worker\.xml/;
    open(outHndl,">".$outputFile);
    print outHndl $xmlOutput->XMLout($parameters);
    close(outHndl);
}

exit;
