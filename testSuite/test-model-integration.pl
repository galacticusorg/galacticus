#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "../";
 $ENV{"GALACTICUS_ROOT_V094"} = getcwd()."/../";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
require Galacticus::HDF5;
require Stats::Histograms;

# Run a default Galacticus model as an integration test.
# Andrew Benson (05-December-2014)

# Run the model.
system("cd ..; scripts/aux/launch.pl testSuite/test-model-integration.xml");

# Create mass function bins.
my $xBins = pdl sequence(10)/2.0+8.0;

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file' } = "testSuite/outputs/test-model-integration/galacticus_0:1/galacticus.hdf5";
$dataSet->{'store'} = 0;
$dataSet->{'tree' } = "all";
&HDF5::Get_Parameters($dataSet    );
&HDF5::Get_Times     ($dataSet    );
&HDF5::Select_Output ($dataSet,0.0);
&HDF5::Count_Trees   ($dataSet    );
&HDF5::Get_Dataset($dataSet,['mergerTreeWeight','diskMassStellar','spheroidMassStellar']);
my $dataSets               = $dataSet->{'dataSets'};
my $logarithmicMassStellar = log10($dataSets->{'diskMassStellar'}+$dataSets->{'spheroidMassStellar'});
my $weight                 = $dataSets->{'mergerTreeWeight'};
(my $massFunction, my $error) = &Histograms::Histogram($xBins,$logarithmicMassStellar,$weight,differential => 1);
$massFunction /= log(10.0);
$error        /= log(10.0);

# Get Mercurial revision.
my $hgRevision = "Unknown";
open(my $hgHndl,"hg tip|");
while ( my $line = <$hgHndl> ) {
    if ( $line =~ m/^changeset:\s+(\d+)/ ) {
	$hgRevision = $1;
    }
}
close($hgHndl);

# Output the mass function.
system("mkdir -p testSuite/archive/integration");
open(my $outputFile,">testSuite/archive/integration/stellarMassFunctionZ0.000_r".$hgRevision.".txt");
for(my $i=0;$i<nelem($xBins);++$i) {
    print $outputFile $i."\t".$xBins->(($i))."\t".$massFunction->(($i))."\t".$error->(($i))."\n";
}
close($outputFile);

exit;
