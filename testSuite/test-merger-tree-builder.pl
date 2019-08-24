#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use File::Slurp qw(slurp);
use List::Util;

# Compare conditional mass functions with a reference set.
# Andrew Benson (06-August-2019)

# Specify models.
my @models =
    (
     {
	 type                => 'reference'                                                              ,
	 fileName            => 'data/mergerTreeBuilding/test-merger-tree-builder-reference.hdf5'        ,
	 color               => 'redYellow'
     },
     {
	 type                => 'referenceGeneric'                                                       ,
	 fileName            => 'data/mergerTreeBuilding/test-merger-tree-builder-reference-generic.hdf5',
	 color               => 'blackGray'
     },
     {
     	 type                => 'Cole et al. (2000); intervalStep=true'                                  ,
     	 fileName            => 'outputs/mergerTreeBuilderCole2000_intervalStepTrue.hdf5'                ,
     	 parameterFile       => 'parameters/mergerTreeBuilderCole2000_intervalStepTrue.xml'              ,
     	 reference           => 'reference'                                                              ,
     	 color               => 'cornflowerBlue'                                                         ,
     	 toleranceSigma      => 6.0                                                                      ,
     	 toleranceFractional => 1.0e-3
     },
     {
     	 type                => 'Cole et al. (2000); intervalStep=false'                                 ,
     	 fileName            => 'outputs/mergerTreeBuilderCole2000_intervalStepFalse.hdf5'               ,
     	 parameterFile       => 'parameters/mergerTreeBuilderCole2000_intervalStepFalse.xml'             ,
     	 color               => 'mediumSeaGreen'                                                         ,
     	 reference           => 'reference'                                                              ,
     	 toleranceSigma      => 6.0                                                                      ,
     	 toleranceFractional => 1.0e-3
     },
     {
     	 type                => 'genericLinearBarrier'                                                   ,
     	 fileName            => 'outputs/mergerTreeBuilderGenericLinearBarrier.hdf5'                     ,
     	 parameterFile       => 'parameters/mergerTreeBuilderGenericLinearBarrier.xml'                   ,
     	 color               => 'peachPuff'                                                              ,
     	 reference           => 'referenceGeneric'                                                       ,
     	 toleranceSigma      => 6.0                                                                      ,
     	 toleranceFractional => 1.0e-3
     },
     {
     	 type                => 'genericSolver'                                                          ,
     	 fileName            => 'outputs/mergerTreeBuilderGenericSolver.hdf5'                            ,
     	 parameterFile       => 'parameters/mergerTreeBuilderGenericSolver.xml'                          ,
     	 color               => 'orange'                                                                 ,
     	 reference           => 'referenceGeneric'                                                       ,
     	 toleranceSigma      => 9.0                                                                      ,
     	 toleranceFractional => 1.0e-3
     }
    );

# Run models.
my $document = << 'OPENER';
<parameterGrid>
  <emailReport>no</emailReport>
  <doAnalysis>no</doAnalysis>
  <modelRootDirectory>testSuite/outputs/test-merger-tree-builder</modelRootDirectory>
  <baseParameters>testSuite/parameters/mergerTreeBuilderCole2000_intervalStepTrue.xml</baseParameters>
  <launchMethod>pbs</launchMethod>
  <pbs>
    <mpiLaunch>no</mpiLaunch>
    <ompThreads>16</ompThreads>
    <maxJobsInQueue>100</maxJobsInQueue>
    <postSubmitSleepDuration>1</postSubmitSleepDuration>
    <jobWaitSleepDuration>10</jobWaitSleepDuration>
  </pbs>
OPENER
my $i = -1;
foreach my $model ( @models ) {
    if ( exists($model->{'parameterFile'}) ) {
	$document .= slurp($model->{'parameterFile'});
	++$i;
	$model->{'tmpName'} = "galacticus_".$i.":1";
    }
}
$document .= "</parameterGrid>\n";
open(my $launchFile,">outputs/test-merger-tree-builder.xml");
print $launchFile $document;
close($launchFile);
system("cd ..; mkdir -p testSuite/outputs/test-merger-tree-builder; scripts/aux/launch.pl testSuite/outputs/test-merger-tree-builder.xml");

# Check for failed models.
system("grep -q -i fatal outputs/test-merger-tree-builder/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-merger-tree-builder/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS: model run\n";
}

# Rename models.
foreach my $model ( @models ) {
    if ( exists($model->{'parameterFile'}) ) {
	system("mv outputs/test-merger-tree-builder/".$model->{'tmpName'}."/galacticus.hdf5 ".$model->{'fileName'});
    }
}

# Read test and reference data.
my $modelData;
foreach my $model ( @models ) {
    $modelData->{$model->{'type'}}->{'file' } = new PDL::IO::HDF5($model->{'fileName'});
    $modelData->{$model->{'type'}}->{'group'} = $modelData->{$model->{'type'}}->{'file'}->group("conditionalMassFunction");
    foreach my $datasetName ( 'massParent', 'massRatio', 'redshiftParent', 'redshiftProgenitor', 'conditionalMassFunction', 'conditionalMassFunctionError' ) {
	$modelData->{$model->{'type'}}->{$datasetName} = $modelData->{$model->{'type'}}->{'group'}->dataset($datasetName)->get();
    }
    $modelData->{$model->{'type'}}->{'matchMeasureMaximum'       } = pdl 0.0;
    $modelData->{$model->{'type'}}->{'deviationFractionalMaximum'} = pdl 0.0;
}

# Iterate over redshifts.
system("mkdir -p outputs/mergerTreeBuilder");
for(my $iRedshift=0;$iRedshift<nelem($modelData->{'reference'}->{'redshiftParent'});++$iRedshift) {
    # Iterate over parent masses.
    for(my $iMass=0;$iMass<nelem($modelData->{'reference'}->{'massParent'});++$iMass) {
	# Begin construction of plot.
	my $plot;
	my $plotFilePDF = "outputs/mergerTreeBuilder/conditionalMassFunctions_".$iRedshift."_".$iMass.".pdf";
	(my $plotFileTeX = $plotFilePDF) =~ s/\.pdf$/.tex/;
	open(my $gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
	print $gnuPlot "set output '".$plotFileTeX."'\n";
	print $gnuPlot "set title '\$\\log_{10}(M_\\mathrm{par}/\\mathrm{M}_\\odot)=".sprintf("%5.2f",$modelData->{'reference'}->{'massParent'}->log10()->(($iMass)))."; z_\\mathrm{pro}=".sprintf("%4.2f",$modelData->{'reference'}->{'redshiftProgenitor'}->(($iRedshift)))."\$' offset 0,-1\n";
	print $gnuPlot "set xlabel 'Progenitor mass ratio; \$x = M_\\mathrm{pro}/M_\\mathrm{par}\$'\n";
	print $gnuPlot "set ylabel 'Conditional mass function; \$\\mathrm{d}f/\\mathrm{d}\\log_{10} x\$' offset -0.5,0\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.12,0.67\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [1.0e-3:2.0]\n";
	print $gnuPlot "set yrange [1.0e-2:1.0e1]\n";
	print $gnuPlot "set pointsize 1.0\n";
	foreach my $model ( @models ) {
	    # Construct measure of how well the model matches the reference model.
	    if ( exists($model->{'reference'}) ) {
		my $select = which
		    (
		     ($modelData->{$model->{'reference'}}->{'conditionalMassFunction'}->(($iRedshift),($iMass),:) > 0.0)
		     &
		     ($modelData->{$model->{'type'     }}->{'conditionalMassFunction'}->(($iRedshift),($iMass),:) > 0.0)
		    );
		my $matchMeasure        = +   ($modelData->{$model->{'type'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),$select)   -$modelData->{$model->{'reference'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),$select)   )**2
		    /   ($modelData->{$model->{'type'}}->{'conditionalMassFunctionError'}->(($iRedshift),($iMass),$select)**2+$modelData->{$model->{'reference'}}->{'conditionalMassFunctionError'}->(($iRedshift),($iMass),$select)**2)   ;
		my $deviationFractional = +abs($modelData->{$model->{'type'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),$select)   -$modelData->{$model->{'reference'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),$select)   )
		    /max($modelData->{$model->{'type'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),$select)   ,$modelData->{$model->{'reference'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),$select)   )   ;
		$modelData->{$model->{'type'}}->{'matchMeasureMaximum'       } .= $matchMeasure       ->maximum()
		    if ( nelem($matchMeasure) > 0 && $matchMeasure       ->maximum() > $modelData->{$model->{'type'}}->{'matchMeasureMaximum'       } );
		$modelData->{$model->{'type'}}->{'deviationFractionalMaximum'} .= $deviationFractional->maximum()
		    if ( nelem($matchMeasure) > 0 && $deviationFractional->maximum() > $modelData->{$model->{'type'}}->{'deviationFractionalMaximum'} );
	    }
	    &GnuPlot::PrettyPlots::Prepare_Dataset
		(
		 \$plot                                                                                                  ,
		 $modelData             ->{$model->{'type'}}->{'massRatio'                   }                           ,
		 $modelData             ->{$model->{'type'}}->{'conditionalMassFunction'     }->(($iRedshift),($iMass),:),
		 errorDown => $modelData->{$model->{'type'}}->{'conditionalMassFunctionError'}->(($iRedshift),($iMass),:),
		 errorUp   => $modelData->{$model->{'type'}}->{'conditionalMassFunctionError'}->(($iRedshift),($iMass),:),
		 style     => "point"                                                                                    ,
		 symbol    => [6,7]                                                                                      ,
		 weight    => [3,1]                                                                                      ,
		 color     => $GnuPlot::PrettyPlots::colorPairs{$model->{'color'}},
		 title     =>                                   $model->{'type'}
		);
	}
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
    }
}
# Report on the maximum deviations.
my $typeLength = &List::Util::max(map {length($_->{'type'})} @models);
foreach my $model ( @models ) {
    next
	unless ( exists($model->{'reference'}) );
    my $status = ($modelData->{$model->{'type'}}->{'matchMeasureMaximum'} > $model->{'toleranceSigma'} && $modelData->{$model->{'type'}}->{'deviationFractionalMaximum'} > $model->{'toleranceFractional'}) ? "FAILED" : "success";
    print $status.": ".$model->{'type'}.(" " x ($typeLength-length($model->{'type'}))).": Δσₘₐₓ = ".sprintf("%3.1f",$modelData->{$model->{'type'}}->{'matchMeasureMaximum'})." {limit: ".sprintf("%3.1f",$model->{'toleranceSigma'})."}; Δ(log f)ₘₐₓ = ".sprintf("%7.1e",$modelData->{$model->{'type'}}->{'deviationFractionalMaximum'})." {limit: ".sprintf("%7.1e",$model->{'toleranceFractional'})."}\n";
}

exit 0;
