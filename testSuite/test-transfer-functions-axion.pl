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

my @models =
    (
     {
	 label                   => "CDM"                      ,
	 plotLabel               => "CDM"                      ,
	 testPeaks               => "no"
     },
     {
	 label                   => "AxionCAMB"                ,
	 plotLabel               => "AxionCAMB"                ,
	 color                   => "redYellow"                ,
	 testPeaks               => "no"
     },
     {
	 label                   => "AxionMurgia2017"          ,
	 plotLabel               => "Murgia et al. (2017)"     ,
	 color                   => "peachPuff"                ,
	 testPeaks               => "no"                       ,
	 toleranceCutOff         => 0.10
     },
     {
	 label                   => "AxionHu2000"              ,
	 plotLabel               => "Hu et al. (2000)"         ,
	 color                   => "mediumSeaGreen"           ,
	 testPeaks               => "yes"                      ,
	 toleranceCutOff         => 0.160                      ,
	 tolerancePeakWavenumber => 0.090                      ,
	 tolerancePeakAmplitude  => 0.800
     },
     {
	 label                   => "AxionPassaglia2022"       ,
	 plotLabel               => "Passaglia \\\\& Hu (2022)",
	 color                   => "cornflowerBlue"           ,
	 testPeaks               => "yes"                      ,
	 toleranceCutOff         => 0.025                      ,
	 tolerancePeakWavenumber => 0.015                      ,
	 tolerancePeakAmplitude  => 0.150
     }
    );
my $transferFunctions;
foreach my $model ( @models ) {
    # Run models.
    system("mkdir -p outputs; cd ..; ./Galacticus.exe testSuite/parameters/powerSpectrum".$model->{'label'}.".xml");
    unless ( $? == 0 ) {
    	print "FAIL: failed to run model 'powerSpectrum".$model->{'label'}.".xml'\n";
    	exit;
    }
    # Read data.
    my $galacticus = new PDL::IO::HDF5("outputs/powerSpectrum".$model->{'label'}.".hdf5");
    my $output     = $galacticus->group('Outputs/Output1');
    $transferFunctions->{$model->{'label'}}->{$_} = $output->dataset($_)->get()
	foreach ( "wavenumber", "transferFunction" );
    # Normalize the transfer function.
    $transferFunctions->{$model->{'label'}}->{'transferFunction'} /=     $transferFunctions->{$model->{'label'}}->{'transferFunction'}->((0));
    $transferFunctions->{$model->{'label'}}->{'transferFunction'} .= abs($transferFunctions->{$model->{'label'}}->{'transferFunction'}      );
    $transferFunctions->{$model->{'label'}}->{'transferFunction'} /=     $transferFunctions->{         'CDM'   }->{'transferFunction'}
        unless ( $model->{'label'} eq "CDM" );
}

# Create a plot of the window function.
my $plot;
my $gnuPlot;
my $plotFileTeX = "outputs/transferFunctionsAxion.tex";
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 3.5in,3.5in\n";
print $gnuPlot "set output '".$plotFileTeX."'\n";
print $gnuPlot "set lmargin screen 0.20\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.9,0.93\n";
print $gnuPlot "set xlabel '\$ k \$ Mpc\$^{-1}\$'\n";
print $gnuPlot "set ylabel '\$ (k/10\\mathrm{Mpc}^{-1})^6 |T(k)| \$\n";
print $gnuPlot "set logscale y\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [3.0e+0:15.0]\n";
print $gnuPlot "set yrange [1.0e-3:3.0e-1]\n";
print $gnuPlot "set pointsize 1.0\n";
foreach my $model ( @models ) {
    next
	if ( $model->{'label'} eq "CDM" );
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	    \$plot                                                           ,
	    +  $transferFunctions->{$model->{'label'}}->{'wavenumber'      } ,
	    +  $transferFunctions->{$model->{'label'}}->{'transferFunction'}
	    *(
	      +$transferFunctions->{$model->{'label'}}->{'wavenumber'      }
	      /10.0
	    )**6,
	    style  => "line"                                                 ,
	    weight => [2,1]                                                  ,
	    color  => $GnuPlot::PrettyPlots::colorPairs{$model->{'color'}}   ,
	    title  => $model->{'plotLabel'}
	);
}
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);

# Perform tests.
foreach my $model ( @models ) {
    next
	if ( $model->{'label'} eq "CDM" );
    # Find locations of peaks.
    @{$transferFunctions->{$model->{'label'}}->{'peaks'}} = ();
    for(my $i=1;$i<nelem($transferFunctions->{$model->{'label'}}->{'wavenumber'})-1;++$i) {
	push(
	    @{$transferFunctions->{$model->{'label'}}->{'peaks'}},
	    {
		wavenumber       => $transferFunctions->{$model->{'label'}}->{'wavenumber'      }->(($i)),
		transferFunction => $transferFunctions->{$model->{'label'}}->{'transferFunction'}->(($i))
	    }
	    )
	    if (
		$transferFunctions->{$model->{'label'}}->{'transferFunction'}->(($i)) > $transferFunctions->{$model->{'label'}}->{'transferFunction'}->(($i-1))
		&&
		$transferFunctions->{$model->{'label'}}->{'transferFunction'}->(($i)) > $transferFunctions->{$model->{'label'}}->{'transferFunction'}->(($i+1))
	    );
    }
    # Test amplitudes at ~4 Mpc⁻¹ against AxionCAMB.
    if ( exists($model->{'toleranceCutOff'}) ) {
	my $wavenumberTarget = pdl 4.0;
	my $indexTarget      = minimum_ind(abs($transferFunctions->{$model->{'label'}}->{'wavenumber'}-$wavenumberTarget));
	my $error =
	    abs(
		+$transferFunctions->{$model->{'label'    }}->{'transferFunction'}->(($indexTarget))
		-$transferFunctions->{         'AxionCAMB' }->{'transferFunction'}->(($indexTarget))
	    )
	    /    $transferFunctions->{         'AxionCAMB' }->{'transferFunction'}->(($indexTarget));
	my $errorPercent = sprintf("%4.1f%%",100.0*$error);
	my $status       = $error > $model->{'toleranceCutOff'} ? "FAIL" : "SUCCESS";
	print $model->{'label'}." T(k=4 Mpc⁻¹): ".$errorPercent." (".$status.")\n";
    }
    # Test locations and amplitudes of peaks.
    if ( exists($model->{'tolerancePeakWavenumber'}) ) {
	my @subscripts = ("₀","₁","₂","₃","₄","₅","₆","₇","₈","₉");
	foreach(my $peak=0;$peak<2;++$peak) {
	    my $errorWavenumber =
		abs(
		    +$transferFunctions->{$model->{'label'    }}->{'peaks'}->[$peak]->{'wavenumber'      }
		    -$transferFunctions->{         'AxionCAMB' }->{'peaks'}->[$peak]->{'wavenumber'      }
		)
		/    $transferFunctions->{         'AxionCAMB' }->{'peaks'}->[$peak]->{'wavenumber'      };
	    my $errorAmplitude =
		abs(
		    +$transferFunctions->{$model->{'label'    }}->{'peaks'}->[$peak]->{'transferFunction'}
		    -$transferFunctions->{         'AxionCAMB' }->{'peaks'}->[$peak]->{'transferFunction'}
		)
		/    $transferFunctions->{         'AxionCAMB' }->{'peaks'}->[$peak]->{'transferFunction'};
	    my $errorWavenumberPercent = sprintf("%4.1f%%",100.0*$errorWavenumber);
	    my $errorAmplitudePercent  = sprintf("%4.1f%%",100.0*$errorAmplitude );
	    my $statusWavenumber       = $errorWavenumber > $model->{'tolerancePeakWavenumber'} ? "FAIL" : "SUCCESS";
	    my $statusAmplitude        = $errorAmplitude  > $model->{'tolerancePeakAmplitude' } ? "FAIL" : "SUCCESS";
	    print $model->{'label'}." k".$subscripts[$peak+1].": ".$errorWavenumberPercent." (".$statusWavenumber.")\n";
	    print $model->{'label'}." T".$subscripts[$peak+1].": ".$errorAmplitudePercent ." (".$statusAmplitude .")\n";
	}
    }
}

exit;

