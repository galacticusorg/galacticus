#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use PDL::Ufunc;

# Test tidally truncated NFW profile fits.
# Andrew Benson (15-August-2023)

# Run the model.
system("mkdir -p outputs");
system("cd ..; ./Galacticus.exe testSuite/parameters/tidallyTruncatedNFWFit.xml");
die("FAILED: model failed to run")
   unless ( $? == 0 );

# Read model data.
my $model = new PDL::IO::HDF5("outputs/tidallyTruncatedNFWFit.hdf5");
my $nodes = $model->group('Outputs/Output1/nodeData');
my $data;
$data->{$_} = $nodes->dataset($_)->get()
    foreach ( "basicMass", "darkMatterProfileScale", "darkMatterOnlyRadiusVirial", "radiusTidalTruncationNFW", "densityProfile", "densityProfileRadius" );
$data->{'concentration'       } = +$data->{'darkMatterOnlyRadiusVirial'}
                                  /$data->{'darkMatterProfileScale'    };
$data->{'densityNormalization'} = +$data->{'basicMass'                 }
                                  /$data->{'darkMatterProfileScale'    }**3
                                  /4.0/PI
                                  /(log(1.0+$data->{'concentration'})-$data->{'concentration'}/(1.0+$data->{'concentration'}));
$data->{'metricUntruncated'   } = pdl zeros($data->{'basicMass'});
$data->{'metricTruncated'     } = pdl zeros($data->{'basicMass'});

# Iterate over all halos/subhalos.
for(my $i=0;$i<nelem($data->{'basicMass'});++$i) {
    my $radii              =        $data->{'densityProfileRadius'    }->(:,($i));
    my $densityTarget      =        $data->{'densityProfile'          }->(:,($i));
    my $xs                 = $radii/$data->{'darkMatterProfileScale'  }->(  ($i));
    my $xt                 = $radii/$data->{'radiusTidalTruncationNFW'}->(  ($i));
    my $densityUntruncated = $data->{'densityNormalization'}->(($i))/$xs/(1.0+$xs)**2             ;
    my $densityTruncated   = $data->{'densityNormalization'}->(($i))/$xs/(1.0+$xs)**2/(1.0+$xt**2);
    $data->{'metricUntruncated'}->(($i)) .= sum(log10($densityTarget/$densityUntruncated)**2)/nelem($radii);
    $data->{'metricTruncated'  }->(($i)) .= sum(log10($densityTarget/$densityTruncated  )**2)/nelem($radii);
}
my $status = $data->{'metricTruncated'}->median() < 0.0125 ? "succeeded" : "FAILED";
print $status.": tidally truncated NFW fit\n";

exit;
