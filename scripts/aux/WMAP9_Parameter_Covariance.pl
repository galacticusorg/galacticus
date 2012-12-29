#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use Data::Dumper;
use XML::Simple;
use DateTime;

# Compute the cosmological parameter covariance matrix from the WMAP-9 year results.
# Andrew Benson (28-December-2012)

# Create a working directory.
my $workDirectory = $galacticusPath."/aux/WMAP-9";
system("mkdir -p ".$workDirectory);

# Download the Monte Carlo Markov Chains.
my $chainsURL = "http://lambda.gsfc.nasa.gov/data/map/dr5/dcp/chains/wmap_lcdm_wmap9_spt_act_snls3_chains_v5.tar.gz";
system("wget ".$chainsURL." -O ".$workDirectory."/chains.tar.gz")
    unless ( -e $workDirectory."/chains.tar.gz" );

# Unpack the chains.
system("cd ".$workDirectory."; tar xvfz chains.tar.gz")
    unless ( -e $workDirectory."/description.txt" );

# Read the weights.
(my $weights) = rcols($workDirectory."/weight_including_2012bao_h0",1);
$weights /= $weights->dsum();

# Read the datasets.
my @datasets =
    (
     {
	 label       => "omega_B",
	 file        => "omegabh2",
	 units       => "none",
	 description => "Baryon density parameter, Omega_b*h^2."
     },
     {
	 label       => "omega_M",
	 file        => "omegamh2",
	 units       => "none",
	 description => "Matter density parameter, Omega_M*h^2."
     },
     {
      	 label       => "tau",
      	 file        => "tau",
	 units       => "none",
	 description => "Optical depth to reionization."
     },
     {
      	 label       => "H_0",
      	 file        => "H0",
	 units       => "km/s/Mpc",
	 description => "Hubble parameter."
     },
     {
      	 label       => "n_s",
      	 file        => "ns002",
	 units       => "none",
	 description => "Primordial power spectrum spectral index."
     },
     {
      	 label       => "sigma_8",
      	 file        => "sigma8",
	 units       => "none",
	 description => "Root-variance of mass fluctuations in 8Mpc/h radius top hat spheres."
     }
    );
foreach my $dataset ( @datasets ) {
    my $file = $workDirectory."/".$dataset->{'file'};
    ($dataset->{'data'}) = rcols($file,1);
    $dataset->{'mean'} = sum($dataset->{'data'}*$weights);
}

# Compute covariances.
for(my $i=0;$i<scalar(@datasets);++$i) {
    my $datasetI = $datasets[$i];
    for(my $j=0;$j<=$i;++$j) {
	my $datasetJ = $datasets[$j];
	$datasetI->{'covariance'}->{$datasetJ->{'label'}} = 
	    dsum(
		($datasetI->{'data'}-$datasetI->{'mean'})
		*
		($datasetJ->{'data'}-$datasetJ->{'mean'})
		*
		$weights
	    );
    }
}

# Compute correlations.
for(my $i=0;$i<scalar(@datasets);++$i) {
    my $datasetI = $datasets[$i];
    for(my $j=0;$j<=$i;++$j) {
	my $datasetJ = $datasets[$j];
	$datasetI->{'correlation'}->{$datasetJ->{'label'}} = 
	    $datasetI->{'covariance'}->{$datasetJ->{'label'}}
	/
	    sqrt
	    (
	     $datasetI->{'covariance'}->{$datasetI->{'label'}}
	     *
	     $datasetJ->{'covariance'}->{$datasetJ->{'label'}}
	    );
	
    }
}

# Construct output data structure.
my $output;
for(my $i=0;$i<scalar(@datasets);++$i) {
    my $datasetI = $datasets[$i];
    my @parameters;
    for(my $j=0;$j<$i;++$j) {
	my $datasetJ = $datasets[$j];
	push(
	    @parameters,
	    {
		label       =>                              $datasetJ->{'label'} ,
		covariance  => $datasetI->{'covariance' }->{$datasetJ->{'label'}},
		correlation => $datasetI->{'correlation'}->{$datasetJ->{'label'}}
	    }
	    );
    }
    push(
	@{$output->{'parameter'}},
	{
	    label             =>      $datasetI->{'label'      }                         ,
	    description       =>      $datasetI->{'description'}                         ,
	    units             =>      $datasetI->{'units'      }                         ,
	    mean              =>      $datasetI->{'mean'       }                         ,
	    standardDeviation => sqrt($datasetI->{'covariance' }->{$datasetI->{'label'}}),
	    parameter         => \@parameters
	}
    );
}

# Add metadata.
@{$output->{'description'}} = 
    (
     "WMAP-9 Lambda CDM cosmological parameter best fit values with covariance and correlation matrices.",
     "Uses: wmap9+spt+act+snls3+bao+h0 constraints."
    );
@{$output->{'url'        }} =
    (
     "http://lambda.gsfc.nasa.gov/product/map/dr5/params/lcdm_wmap9_spt_act_snls3.cfm",
     "http://lambda.gsfc.nasa.gov/data/map/dr5/dcp/chains/wmap_lcdm_wmap9_spt_act_snls3_chains_v5.tar.gz"
    );
@{$output->{'createdBy'  }} =
    (
     "Galacticus",
     "scripts/aux/WMAP9_Parameter_Covariance.pl"
    );
$output->{'source'     } =
    "Computed from Monte Carlo Markov Chains.";
my $dt = DateTime->now->set_time_zone('local');
(my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
$output->{'timestamp'  } = $now;

# Serialize data to XML.
my $xml = new XML::Simple(NoAttr=>1, RootName=>"parameters");
open(oHndl,">".$galacticusPath."/cosmology/Cosmological_Parameters_WMAP-9.xml");
print oHndl $xml->XMLout($output);
close(oHndl);

exit;
