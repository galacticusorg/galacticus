# Contains a Perl module which implements calculations of probabilistic
# gravitational lensing amplification for galaxies.

package LensingAmplification;
use strict;
use warnings;
use PDL;
use PDL::Math;
use PDL::NiceSlice;
use PDL::GSL::INTEG;
use PDL::GSL::INTERP;
require PDL::Fit::LM;
use PDL::IO::HDF5;
use Galacticus::HDF5;
use Astro::Cosmology;
use DateTime;
use Data::Dumper;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^lensingAmplification\$" => \&LensingAmplification::Get_Amplification
    );

# Module variables.
my $hubble0;
my $omegaM0;
my $omegaDE0;
my $sourceRedshift;
my $cosmology;
my $kappaEmpty;

my $status = 1;
$status;

sub Get_Amplification {
    # Get the data structure.
    my $dataBlock   = shift;

    # Get the name of the requested dataset.
    my $dataSetName = $_[0];

    # Extract cosmological parameters.
    $hubble0  = $dataBlock->{'parameters'}->{'H_0'         };
    $omegaM0  = $dataBlock->{'parameters'}->{'Omega_Matter'};
    $omegaDE0 = $dataBlock->{'parameters'}->{'Omega_DE'    };
    $cosmology = Astro::Cosmology->new(
	omega_matter => $omegaM0,
	omega_lambda => $omegaDE0,
	H0           => $hubble0
	);

    # Construct the name of a file to which the lensing PDFs can be stored.
    my $lensingPdfDirectory = "data/gravitationalLensingPDFs";
    system("mkdir -p ".$lensingPdfDirectory);
    my $omegaM0fixedPrecision  = sprintf("%6.4f",$omegaM0 );
    my $omegaDE0fixedPrecision = sprintf("%6.4f",$omegaDE0);
    my $hubble0fixedPrecision  = sprintf("%6.3f",$hubble0 );
    my $lensingPdfFile      = $lensingPdfDirectory."/lensingPDF:OmegaMatter".$omegaM0fixedPrecision.":OmegaDE".$omegaDE0fixedPrecision.":H0".$hubble0fixedPrecision.".hdf5";

    # Initiate variables to hold the amplification PDFs and associated data.
    my $gridRedshift;
    my $gridProbability;
    my $amplificationPDFs;

    # Check for a pre-existing file of amplification PDFs.
    if ( -e $lensingPdfFile ) {
        # File exists, read the PDFs from file.
	my $hdfFile = new PDL::IO::HDF5($lensingPdfFile);
	$gridRedshift      = $hdfFile->dataset("redshift"        )->get();
	$gridProbability   = $hdfFile->dataset("probability"     )->get();
	$amplificationPDFs = $hdfFile->dataset("amplificationPDF")->get();
    } else {
	# File does not exist, so compute and store the PDFs.
	print "Galacticus::LensingAmplification: computing amplification PDFs (will be stored for future reference)\n";
	
        # Create a grid of redshifts on which we will build lensing amplification CDFs.
	my $redshiftLow   = pdl  0.1;
	my $redshiftHigh  = pdl 20.0;
	my $redshiftCount = 100;
	$gridRedshift     = pdl sequence($redshiftCount)*($redshiftHigh-$redshiftLow)/($redshiftCount-1)+$redshiftLow;
	
	# Get associated RMS convergence at each redshift.
	my $gridKappaRMS  = &Convergence_RMS($gridRedshift);
	
	# Define a grid of amplifications.
	my $amplificationLow   = 1.0e-1;
	my $amplificationHigh  = 1.0e+2;
	my $amplificationCount = 10000;
	my $gridAmplification  = pdl exp(sequence($amplificationCount)*log($amplificationHigh/$amplificationLow)/($amplificationCount-1)+log($amplificationLow));
	
	# Define a grid of probabilities. We construct this to have high resolution in the tails of the distribution to help catch rare events.
	my $probabilityLow   = 1.0e-6;
	my $probabilityHigh  = 1.0e+0;
	my $probabilityCount = 10000;
	$gridProbability = pdl [0.0];
	$gridProbability = $gridProbability->append(    exp(sequence($probabilityCount/2-1)*log(0.5/$probabilityLow)/($probabilityCount/2-1)+log($probabilityLow)));
	$gridProbability = $gridProbability->append(1.0-exp(sequence($probabilityCount/2-1)*log($probabilityLow/0.5)/($probabilityCount/2-1)+log(0.5            )));
	$gridProbability = $gridProbability->append(1.0);

	# Find the convergence corresponding to the amplification;
	my $gridConvergence  = 1.0-1.0/sqrt($gridAmplification);
	
	# Create an array to hold all PDFs.
	$amplificationPDFs = pdl zeroes($redshiftCount+1,$amplificationCount);

	# Create arrays to hold the fit parameters.
	my $omegaKappas = pdl zeroes(nelem($gridRedshift));
	my $aKappas     = pdl zeroes(nelem($gridRedshift));
	
	# Compute kappa_empty for each redshift.
	for(my $i=0;$i<nelem($gridRedshift);++$i) {
	    # Write a message.
	    my $n = $i+1;
	    my $redshiftFixedPrecision = sprintf("%5.2f",$gridRedshift->index($i));
	    print "-> Computing amplification PDF for z=".$redshiftFixedPrecision." (".$n." of ".nelem($gridRedshift).")\n";

	    # Copy the redshift to the source redshift.
	    $sourceRedshift = $gridRedshift->index($i);
	    
	    # Compute the convergence for an empty beam.
	    ($kappaEmpty,my $abserr,my $ierr) = gslinteg_qag(\&Kappa_Empty_Integrand,0.0,$sourceRedshift,1.0e-3,
							     0.0,1000,6);
	    
	    # Find the parameters of the convergence PDF using a minimization routine.
	    my $moments = pdl [1.0,2.0];
	    my $targets = pdl [-2.0*$gridKappaRMS->index($i)**2,$gridKappaRMS->index($i)**2];
	    
	    # Use a simple grid search to find the parameters which give the
	    # correct moments for the distribution. Ideally we'd use some
	    # minimization algorithm, but implementations in PDL seem to be
	    # lacking (or untrustworthy).
	    my $bestParams = pdl zeroes(2);
	    my $bestDiff   = pdl 1.0e30;
	    my $values     = pdl zeroes(2);
	    my $Dvalues    = pdl zeroes(2,2);
	    my $params     = pdl [0.0,0.0];
	    for(my $k=0;$k<500;++$k) {
		$params->index(0) .= $k*1.0/(500.0-1.0)-0.75;
		for(my $j=0;$j<500;++$j) {
		    $params->index(1) .= $j*1.0/(500.0-1.0)+0.25;
		    &Fit_Parameters($moments,$params,$values,$Dvalues);
		    my $diff = sum(($targets-$values)**2/$gridKappaRMS->index($i)**4);
		    if ( $diff < $bestDiff && $values->index(1) > 0.0 ) {
			$bestDiff   .= $diff;
			$bestParams .= $params;
		    }
		}
	    }
	    my $weight = 1.0;
	    my $pf;
	    if ( $sourceRedshift >= 1.0 ) {
		(my $yf,$pf,my $cf,my $if) = lmfit($moments, $targets, $weight, \&Fit_Parameters, $bestParams, {Maxiter => 3000, Eps => 1.0e-3});
	    } else {
		$pf = $bestParams;
	    }
	    my $omegaKappa .= 10.0**($pf->index(0));
	    my $aKappa     .= 10.0**($pf->index(1));
	    $omegaKappas->index($i) .= $omegaKappa;
	    $aKappas    ->index($i) .= $aKappa;
	    &Fit_Parameters($moments,$pf,$values,$Dvalues);
	    print "              omega_kappa = ".$omegaKappa."\n";
	    print "                  A_kappa = ".$aKappa."\n";
	    print "   Computed moments (1,2) = ".$values."\n";
	    print "     Target moments (1,2) = ".$targets."\n";
	    
	    # Compute the convergence PDF.
	    my $convergencePDF = pdl zeroes($amplificationCount);
	    my $finiteIndices  = which($gridConvergence >= $kappaEmpty);
	    $convergencePDF->index($finiteIndices) .= exp(-(0.5/$omegaKappa**2)*(log(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))+0.5*$omegaKappa**2)**2*(1.0+$aKappa/(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))))/($gridConvergence->index($finiteIndices)+abs($kappaEmpty));
	    my $interpolator = PDL::GSL::INTERP->init('linear',$gridConvergence,$convergencePDF,{Sort => 0}); 
	    my $normalizer   = $interpolator->integ($gridConvergence((0)),$gridConvergence(($amplificationCount-1)),{Extrapolate => 0});
	    $convergencePDF /= $normalizer;
	    
	    # Compute the amplification PDF.
	    my $amplification0          = pdl 3.0;
	    my $convergence0            = 1.0-1.0/sqrt($amplification0);
	    my $convergenceInterpolator = PDL::GSL::INTERP->init('linear',$gridConvergence,$convergencePDF,{Sort => 0}); 
	    my $convergencePDF0         = $interpolator->eval($convergence0,{Extrapolate => 0});
	    my $isAmplified             = which($gridAmplification > 1.0);
	    my $amplificationPDF        = 0.5*(1.0-$gridConvergence)**3*$convergencePDF;
	    $amplificationPDF->index($isAmplified) += 0.5*exp(-0.25/($gridAmplification->index($isAmplified)-1.0)**4)*(1.0-$convergence0)**3*$convergencePDF0/($gridAmplification->index($isAmplified)/$amplification0)**3;
	    my $amplificationInterpolator = PDL::GSL::INTERP->init('linear',$gridAmplification,$amplificationPDF,{Sort => 0}); 
	    my $amplificationNormalizer   = $amplificationInterpolator->integ($gridAmplification((0)),$gridAmplification(($amplificationCount-1)),{Extrapolate => 0});
	    $amplificationPDF /= $amplificationNormalizer;
	    
	    # Compute the cumulative PDF.
	    my $amplificationCDF = pdl zeroes($amplificationCount);
	    $amplificationInterpolator = PDL::GSL::INTERP->init('linear',$gridAmplification,$amplificationPDF,{Sort => 0}); 
	    for(my $m=0;$m<nelem($amplificationPDF);++$m) {
		$amplificationCDF(($m)) .= $amplificationInterpolator->integ($gridAmplification((0)),$gridAmplification(($m)),{Extrapolate => 0});
	    }

	    # Interpolate onto the probability grid.
	    my $uniqueIndices = $amplificationCDF->uniqind();
	    my $probabilityInterpolator = PDL::GSL::INTERP->init('linear',$amplificationCDF->index($uniqueIndices),$gridAmplification->index($uniqueIndices),{Sort => 0}); 
	    $amplificationPDFs(($i+1),:) .= $probabilityInterpolator->eval($gridProbability,{Extrapolate => 1});

	    # Write a message.
	    print "<- done\n";
	}
	
        # Add in a zero redshift amplification PDF which is just a delta function at amplification of unity.
	$amplificationPDFs((0),:)            .= 1.0;
	my $outputRedshift                    = pdl [0.0];
	$gridRedshift                         = $outputRedshift  ->append($gridRedshift);
	my $outputOmegaKappa                  = pdl [0.0];
	$outputOmegaKappa                     = $outputOmegaKappa->append($omegaKappas );
	my $outputAKappa                      = pdl [0.0];
	$outputAKappa                         = $outputAKappa    ->append($aKappas     );

	# Store the computed PDFs to file.
	my $hdfFile = new PDL::IO::HDF5(">".$lensingPdfFile);
	my $pdfDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "amplificationPDF",
	    parent  => $hdfFile,
	    fileObj => $hdfFile
	    );
	$pdfDataSet->set($amplificationPDFs);
	$pdfDataSet->attrSet(
	    'description' => 'The amplification PDF as a function of redshift and amplification.\n'
	    );
	my $redshiftDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "redshift",
	    parent  => $hdfFile,
	    fileObj => $hdfFile
	    );
	$redshiftDataSet->set($gridRedshift);
	$redshiftDataSet->attrSet(
	    'description' => 'Redshifts at which the PDF is tabulated.'
	    );
	my $probabilityDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "probability",
	    parent  => $hdfFile,
	    fileObj => $hdfFile
	    );
	$probabilityDataSet->set($gridProbability);
	$probabilityDataSet->attrSet(
	    'description' => 'Probabilities at which the PDF is tabulated.'
	    );
	# Write the fit parameters to file.
	my $omegaKappaDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "omega_kappa",
	    parent  => $hdfFile,
	    fileObj => $hdfFile
	    );
	$omegaKappaDataSet->set($outputOmegaKappa);
	$omegaKappaDataSet->attrSet(
	    'description' => 'Values of the fitting parameter omega_kappa for each redshift.'
	    );
	my $aKappaDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "A_kappa",
	    parent  => $hdfFile,
	    fileObj => $hdfFile
	    );
	$aKappaDataSet->set($outputAKappa);
	$aKappaDataSet->attrSet(
	    'description' => 'Values of the fitting parameter A_kappa for each redshift.'
	    );
	# Write attributes to the file.
	$hdfFile->attrSet(
	    'description' => 'Cumulative probability density functions for gravitational lensing amplification, computed from the analytic PDF of Takahasi et al. (2011).',
	    'URL'          => 'http://adsabs.harvard.edu/abs/2011arXiv1106.3823T',
	    'source'       => 'Takahasi et al. (2011; arXiv:1106.3823)',
	    'computedDate' => DateTime->now(),
	    'computedBy'   => 'Galacticus::LensingAmplification (part of the Galacticus toolkit)'
	    );
    }

    # Retrieve the redshift of the galaxies.
    &HDF5::Get_Dataset($dataBlock,["redshift"]);

    # Generate uniform random deviates in the 0-1 range for all galaxies.
    my $dataSets = $dataBlock->{'dataSets'};
    my $uniformDeviates = pdl random(nelem($dataSets->{"redshift"}));

    # Interpolate in the probability distributions to get the corresponding magnifications.
    my $interpIndices      = pdl zeroes(2,nelem($dataSets->{"redshift"}));
    my $interpRedshift     = PDL::GSL::INTERP->init('linear',$gridRedshift   ,sequence(nelem($gridRedshift   )),{Sort => 0}); 
    my $interpProbability  = PDL::GSL::INTERP->init('linear',$gridProbability,sequence(nelem($gridProbability)),{Sort => 0}); 
    $interpIndices((0),:) .= $interpRedshift   ->eval($dataSets->{"redshift"},{Extrapolate => 1});
    $interpIndices((1),:) .= $interpProbability->eval($uniformDeviates       ,{Extrapolate => 1});
    $dataSets->{"lensingAmplification"} = $amplificationPDFs->interpND($interpIndices);

}

sub Kappa_Empty_Integrand {
    # Integrand function for computing kappa_empty.
    my ($z) = @_;
    my $a       = 1.0/(1.0+$z);
    my $hubbleZ = $hubble0*sqrt($omegaM0/$a**3+$omegaDE0+(1.0-$omegaM0-$omegaDE0)/$a**2);
    my $comove02Z  = $cosmology->comov_dist($z);
    my $comove02Zs = $cosmology->comov_dist($sourceRedshift);
    my $comoveZ2Zs = $comove02Zs-$comove02Z;
    return (-1.5*$hubble0**2*$omegaM0*(1.0+$z)/$hubbleZ)*$comove02Z*$comoveZ2Zs/$comove02Zs/2.998e5;
}

sub Convergence_RMS {
    # Returns the RMS of the convergence, <kappa^2>^1/2, for
    # lensing. Currently just uses a fit to the results of Takahashi et
    # al. (their Fig. 2) as using an analytic estimate of the nonlinear
    # power spectrum doesn't give very accurate results.

    # Get the redshifts.
    my $redshifts = shift;

    # Set the values of <kappa^2>^1/2 as a function of redshift.
    my $tableRedshifts = pdl ( 0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0 );
    my $tableKappaRMS  = pdl ( 0.0, 0.02297, 0.04405, 0.0600, 0.07919, 0.1011, 0.1159);

    # Interpolate in redshift to get the required <kappa^2>^1/2 values.
    (my $kappaRMS, my $error) = interpolate($redshifts,$tableRedshifts,$tableKappaRMS);

    # Return the computed values.
    return $kappaRMS;
}

sub Fit_Parameters {
    # Used in fitting parameters of the convergence distribution to the measured RMS convergence.

    # Get the arguments.
    my ($x,$par,$ym,$dyda) = @_;

    # $omegaKappa and $aKappa are fit parameters, internal to this function.
    my ($omegaKappa,$aKappa) = map { $par->slice("($_)") } (0..1);
    $omegaKappa = 10.0**$omegaKappa;
    $aKappa     = 10.0**$aKappa;

    # Create grid of convergences.
    my $convergenceLow   = $kappaEmpty;
    my $convergenceHigh  = 4.0;
    my $convergenceCount = 10000;
    my $gridConvergence  = pdl sequence($convergenceCount/2)*(-$convergenceLow-$convergenceLow)/($convergenceCount/2-1)+$convergenceLow;
    $gridConvergence     = $gridConvergence->append((sequence($convergenceCount/2)+1)*($convergenceHigh+$convergenceLow)/($convergenceCount/2)-$convergenceLow);

    # Compute the distribution function.
    my $convergencePDF            = pdl zeroes($convergenceCount);
    my $convergencePDFDomegaKappa = pdl zeroes($convergenceCount);
    my $convergencePDFDaKappa     = pdl zeroes($convergenceCount);
    my $finiteIndices             = which($gridConvergence > $kappaEmpty);
    $convergencePDF->index($finiteIndices) .= exp(-(0.5/$omegaKappa**2)*(log(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))+0.5*$omegaKappa**2)**2*(1.0+$aKappa/(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))))/($gridConvergence->index($finiteIndices)+abs($kappaEmpty));
    my $interpolator = PDL::GSL::INTERP->init('linear',$gridConvergence,$convergencePDF,{Sort => 0}); 
    my $normalizer   = $interpolator->integ($gridConvergence((0)),$gridConvergence(($convergenceCount-1)),{Extrapolate => 0});
    $convergencePDF /= $normalizer if ( $normalizer > 0.0 );
    $convergencePDFDomegaKappa->index($finiteIndices) .= 
	$convergencePDF->index($finiteIndices)
	*(1.0+$aKappa/(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty)))
	*(
	    +(log(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))+0.5*$omegaKappa**2)**2/$omegaKappa**3
	    -(log(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))+0.5*$omegaKappa**2)   /$omegaKappa
	);
    $convergencePDFDaKappa->index($finiteIndices) .= 
	$convergencePDF->index($finiteIndices)
	*(-0.5)
	*(log(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty))+0.5*$omegaKappa**2)**2/$omegaKappa**2
	/(1.0+$gridConvergence->index($finiteIndices)/abs($kappaEmpty));


    # Compute moments of the distribution.
    my (@dy) = map {$dyda -> slice(",($_)") } (0..1);
    my $moments            = pdl zeroes(2);
    my $momentsDomegaKappa = pdl zeroes(2);
    my $momentsDaKappa     = pdl zeroes(2);
    for(my $i=0;$i<nelem($x);++$i) {
	my $integrand                   = $convergencePDF           *($gridConvergence**$x->index($i));
	my $integrandDomegaKappa        = $convergencePDFDomegaKappa*($gridConvergence**$x->index($i));
	my $integrandDaKappa            = $convergencePDFDaKappa    *($gridConvergence**$x->index($i));
	my $interpolator                = PDL::GSL::INTERP->init('linear',$gridConvergence,$integrand           ,{Sort => 0}); 
	my $interpolatorDomegaKappa     = PDL::GSL::INTERP->init('linear',$gridConvergence,$integrandDomegaKappa,{Sort => 0}); 
	my $interpolatorDaKappa         = PDL::GSL::INTERP->init('linear',$gridConvergence,$integrandDaKappa    ,{Sort => 0}); 
	my $moment                      = $interpolator           ->integ($gridConvergence((0)),$gridConvergence(($convergenceCount-1)),{Extrapolate => 0});    
	my $momentDomegaKappa           = $interpolatorDomegaKappa->integ($gridConvergence((0)),$gridConvergence(($convergenceCount-1)),{Extrapolate => 0});    
	my $momentDaKappa               = $interpolatorDaKappa    ->integ($gridConvergence((0)),$gridConvergence(($convergenceCount-1)),{Extrapolate => 0});    
	$moments           ->index($i) .= $moment;
	$momentsDomegaKappa->index($i) .= $momentDomegaKappa;
	$momentsDaKappa    ->index($i) .= $momentDaKappa;
    }
    $ym    .= $moments;
    $dy[0] .= $momentsDomegaKappa;
    $dy[1] .= $momentsDaKappa;
}
