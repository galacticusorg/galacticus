#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use XML::Simple;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Cloudy;
use Clone 'clone';
use Data::Dumper;

# Initialize a globally-visible object in which our line data will be stored.
our $lineData;

# Get arguments.
my %options =
    (
     workspace => "cloudyTable/",
     reprocess => "no"
    );
# Parse options.
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Create workspace.
$options{'workspace'} .= "/"
    unless ( $options{'workspace'} =~ m/\/$/ );
system("mkdir -p ".$options{'workspace'});

# Initalize Cloudy.
(my $cloudyPath, my $cloudyVersion) = &Cloudy::Initialize(%options);

# Set physical constants.
my $plancksConstant          = pdl 6.6260700400000e-34; # J s
my $speedOfLight             = pdl 2.9979245800000e+08; # m s¯¹
my $electronVolt             = pdl 1.6021766340000e-19; # J
my $rydbergEnergy            = pdl 1.3605693122994e+01; # eV
my $angstroms                = pdl 1.0000000000000e-10; # m
my $luminositySolar          = pdl 3.8390000000000e+26; # W
my $wavelengthLymanContinuum = pdl 9.1176000000000e+02; # Å

# Specify abundances and depletion model. This is based upon the work by Gutkin, Charlot & Bruzual (2016;
# https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G).
my %abundancesReference =
    (
     H =>
     {
	 atomicNumber         =>   1,
	 logAbundanceByNumber =>   0.00,
	 undepletedFraction   =>   1.000
     },
     He =>
     {
	 atomicNumber         =>   2,
	 logAbundanceByNumber =>  -1.01,
	 undepletedFraction   =>   1.000,
	 adjustAbundance      => \&adjustAbundanceHelium
     },
     Li =>
     {
	 atomicNumber         =>   3,
	 logAbundanceByNumber => -10.99,
	 undepletedFraction   =>   0.160
     },
     Be =>
     {
	 atomicNumber         =>   4,
	 logAbundanceByNumber => -10.63,
	 undepletedFraction   =>   0.600
     },
     B =>
     {
	 atomicNumber         =>   5,
	 logAbundanceByNumber =>  -9.47,
	 undepletedFraction   =>   0.130
     },
     C =>
     {
	 atomicNumber         =>   6,
	 logAbundanceByNumber =>  -3.53,
	 undepletedFraction   =>   0.500,
	 adjustAbundance      => \&adjustAbundanceCarbon
     },
     N =>
     {
	 atomicNumber         =>   7,
	 logAbundanceByNumber =>  -4.32,
	 undepletedFraction   =>   1.000,
	 adjustAbundance      => \&adjustAbundanceNitrogen
     },
     O =>
     {
	 atomicNumber         =>   8,
	 logAbundanceByNumber =>  -3.17,
	 undepletedFraction   =>   0.700
     },
     F =>
     {
	 atomicNumber         =>   9,
	 logAbundanceByNumber =>  -7.47,
	 undepletedFraction   =>   0.300
     },
     Ne =>
     {
	 atomicNumber         =>  10,
	 logAbundanceByNumber =>  -4.01,
	 undepletedFraction   =>   1.000
     },
     Na =>
     {
	 atomicNumber         =>  11,
	 logAbundanceByNumber =>  -5.70,
	 undepletedFraction   =>   0.250
     },
     Mg =>
     {
	 atomicNumber         =>  12,
	 logAbundanceByNumber =>  -4.45,
	 undepletedFraction   =>   0.200
     },
     Al =>
     {
	 atomicNumber         =>  13,
	 logAbundanceByNumber =>  -5.56,
	 undepletedFraction   =>   0.020
     },
     Si =>
     {
	 atomicNumber         =>  14,
	 logAbundanceByNumber =>  -4.48,
	 undepletedFraction   =>   0.100
     },
     P =>
     {
	 atomicNumber         =>  15,
	 logAbundanceByNumber =>  -6.57,
	 undepletedFraction   =>   0.250
     },
     S =>
     {
	 atomicNumber         =>  16,
	 logAbundanceByNumber =>  -4.87,
	 undepletedFraction   =>   1.000
     },
     Cl =>
     {
	 atomicNumber         =>  17,
	 logAbundanceByNumber =>  -6.53,
	 undepletedFraction   =>   0.500
     },
     Ar =>
     {
	 atomicNumber         =>  18,
	 logAbundanceByNumber =>  -5.63,
	 undepletedFraction   =>   1.000
     },
     K =>
     {
	 atomicNumber         =>  19,
	 logAbundanceByNumber =>  -6.92,
	 undepletedFraction   =>   0.300
     },
     Ca =>
     {
	 atomicNumber         =>  20,
	 logAbundanceByNumber =>  -5.67,
	 undepletedFraction   =>   0.003
     },
     Sc =>
     {
	 atomicNumber         =>  21,
	 logAbundanceByNumber =>  -8.86,
	 undepletedFraction   =>   0.005
     },
     Ti =>
     {
	 atomicNumber         =>  22,
	 logAbundanceByNumber =>  -7.01,
	 undepletedFraction   =>   0.008
     },
     V =>
     {
	 atomicNumber         =>  23,
	 logAbundanceByNumber =>  -8.03,
	 undepletedFraction   =>   0.006
     },
     Cr =>
     {
	 atomicNumber         =>  24,
	 logAbundanceByNumber =>  -6.36,
	 undepletedFraction   =>   0.006
     },
     Mn =>
     {
	 atomicNumber         =>  25,
	 logAbundanceByNumber =>  -6.64,
	 undepletedFraction   =>   0.050
     },
     Fe =>
     {
	 atomicNumber         =>  26,
	 logAbundanceByNumber =>  -4.51,
	 undepletedFraction   =>   0.010
     },
     Co =>
     {
	 atomicNumber         =>  27,
	 logAbundanceByNumber =>  -7.11,
	 undepletedFraction   =>   0.010
     },
     Ni =>
     {
	 atomicNumber         =>  28,
	 logAbundanceByNumber =>  -5.78,
	 undepletedFraction   =>   0.040
     },
     Cu =>
     {
	 atomicNumber         =>  29,
	 logAbundanceByNumber =>  -7.82,
	 undepletedFraction   =>   0.100
     },
     Zn =>
     {
	 atomicNumber         =>  30,
	 logAbundanceByNumber =>  -7.43,
	 undepletedFraction   =>   0.250
     }
    );

# Augment abundance data with names and atomic masses. These are read from Galacticus' atomic data file.
my $xml = new XML::Simple();
my $atomicData = $xml->XMLin($ENV{'GALACTICUS_DATA_PATH'}."/static/abundances/Atomic_Data.xml");
foreach my $elementName ( keys(%{$atomicData->{'element'}}) ) {
    my $element = $atomicData->{'element'}->{$elementName};
    if ( exists($abundancesReference{$element->{'shortLabel'}}) ) {
	# Cloudy uses British spelling (presumably because it originated at Cambridge) - so correct for that here so that we can
	# use the correct names when generating a Cloudy script.
	my $nameBritish = $elementName eq "Sulfur" ? "Sulphur" : $elementName;
	# Store the name and atomic mass for this element.
	$abundancesReference{$element->{'shortLabel'}}->{'atomicMass'} = $element->{'atomicMass'};
	$abundancesReference{$element->{'shortLabel'}}->{'name'      } = $nameBritish;
    }
}

# Get an SSP model.
## Read data from file.
my $stellarPopulationsFileName  = "SSP_Spectra_BC2003_highResolution_imfChabrier.hdf5";
my $stellarPopulations          = new PDL::IO::HDF5($stellarPopulationsFileName);
my $ages                        = $stellarPopulations->dataset('ages'         )->get();
my $logMetallicities            = $stellarPopulations->dataset('metallicities')->get();
my $wavelength                  = $stellarPopulations->dataset('wavelengths'  )->get();
my $spectra                     = $stellarPopulations->dataset('spectra'      )->get();
## Determine energies in Rydbergs.
my $energy                      = $plancksConstant*$speedOfLight/$wavelength/$angstroms/$electronVolt/$rydbergEnergy;
## Construct wavelength intervals.
my $deltaWavelength             = $wavelength->copy();
$deltaWavelength->(0:-2)       .= $wavelength->(1:-1)-$deltaWavelength->(0:-2);
$deltaWavelength->((-1))       .=                     $deltaWavelength->((-2));

# Refine the range of metallicities, by interpolating the spectra to intermediate points.
my $refineMetallicityBy       = 2;
my $countRefinedMetallicities = $refineMetallicityBy*(nelem($logMetallicities)-1)+1;
my $logMetallicitiesRefined   = pdl zeros(                                  $countRefinedMetallicities);
my $spectraRefined            = pdl zeros($spectra->dim(0),$spectra->dim(1),$countRefinedMetallicities);
for(my $i=0;$i<nelem($logMetallicities)-1;++$i) {
    my $deltaLogMetallicity = $logMetallicities->(($i+1))-$logMetallicities->(($i));
    for(my $j=0;$j<$refineMetallicityBy;++$j) {
	(my $zero, my $nonZero) =
	    which_both(
		($spectra->(:,:,($i  )) <= 0.0)
		|
		($spectra->(:,:,($i+1)) <= 0.0)
	    );
	# Spectra are interpolated in log-log space where possible (i.e. where the spectrum is non-zero at both endpoints of the
	# interpolation), and in lin-log space otherwise. Metallicities are always interpolated in log space.
	$spectraRefined         ->(:,:,($i*$refineMetallicityBy+$j))->flat()->(   $zero) .=     +$spectra            ->(:,:,($i  ))->flat()->(   $zero)       *(1.0-($j/$refineMetallicityBy))
	                                                                                        +$spectra            ->(:,:,($i+1))->flat()->(   $zero)       *(0.0+($j/$refineMetallicityBy));
	$spectraRefined         ->(:,:,($i*$refineMetallicityBy+$j))->flat()->($nonZero) .= exp(
	  	                                                                                +$spectra            ->(:,:,($i  ))->flat()->($nonZero)->log()*(1.0-($j/$refineMetallicityBy))
	 	                                                                                +$spectra            ->(:,:,($i+1))->flat()->($nonZero)->log()*(0.0+($j/$refineMetallicityBy))
	    );
	$logMetallicitiesRefined->(    ($i*$refineMetallicityBy+$j))                     .=     +$logMetallicities   ->(    ($i  ))
	                                                                                        +$deltaLogMetallicity                                         *     ($j/$refineMetallicityBy) ;
    }
}
$logMetallicitiesRefined->(    (-1)) .= $logMetallicities->(    (-1));    
$spectraRefined         ->(:,:,(-1)) .= $spectra         ->(:,:,(-1));    
undef($logMetallicities);
undef($spectra         );
$logMetallicities = $logMetallicitiesRefined;
$spectra          = $spectraRefined         ;

# Evaluate the number of Lyman continuum photons emitted per second for each population (age,metallicity).
## Construct the integrand. Spectra are in units of L☉ Hz¯¹. We want to evaluate Qₕ = ∫_Eₕ^∞ dν S_ν/hν = ∫₀^λₕ dλ S_ν / hλ.
my $integrand                   = $spectra*$deltaWavelength/$wavelength*$luminositySolar/$plancksConstant;
## Find the range to include in each integral.
my $hydrogenContinuum           = which($wavelength < $wavelengthLymanContinuum);
## Evaluate the integrals.
my $ionizingLuminosityPerMass   = $integrand->($hydrogenContinuum,:,:)->sumover();

# Define ionizing luminosities, Qₕ, to tabulate.
my $logHydrogenLuminosities = pdl [ 48.0, 49.0, 50.0, 51.0, 52.0 ];

# Define hydrogen densities, nₕ, to tabulate.
my $logHydrogenDensities    = pdl [  1.0,  1.5,  2.0,  2.5,  3.0 ];

# Define the list of lines to extract. The following dictionary contains keys which match the line names in the Cloudy emission
# lines output file, and values which are our internal names for these lines.
my %lineList =
    (
     ## WORKAROUND     
     ## Currently using in-air wavelengths as Cloudy has a bug that prevents us from outputting in vacuum wavelengths. See these
     ## two error reports:
     ##   https://cloudyastrophysics.groups.io/g/Main/topic/101921207#5396
     ##   https://cloudyastrophysics.groups.io/g/Main/topic/102424985#5431
     # In air wavelengths.
     "H  1                6562.80A" => "balmerAlpha6565"  ,
     "H  1                4861.32A" => "balmerBeta4863"   ,
     "H  1                4340.46A" => "balmerGamma4342"  ,
     "H  1                4101.73A" => "balmerDelta4103"  ,
     "H  1                1.87510m" => "paschenAlpha18756",
     "H  1                1.28181m" => "paschenBeta12822" ,
     "O  2                3726.03A" => "oxygenII3727"     ,
     "O  2                3728.81A" => "oxygenII3730"     ,
     "O  3                4958.91A" => "oxygenIII4960"    ,
     "O  3                5006.84A" => "oxygenIII5008"    ,
     "O  3                4931.23A" => "oxygenIII4933"    ,
     "N  2                6583.45A" => "nitrogenII6585"   ,
     "N  2                6548.05A" => "nitrogenII6550"   ,
     "S  2                6730.82A" => "sulfurII6733"     ,
     "S  2                6716.44A" => "sulfurII6718"     ,
     "S  3                9068.62A" => "sulfurIII9071"    ,
     "S  3                9530.62A" => "sulfurIII9533"
     # In vacuum wavelengths.
     # "H  1                6564.62A" => "balmerAlpha6565"  ,
     # "H  1                4862.69A" => "balmerBeta4863"   ,
     # "H  1                4341.68A" => "balmerGamma4342"  ,
     # "H  1                4102.89A" => "balmerDelta4103"  ,
     # "H  1                1.87561m" => "paschenAlpha18756",
     # "H  1                1.28215m" => "paschenBeta12822" ,
     # "O  2                3727.09A" => "oxygenII3727"     ,
     # "O  2                3729.88A" => "oxygenII3730"     ,
     # "O  3                4960.29A" => "oxygenIII4960"    ,
     # "O  3                5008.24A" => "oxygenIII5008"    ,
     # "O  3                4932.60A" => "oxygenIII4933"    ,
     # "N  2                6585.27A" => "nitrogenII6585"   ,
     # "N  2                6549.86A" => "nitrogenII6550"   ,
     # "S  2                6732.67A" => "sulfurII6733"     ,
     # "S  2                6718.29A" => "sulfurII6718"     ,
     # "S  3                9071.11A" => "sulfurIII9071"    ,
     # "S  3                9533.23A" => "sulfurIII9533"
    );


if ( $options{'reprocess'} eq "yes" ) {
    # Reprocess output files. This can be useful if some previous processing of Cloudy output files failed (we often have tens of
    # thousands of these so some intermittment failures can occur).
    my $tableFile         = new PDL::IO::HDF5("cloudyTable.hdf5");
    my $lineGroup         = $tableFile->group('lines');
    $lineData->{'status'} = $lineGroup->dataset('status')->get();
    foreach my $lineIdentifier ( keys(%lineList) ) {
	my $lineName = $lineList{$lineIdentifier};
	$lineData->{$lineName}->{'luminosity'} = $lineGroup->dataset($lineName)->get();
    }
    my $jobNumber = -1;
    for(my $iAge=0;$iAge<$ionizingLuminosityPerMass->dim(0);++$iAge) {
	for(my $iMetallicity=0;$iMetallicity<$ionizingLuminosityPerMass->dim(1);++$iMetallicity) {
	    for(my $iLogHydrogenLuminosity=0;$iLogHydrogenLuminosity<nelem($logHydrogenLuminosities);++$iLogHydrogenLuminosity) {
		for(my $iLogHydrogenDensity=0;$iLogHydrogenDensity<nelem($logHydrogenDensities);++$iLogHydrogenDensity) {
		    ++$jobNumber;
		    next
			if ( $lineData->{'status'}->(($iAge),($iMetallicity),($iLogHydrogenLuminosity),($iLogHydrogenDensity)) == 0 );
		    # Reset the status before attempting to reprocess.
		    my $statusOld = $lineData->{'status'}->(($iAge),($iMetallicity),($iLogHydrogenLuminosity),($iLogHydrogenDensity))->sclr();
		    $lineData->{'status'}->(($iAge),($iMetallicity),($iLogHydrogenLuminosity),($iLogHydrogenDensity)) .= 0;
		    &linesParse(
			"lines"        .$jobNumber.".out",
			"continuum"    .$jobNumber.".out",
			"emissionLines".$jobNumber.".sh" ,
			"emissionLines".$jobNumber.".log",
			"cloudyInput"  .$jobNumber.".txt",
			$iAge                            ,
			$iMetallicity                    ,
			$iLogHydrogenLuminosity          ,
			$iLogHydrogenDensity
			);
		    print "Reprocess job number ".$jobNumber." (status = ".$statusOld." ==> ".$lineData->{'status'}->(($iAge),($iMetallicity),($iLogHydrogenLuminosity),($iLogHydrogenDensity)).")\n";
		}
	    }
	}
    }
} else {
    # Compute the table.
    # Initialize the line luminosity tables.
    $lineData->{$lineList{$_}}->{'luminosity'} =
	pdl zeroes
	(
	 nelem($ages                   ),
	 nelem($logMetallicities       ),
	 nelem($logHydrogenLuminosities),
	 nelem($logHydrogenDensities   )
	)
	foreach ( keys(%lineList) );
    $lineData->{'status'} =
	pdl long zeroes
	(
	 nelem($ages                   ),
	 nelem($logMetallicities       ),
	 nelem($logHydrogenLuminosities),
	 nelem($logHydrogenDensities   )
	);
    # Begin constructing Cloudy jobs.
    ## Iterate over (age,metallicity).
    my @jobs;
    my $jobNumber = -1;
    my $jobCount  = +      $ionizingLuminosityPerMass->dim(0)
                    *      $ionizingLuminosityPerMass->dim(1)
                    *nelem($logHydrogenLuminosities          )
                    *nelem($logHydrogenDensities             );
    for(my $iAge=0;$iAge<$ionizingLuminosityPerMass->dim(0);++$iAge) {
	for(my $iMetallicity=0;$iMetallicity<$ionizingLuminosityPerMass->dim(1);++$iMetallicity) {
	    # Normalize the spectrum - this is a convenience only as the normalization will be recomputed by Cloudy.
	    my $normalizer = $spectra->(:,($iAge),($iMetallicity))->maximum();
	    $spectra->(:,($iAge),($iMetallicity)) /= $normalizer;	
	    # Get abundances for this metallicity.
	    ## Compute metallicity relative to Solar.
	    my $metallicity        = pdl 10.0**$logMetallicities->(($iMetallicity));
	    ## Specify a dust-to-metals ratio. The following corresponds to the dust-to-metals ratio for the reference model of
	    # Gutkin, Charlot & Bruzual (2016; https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; table 1). It differs slightly
	    # from their stated value, but agrees with out internal calculation of this value from their data.
	    my $dustToMetals       = pdl 0.401;
	    my %abundances         = &adjustAbundances(\%abundancesReference,$metallicity,$dustToMetals);
	    # Create a depletion file that can be read by Cloudy. Cloudy does not permit digits in these file names. To work around
	    # this, we translate our numerical age and metallicity indices into ASCII characters, by mapping digits (0→A, 1→B,
	    # etc.).
	    my $encodedAge         = join("",map {chr($_+17)} unpack("c*","$iAge"        ));
	    my $encodedMetallicity = join("",map {chr($_+17)} unpack("c*","$iMetallicity"));
	    my $depletionFileName = "grains_".$encodedAge."_".$encodedMetallicity.".dpl";
	    open(my $depletionFile,">",$options{'workspace'}.$depletionFileName);
	    foreach my $element ( sort(keys(%abundances)) ) {
		print $depletionFile $abundances{$element}->{'name'}." ".$abundances{$element}->{'undepletedFraction'}."\n";
	    }
	    close($depletionFile);
	    ## Iterate over Qₕ.
	    for(my $iLogHydrogenLuminosity=0;$iLogHydrogenLuminosity<nelem($logHydrogenLuminosities);++$iLogHydrogenLuminosity) {
		## Iterate over nₕ.
		for(my $iLogHydrogenDensity=0;$iLogHydrogenDensity<nelem($logHydrogenDensities);++$iLogHydrogenDensity) {
		    # Generate a Cloudy parameter file.
		    ++$jobNumber;
		    print "Generating model ".$jobNumber." of ".$jobCount."\n"
			if ( $jobNumber % 100 == 0 );
		    my $cloudyScript;
		    $cloudyScript .= "title emission line job number ".$jobNumber."\n";
		    $cloudyScript .= "# [".$iAge                  ."] age     = ".$ages                   ->(($iAge                  ))."\n";
		    $cloudyScript .= "# [".$iMetallicity          ."] log Z   = ".$logMetallicities       ->(($iMetallicity          ))."\n";
		    $cloudyScript .= "# [".$iLogHydrogenLuminosity."] log Q_H = ".$logHydrogenLuminosities->(($iLogHydrogenLuminosity))."\n";
		    $cloudyScript .= "# [".$iLogHydrogenDensity   ."] log n_H = ".$logHydrogenDensities   ->(($iLogHydrogenDensity   ))."\n";
		    ## Set the input spectrum for Cloudy.
		    my $counter = -1;
		    for(my $iWavelength=nelem($wavelength)-1;$iWavelength>=0;--$iWavelength) {
			++$counter;
			if ( $counter % 3 == 0 ) {
			    if ( $counter == 0 ) {
				$cloudyScript .=  "interpolate";
			    } else {
				$cloudyScript .=  "\ncontinue";
			    }
			}
			$cloudyScript .=  " (".$energy->(($iWavelength))." ".$spectra->(($iWavelength),($iAge),($iMetallicity))->log10().")";
		    }
		    $cloudyScript .= "\n";
		    ## Set normalization of the spectrum.
		    $cloudyScript .= "q(h) = ".$logHydrogenLuminosities->(($iLogHydrogenLuminosity))."\n";
		    # Set the chemical composition of the HII region.
		    # Start with Cloudy's HII region abundances - but no grains - we will add these later.
		    $cloudyScript .= "abundances HII region no grains\n";
		    # Set individual element abundances.
		    foreach my $element ( sort(keys(%abundances)) ) {
			# Skip hydrogen.
			next
			    if ( $abundances{$element}->{'atomicNumber'} < 2 );
			# Set the abundance for this element.
			$cloudyScript .= "element abundances ".lc($abundances{$element}->{'name'})." ".$abundances{$element}->{'logAbundanceByNumber'}."\n";
		    }
		    # Specify Cloudy's default ISM grains, but with abundance reduced in proportion to the metallicity to retain a fixed
		    # dust-to-metals ratio.
		    $cloudyScript .= "grains _ISM ".$logMetallicities->(($iMetallicity))." _log\n";
		    # Deplete metals into grains using our custom depletions file. Note that these depletions differ from those assumed by
		    # the "grains _ISM" model above. This seems to be at the ~10% level, so we do not worry too much about this.
		    $cloudyScript .= "metals deplete \"".$depletionFileName."\"\n";
		    # Set HII region density - this is log₁₀(nₕ/cm¯³).
		    $cloudyScript .= "hden ".$logHydrogenDensities   ->(($iLogHydrogenDensity   ))."\n";
		    # Set other HII region properties.
		    $cloudyScript .= "sphere expanding\n";
		    $cloudyScript .= "radius 16.0\n";
		    # Set cosmic rays (needed to avoid problems in Cloudy in neutral gas).
		    $cloudyScript .= "cosmic rays background\n";
		    # Set stopping criteria.
		    $cloudyScript .= "stop temperature 1000 k\n";
		    $cloudyScript .= "iterate to convergence\n";
		    ## Output the continuum for refernce.
		    $cloudyScript .= "punch continuum \"continuum".$jobNumber.".out\"\n";
		    ## Set line output options.
		    $cloudyScript .= "print lines faint _off\n";
		    ## WORKAROUND     
		    ## Currently disabled as Cloudy has a bug that prevents us from outputting in vacuum wavelengths. See these
		    ## two error reports:
		    ##   https://cloudyastrophysics.groups.io/g/Main/topic/101921207#5396
		    ##   https://cloudyastrophysics.groups.io/g/Main/topic/102424985#5431
		    #$cloudyScript .= "print line vacuum\n";
		    $cloudyScript .= "save lines, array \"lines".$jobNumber.".out\"\n";
		    ## Write the Cloudy script to file.
		    my $cloudyScriptFileName = "cloudyInput".$jobNumber.".txt";
		    open(my $cloudyScriptFile,">".$options{'workspace'}.$cloudyScriptFileName);
		    print $cloudyScriptFile $cloudyScript;
		    close($cloudyScriptFile);
		    ## Construct a job to run this parameter file.
		    my %job = 
			(
			 launchFile   => $options{'workspace'}."emissionLines".$jobNumber.".sh" ,
			 label        =>                       "emissionLines".$jobNumber       ,
			 logFile      => $options{'workspace'}."emissionLines".$jobNumber.".log",
			 command      => 
		                         "cd ".$options{'workspace'}."; ulimit -c 0\n"                 .
		                         $cloudyPath."/source/cloudy.exe < ".$cloudyScriptFileName."\n".
		                         "if [ $? != 0 ]; then\necho CLOUDY FAILED\nfi\n"              ,
			 ppn          => 1,
			 walltime     => "10:00:00",
			 onCompletion => 
			 {
			     function  => \&linesParse,
			     arguments => 
				 [
				  "lines"        .$jobNumber.".out",
				  "continuum"    .$jobNumber.".out",
				  "emissionLines".$jobNumber.".sh" ,
				  "emissionLines".$jobNumber.".log",
				  $cloudyScriptFileName            ,
				  $iAge                            ,
				  $iMetallicity                    ,
				  $iLogHydrogenLuminosity          ,
				  $iLogHydrogenDensity
				 ]
			 }
			);
		    ## Push job to job list.
		    push(@jobs,\%job);
		}
	    }
	}
    }
    # Launch all jobs.
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs);
}

# Validate results.
&validateResults($ages,$logMetallicities,$logHydrogenLuminosities,$logHydrogenDensities,$ionizingLuminosityPerMass,$lineData);

# Write the line data to file.
my $tableFile = new PDL::IO::HDF5(">cloudyTable.hdf5");
# Write parameter grid points and attributes.
$tableFile->dataset('age'                                 )->    set(               $ages                                                           );
$tableFile->dataset('age'                                 )->attrSet(description => "Age of the stellar population."                                );
$tableFile->dataset('age'                                 )->attrSet(units       => "Gyr"                                                           );
$tableFile->dataset('age'                                 )->attrSet(unitsInSI   => 3.15576e16                                                      );
$tableFile->dataset('metallicity'                         )->    set(               10.0**$logMetallicities                                         );
$tableFile->dataset('metallicity'                         )->attrSet(description => "Metallicity relative to Solar."                                );
$tableFile->dataset('ionizingLuminosityHydrogen'          )->    set(               10.0**$logHydrogenLuminosities                                  );
$tableFile->dataset('ionizingLuminosityHydrogen'          )->attrSet(description => "Hydrogen ionizing photon emission rate."                       );
$tableFile->dataset('ionizingLuminosityHydrogen'          )->attrSet(units       => "photons s¯¹"                                                   );
$tableFile->dataset('ionizingLuminosityHydrogen'          )->attrSet(unitsInSI   => 1.0                                                             );
$tableFile->dataset('densityHydrogen'                     )->    set(               10.0**$logHydrogenDensities                                     );
$tableFile->dataset('densityHydrogen'                     )->attrSet(description => "Hydrogen density."                                             );
$tableFile->dataset('densityHydrogen'                     )->attrSet(units       => "cm¯³"                                                          );
$tableFile->dataset('densityHydrogen'                     )->attrSet(unitsInSI   => 1.0e6                                                           );

# Write table of ionizing rates per unit mass of stars formed.
$tableFile->dataset('ionizingLuminosityHydrogenNormalized')->    set(               $ionizingLuminosityPerMass                                      );
$tableFile->dataset('ionizingLuminosityHydrogenNormalized')->attrSet(description => "Hydrogen ionizing photon emission rate per unit mass of stars.");
$tableFile->dataset('ionizingLuminosityHydrogenNormalized')->attrSet(units       => "photons s¯¹ M☉¯¹"                                              );
$tableFile->dataset('ionizingLuminosityHydrogenNormalized')->attrSet(unitsInSI   => 5.025125628140704e-31                                           );

# Write line data.
my $lineGroup = $tableFile->group('lines');
$lineGroup->dataset('status')->set($lineData->{'status'});
$lineGroup->dataset('status')->attrSet(description => "Cloudy model status: 0 = success; 1 = disaster; 2 = non-zero exit status; 3 = missing output file; 4 = missing emission lines");
foreach ( keys(%lineList) ) {
    my $lineName = $lineList{$_};
    $lineGroup->dataset($lineName)->    set(               $lineData->{$lineName}->{'luminosity'});
    $lineGroup->dataset($lineName)->attrSet(description => "Luminosity of the line."             );
    $lineGroup->dataset($lineName)->attrSet(units       => "erg s¯¹"                             );
    $lineGroup->dataset($lineName)->attrSet(unitsInSI   => 1.0e-7                                );
    $lineGroup->dataset($lineName)->attrSet(wavelength  => $lineData->{$lineName}->{'wavelength'});
}

exit;

sub linesParse {
    # Parse output from a Cloudy job to extract line data.
    my $linesFileName          = shift();
    my $continuumFileName      = shift();
    my $launchFileName         = shift();
    my $logFileName            = shift();
    my $cloudyScriptFileName   = shift();
    my $iAge                   = shift();
    my $iMetallicity           = shift();
    my $iLogHydrogenLuminosity = shift();
    my $iLogHydrogenDensity    = shift();
    # Check for successful completion.
    my $status = $lineData
	->{'status'}
        ->(
	   ($iAge                  ),
	   ($iMetallicity          ),
	   ($iLogHydrogenLuminosity),
	   ($iLogHydrogenDensity   )
	);
    my $label = $iAge." ".$iMetallicity." ".$iLogHydrogenLuminosity." ".$iLogHydrogenDensity;
    system("grep -q DISASTER ".$options{'workspace'}.$logFileName);
    if ( $? == 0 ) {
	print "FAIL (Cloudy failed disasterously): ".$label." see ".$logFileName."\n";
	$status .= 1
	    if (
		$status == 0
	    );
    }
    system("grep -q FAILED ".$options{'workspace'}.$logFileName);
    if ( $? == 0 ) {
	print "FAIL [Cloudy exited with non-zero status]: ".$label." see ".$logFileName."\n";
    $status .= 2
	if (
	    $status == 0 
	);
    }
    # Allow multiple attempts to read the lines file in case the file is written over NFS and we must wait for it to catch up.    
    my $attemptsMaximum = 10;
    my @linesFound;
    for(my $attempt=0;$attempt<$attemptsMaximum;++$attempt) {
	# Read the lines file.
	my $badFile = 0;
	open(my $linesFile,$options{'workspace'}.$linesFileName);
	while ( my $line = <$linesFile> ) {
	    unless ( $line =~ m/^#/ ) {
		my @columns        = split(/\t/,$line);
		$badFile = 1
		    unless ( scalar(@columns) == 5 );
		my $lineEnergy     = $columns[0];
		my $lineLabel      = $columns[1];
		my $lineLuminosity = $columns[3];
		if ( exists($lineList{$lineLabel}) ) {
		    my $lineName = $lineList{$lineLabel};
		    $lineLuminosity =~ s/^\s+//;
		    $lineLuminosity =~ s/\s+$//;
		    $lineEnergy     =~ s/^\s+//;
		    $lineEnergy     =~ s/\s+$//;
		    my $lineWavelength = $plancksConstant*$speedOfLight/$rydbergEnergy/$electronVolt/$lineEnergy/$angstroms;
		    $lineData
			->{$lineName}
		        ->{'luminosity'}
		        ->(
		    	    ($iAge                  ),
			    ($iMetallicity          ),
			    ($iLogHydrogenLuminosity),
			    ($iLogHydrogenDensity   )
			    )
			.= 10.0**$lineLuminosity;
		    $lineData
			->{$lineName}
		        ->{'wavelength'}
		        = $lineWavelength;
		    push(@linesFound,$lineName);
		}
	    }
	}
	close($linesFile);
	if ( $badFile == 0 ) {
	    last;
	} else {
	    if ( $attempt == $attemptsMaximum-1 ) {
		print "FAIL [unable to find Cloudy output lines file]: ".$label." see ".$logFileName."\n";
		$status .= 3
		    if (
			$status == 0
		    );
	    } else {
		sleep(10);
	    }
	}
    }
    # Check that we found all lines.
    foreach my $lineName ( keys(%lineList) ) {	
	unless ( grep {$_ eq $lineList{$lineName}} @linesFound ) {
	    print "FAIL [some emission lines missing]: ".$label." '".$lineName."' see ".$logFileName."\n";
	    $status .= 4
		if (
		    $status == 0
		);
	}
    }
    # # Clean up.
    unlink($options{'workspace'}.$linesFileName,$options{'workspace'}.$continuumFileName,$options{'workspace'}.$launchFileName,$options{'workspace'}.$logFileName,$options{'workspace'}.$cloudyScriptFileName)
      	if (
      	    $status == 0
      	);
}

sub adjustAbundances {
    # Copy and modify the provided reference abundances to create abundances suitable for the given metallicity (linear, relative
    # to Solar), and dust-to-metals ratio.
    my %abundancesReference = %{shift()};
    my $metallicity         =   shift() ;
    my $dustToMetalsRatio   =   shift() ;
    # Make a copy.
    my %abundances          = %{clone(\%abundancesReference)};
    # Determine metallicity and dust-to-metals ratio for reference abundances.
    my $i;
    my @elements                 = sort {$abundancesReference{$a}->{'atomicNumber'} <=> $abundancesReference{$b}->{'atomicNumber'}} keys(%abundancesReference);
    $i = -1;
    my $isMetal                  = pdl map {++$i;$abundancesReference{$_}->{'atomicNumber'} >  2 ? $i : ()                                                                                               } @elements;
    $i = -1;
    my $isNotMetal               = pdl map {++$i;$abundancesReference{$_}->{'atomicNumber'} <= 2 ? $i : ()                                                                                               } @elements;
    my $abundancesByMass         = pdl map {     $abundancesReference{$_}->{'atomicMass'}*10.0**$abundancesReference{$_}->{'logAbundanceByNumber'}                                                       } @elements;
    my $dustByMass               = pdl map {     $abundancesReference{$_}->{'atomicMass'}*10.0**$abundancesReference{$_}->{'logAbundanceByNumber'}*(1.0-$abundancesReference{$_}->{'undepletedFraction'})} @elements;
    my $metallicityReference     = +$abundancesByMass->($isMetal)->sum()
	                           /$abundancesByMass            ->sum();
    my $dustMetalsRatioReference = +$dustByMass      ->($isMetal)->sum()
                                   /$abundancesByMass->($isMetal)->sum();
    # Adjust abundance of all metals with the overall metallicity.
    foreach my $element ( @elements ) {
	next
	    unless ( $abundances{$element}->{'atomicNumber'} > 2 );
	$abundances{$element}->{'logAbundanceByNumber'} += $metallicity->log10();
    }
    # Apply any custom adjustments to individual elements.
    foreach my $element ( @elements ) {
	&{$abundances{$element}->{'adjustAbundance'}}(\%abundances,$metallicity*$metallicityReference)
	    if ( exists($abundances{$element}->{'adjustAbundance'}) );
    }
    # Renormalize to keep the total metallicity fixed.
    my $abundancesByMassNew   = pdl map {$abundances{$_}->{'atomicMass'}*10.0**$abundances{$_}->{'logAbundanceByNumber'}} @elements;
    my $renormalizationFactor =      +$metallicity*$metallicityReference
	                        /(1.0-$metallicity*$metallicityReference)
                                *$abundancesByMassNew->($isNotMetal)->sum()
				/$abundancesByMassNew->($isMetal   )->sum();
    foreach my $element ( @elements ) {
	# Skip non-metals.
	next
	    unless ( $abundances{$element}->{'atomicNumber'} > 2 );
	# Renormalize.
	$abundances{$element}->{'logAbundanceByNumber'} += $renormalizationFactor->log10();
    }
    # Adjust depletion factors. Here we linearly interpolate between the reference depletion fractions, f☉, and dust-to-metals
    # ratio ξ☉, and the limiting cases of (f,ξ)=(0,0) and (f,ξ)=(1,1) as per Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; §2.3.2).
    foreach my $element ( @elements ) {
	my $depletionFractionReference = 1.0-$abundancesReference{$element}->{'undepletedFraction'};
	my $depletionFraction;
	if ( $dustToMetalsRatio < $dustMetalsRatioReference ) {
	    $depletionFraction = $depletionFractionReference*$dustToMetalsRatio/$dustMetalsRatioReference;
	} else {
	    $depletionFraction = $depletionFractionReference+($dustToMetalsRatio-$dustMetalsRatioReference)/(1.0-$dustMetalsRatioReference)*(1.0-$depletionFractionReference);
	}
	my $undepletedFraction = 1.0-$depletionFraction;
	$abundances{$element}->{'undepletedFraction'} = $undepletedFraction;
    }
    # Return the modified set of abundances.
    return %abundances;
}

sub adjustAbundanceNitrogen {
    # Adjust the abundance of nitrogen following the model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; eqn. 11).
    my %abundances                             = %{shift()};
    my $metallicity                            =   shift() ;
    my $logOH                                  = $abundances{'O'}->{'logAbundanceByNumber'};
    my $OH                                     = 10.0**$logOH;
    my $NH                                     = 0.41*$OH*(10.0**(-1.6)+10.0**(2.33+$logOH));
    my $logNH                                  = log10($NH);
    $abundances{'N'}->{'logAbundanceByNumber'} = $logNH;
}

sub adjustAbundanceCarbon {
    # Adjust the abundance of nitrogen following the model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G).
    my %abundances                             = %{shift()};
    my $metallicity                            =   shift() ;
    my $logOH                                  = $abundances{'O'}->{'logAbundanceByNumber'};
    $abundances{'C'}->{'logAbundanceByNumber'} = log10(0.44)+$logOH;
}

sub adjustAbundanceHelium {
    # Adjust the abundance of nitrogen following the model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; eqn. 12).
    my %abundances                             = %{shift()};
    my $metallicity                            =   shift() ;
    my $massFractionHe                         = 0.2485+1.7756*$metallicity;
    my $massFractionH                          = 1.0-$massFractionHe-$metallicity;
    my $HeH                                    = +($massFractionHe/$abundances{'He'}->{'atomicMass'})
	                                         /($massFractionH /$abundances{'H' }->{'atomicMass'});
    my $logHeH                                 = log10($HeH);
    $abundances{'He'}->{'logAbundanceByNumber'} = $logHeH;
}

sub validateResults {
    # Validate the results of the Cloudy calculations.
    my $ages                      = shift();
    my $logMetallicities          = shift();
    my $logHydrogenLuminosities   = shift();
    my $logHydrogenDensities      = shift();
    my $ionizingLuminosityPerMass = shift();
    my $lineData                  = shift();
    
    # Validate the ratio L_Hα/Qₕ. From Osterbrock & Ferland (2006; https://ui.adsabs.harvard.edu/abs/2006agna.book.....O) we
    # expect:
    #
    #    (L_Hα/ergs s¯¹) / (Qₕ/photons s¯¹) = (α^eff_Hα/α_B) hν_Hα = 1.37 x 10¯¹²
    #    
    # where α^eff_Hα is the effective recombination coefficient for Hα, α_B is the case B recombination coefficient, and hν_Hα is
    # the energy of an Hα photon. The numerical value above is for a pure hydrogen/helium nebula, with temperature of 10,000K and
    # density n_e = 100 cm¯³.
    my $selectAge         = which($ages                                                    <  0.01                    ); # Consider young (< 10 Myr) populations.
    my $selectMetallicity = which($logMetallicities                                        == $logMetallicities->((0))); # Select primordial metallicity.
    my $selection         = which($lineData->{'status'}->($selectAge,$selectMetallicity,:) == 0                       ); # Select successful models.
    # Compute the ratio.
    my $ratio = $lineData->{'balmerAlpha6565'}->{'luminosity'}->copy();
    for(my $i=0;$i<$ratio->dim(2);++$i) {
	$ratio->(:,:,($i)) /= 10.0**$logHydrogenLuminosities->(($i));
    }
    # Find the minimum and maximum ratios.
    my $ratioMinimum = $ratio->($selectAge,$selectMetallicity,:)->flat()->($selection)->minimum();
    my $ratioMaximum = $ratio->($selectAge,$selectMetallicity,:)->flat()->($selection)->maximum();
    # Report.
    my $ratioTarget = pdl 1.37e-12;
    print "Validation: ratio L_Hα/Qₕ\n";
    print "   Expected: ".$ratioTarget."\n";
    print "   Found range: ".sprintf("%8.3e",$ratioMinimum)." to ".sprintf("%8.3e",$ratioMaximum)."\n\n";

    # Validate the ionizing photon luminosity per star formation rate. Compute the ratio of star formation rate to ionizing
    # luminosity. This assumes a constant star formation rate of φ=1 M☉/yr. The ratio is then:
    #
    #    (SFR/M☉ yr¯¹)/(Qₕ/photons s¯¹) = φ / ∫₀᪲ φ Qₕ(t) dt
    #
    # a value of 7.4 x 10¯⁵⁴ is expected for a Kroupa IMF (Osterbrock & Ferland, 2006;
    # https://ui.adsabs.harvard.edu/abs/2006agna.book.....O). For a Chabrier IMF a lower value of around 4.7 x 10¯⁵⁴ is expected.
    my $ageStep       =  $ages   ->copy();
    $ageStep->(0:-2) .= +$ages   ->(1:-1)
	                -$ages   ->(0:-2);
    $ageStep->((-1)) .=  $ageStep->((-2));
    my $giga          = pdl 1.0e9;
    my $starFormtionRateToIonizingLuminosity = 1.0/sum($ionizingLuminosityPerMass->(:,(0))*$ageStep*$giga);
    # Report.
    print "Validation: ratio SFR/Qₕ\n";
    print "   Expected (Kroupa IMF): 7.54e-54\n";
    print "   Expected (Chabrier IMF): 4.74e-54\n";
    print "   Found: ".$starFormtionRateToIonizingLuminosity."\n";
    
    # Validate model success.
    my $selectSuccess       = which($lineData->{'status'} == 0);
    my $selectDisaster      = which($lineData->{'status'} == 1);
    my $selectExit          = which($lineData->{'status'} == 2);
    my $selectMissingOutput = which($lineData->{'status'} == 3);
    my $selectMissingLine   = which($lineData->{'status'} == 4);
    my $modelsTotal         = nelem($lineData->{'status'}     );
    my $modelsSuccess       = nelem($selectSuccess           );
    my $modelsDisaster      = nelem($selectDisaster          );
    my $modelsExit          = nelem($selectExit              );
    my $modelsMissingOutput = nelem($selectMissingOutput     );
    my $modelsMissingLine   = nelem($selectMissingLine       );
    print "Validation: model success\n";
    print "      Total models: ".$modelsTotal        ."\n";
    print "        Successful: ".$modelsSuccess      ."\n";
    print "         Disasters: ".$modelsDisaster     ."\n";
    print "    Bad exit codes: ".$modelsExit         ."\n";
    print "   Missing outputs: ".$modelsMissingOutput."\n";
    print "     Missing lines: ".$modelsMissingLine  ."\n";
}
