# Contains a Perl module which implements calculation of cross-sections in the ISM.

package ISMCrossSections;
use strict;
use warnings;
use utf8;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use XML::Simple;

# Variables to hold fitting function data.
our $mm83Energy;
our $mm83C0;
our $mm83C1;
our $mm83C2;
our $w00Energy;
our $w00Metallicity;
our $w00CrossSection;

sub Cross_Sections {
    # Return a PDL of cross sections at the specified energies. Allows optional specification of the metallicity and model to use.
    my $energy = shift;
    # Extract any remaining options.
    my (%options) = @_
	if ( scalar(@_) >= 1 );

    # Select model and metallicity.
    my $model       = "Morrison-McCammon1983";
    $model          = $options{'model'      }
       if ( exists($options{'model'      }) );
    my $metallicity = pdl 0.012000314;
    $metallicity    = $options{'metallicity'}
       if ( exists($options{'metallicity'}) );

    # Define constants.
    my $thomsonCrossSection = pdl 6.65245854533e-25; # cm^2
    
    # Initialize the cross-section PDL.
    my $crossSection;

    # Branch on model selection.
    if ( $model eq "Morrison-McCammon1983" ) {
	unless ( defined($mm83Energy) ) {
	    my $xml           = new XML::Simple;
	    my $absorptionFit = $xml->XMLin("data/atomic/Interstellar_Photoelectric_Absorption_Morrison_McCammon_1983.xml");
	    $mm83Energy       = pdl @{$absorptionFit->{'energy'}->{'datum'}};
	    $mm83C0           = pdl @{$absorptionFit->{'c0'    }->{'datum'}};
	    $mm83C1           = pdl @{$absorptionFit->{'c1'    }->{'datum'}};
	    $mm83C2           = pdl @{$absorptionFit->{'c2'    }->{'datum'}};
	}    
	my $energyRanges  = vsearch($energy,$mm83Energy);
	# Photo-electric absorption cross-section (cm^2).
	$crossSection    = 
	    (
	     ($mm83C0->index($energyRanges))           + 
	     ($mm83C1->index($energyRanges))*$energy   + 
	     ($mm83C2->index($energyRanges))*$energy**2
	    )*($energy**(-3))*1.0e-24;
	
	# Add Thompson cross-section assuming 1.2 electrons per hydrogen atom.
	$crossSection += 1.2*$thomsonCrossSection;
    } elsif ( $model eq "Wilms2000" ) {
	# Specify the files to download.
	my @files = (
	    {
		name => "dotbvabs.f",
		url  => "ftp://heasarc.gsfc.nasa.gov/software/lheasoft/release/current/Xspec/src/XSFunctions/"
	    },
	    {
		name => "gphoto.f",
		url  => "ftp://heasarc.gsfc.nasa.gov/software/lheasoft/release/current/Xspec/src/XSFunctions/"
	    },
	    {
		name => "phfit2.f",
		url  => "ftp://heasarc.gsfc.nasa.gov/software/lheasoft/release/current/Xspec/src/XSFunctions/"
	    },
	    {
		name => "j4save.f",
		url  => "http://www.netlib.org/alliant/quad/"
	    }
	    );
	# Download the files that we need to compute the absorption.
	foreach my $file ( @files ) {
	    system("mkdir -p aux/XSpec; wget ".$file->{'url'}.$file->{'name'}." -O aux/XSpec/".$file->{'name'})
		unless ( -e "aux/XSpec/".$file->{'name'} );
	}
	# Build the wrapper code that will be used to generate the table of absorptions.
	system("make XRay_Absorption_ISM_Wilms2000.exe")
	    unless ( -e "XRay_Absorption_ISM_Wilms2000.exe" );
	die("Unable to build XRay_Absorption_ISM_Wilms2000.exe")
	    unless ( -e "XRay_Absorption_ISM_Wilms2000.exe" );
	# Build the absorption table.
	system("XRay_Absorption_ISM_Wilms2000.exe")
	    unless ( -e "data/atomic/Interstellar_Absorption_Wilms_2000.hdf5" );
	# Read the data file.
	unless ( defined($w00Energy) ) {
	    my $absorptionFile = new PDL::IO::HDF5("data/atomic/Interstellar_Absorption_Wilms_2000.hdf5");
	    $w00Energy       = $absorptionFile->dataset('energy'      )->get();
	    $w00Metallicity  = $absorptionFile->dataset('metallicity' )->get();
	    $w00CrossSection = $absorptionFile->dataset('crossSection')->get();
	}
	# Interpolate in metallicity and energy to get cross-sections.
	(my $energyIndices   ,my $energyError     ) = interpolate($energy     ,$w00Energy     ,sequence(nelem($w00Energy     )));
	(my $metallicityIndex,my $metallicityError) = interpolate($metallicity,$w00Metallicity,sequence(nelem($w00Metallicity)));
	my $interpIndices = pdl zeroes(2,nelem($energy));
	$interpIndices->((0),:) .= $energyIndices;
	$interpIndices->((1),:) .= $metallicityIndex;
	$crossSection = $w00CrossSection->interpND($interpIndices);
    } else {
	die('ISMCrossSections::Cross_Sections(): unknown model')
    }

    # Return the computed cross-sections.
    return $crossSection;
}

1;
