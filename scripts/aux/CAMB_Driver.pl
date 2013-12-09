#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V093'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V093'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use XML::Simple;
use Data::Dumper;
use File::Copy;
use IO::Prompt;
use IO::Interactive qw( is_interactive );

# Driver script for CAMB.
# Andrew Benson (2-December-2013)

# Get arguments.
die "Usage: CAMB_Driver.pl <parameterFile> <transferFunctionFile> <kMax> <fileFormatVersion>"
    unless ( scalar(@ARGV) == 4 );
my $parameterFile        = $ARGV[0];
my $transferFunctionFile = $ARGV[1];
my $kMax                 = $ARGV[2];
my $fileFormat           = $ARGV[3];

# Ensure the requested file format version is compatible.
my $fileFormatCurrent = 1;
die('CAMB_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Download the code.
unless ( -e $galacticusPath."aux/CAMB.tar.gz" ) {
    if (is_interactive()) {
	print "Have you read and accepted the CAMB license at http://camb.info/CAMBsubmit.html ? (y/n)\n";
	my $key = prompt('',-ty1);
	die("You must read and accept the CAMB license before it can be used.")
	    unless ( lc($key) eq "y" );
    } else {
	print "WARNING: About to download CAMB, note that you should have read accepted the CAMB license at http://camb.info/CAMBsubmit.html\n";
    }
    system("wget http://camb.info/CAMB.tar.gz -O ".$galacticusPath."aux/CAMB.tar.gz");
    die("CAMB_Driver.pl: unable to download CAMB")
	unless ( -e $galacticusPath."aux/CAMB.tar.gz" );
}

# Unpack the code.
unless ( -e $galacticusPath."aux/camb" ) {
    print "CAMB_Driver.pl: unpacking CAMB code.\n";
    system("tar -x -v -z -C ".$galacticusPath."aux -f ".$galacticusPath."aux/CAMB.tar.gz");
    die("CAMB_Driver.pl: FATAL - failed to unpack CAMB code.")
	unless ( -e $galacticusPath."aux/camb" );
}

# Build the code.
unless ( -e $galacticusPath."aux/camb/camb" ) {
    print "CAMB_Driver.pl: compiling CMBFast code.\n";
    system("cd ".$galacticusPath."aux/camb/; sed -r -i~ s/\"F90C\\s*=\\s*.*\"/\"F90C = gfortran\"/ Makefile; sed -r -i~ s/\"^FFLAGS\\s*=\\s*.*\"/\"FFLAGS = -Ofast -march=native -fopenmp\"/ Makefile; sed -r -i~ /\"F90CRLINK\\s*=\\s*.*\"/d Makefile; make -j1");
    die("CAMB_Driver.pl: FATAL - failed to build CAMB code.")
	unless ( -e $galacticusPath."aux/camb/camb" );
}

# Specify default parameters.
my $Omega_nu          = 0.0;

# Parse the parameter file.
my $xml           = new XML::Simple;
my $data          = $xml->XMLin($parameterFile);
my $parameterHash = $data->{'parameter'};

# Check that required parameters exist.
my @parameters = ( "Omega_b", "Omega_Matter", "Omega_DE", "H_0", "T_CMB", "Y_He" );
foreach my $parameter ( @parameters ) {
    die("CAMB_Driver.pl: FATAL - parameter ".$parameter." can not be found.")
	unless ( exists($data->{'parameter'}->{$parameter}) );
}

# Calculate derived parameters.
my $Omega_c = $parameterHash->{'Omega_Matter'}->{'value'}-$parameterHash->{'Omega_b'}->{'value'};
$kMax = $kMax/($parameterHash->{'H_0'}->{'value'}/100.0);

my $makeFile = 0;
if ( -e $transferFunctionFile ) {
    my $xmlDoc           = new XML::Simple;
    my $transferFunction = $xmlDoc->XMLin($transferFunctionFile);
    if ( exists($transferFunction->{'fileFormat'}) ) { 
	$makeFile = 1
	    unless ( $transferFunction->{'fileFormat'} == $fileFormatCurrent );
    } else {
	$makeFile = 1;
    }
} else {
    $makeFile = 1;
}

# Create the file if necessary.
if ( $makeFile == 1 ) {
   # Create the directory.
   system("mkdir -p `dirname ".$transferFunctionFile."`");
   # Run CAMB.
   open(cambInput,">".$transferFunctionFile.".inp");
   print cambInput <<END;
output_root = camb
get_scalar_cls = F
get_vector_cls = F
get_tensor_cls = F
get_transfer   = T
do_lensing     = F
do_nonlinear = 0
l_max_scalar      = 2200
l_max_tensor      = 1500
k_eta_max_tensor  = 3000
use_physical   = F
omega_baryon   = $parameterHash->{'Omega_b'}->{'value'}
omega_cdm      = $Omega_c
omega_lambda   = $parameterHash->{'Omega_DE'}->{'value'}
omega_neutrino = $Omega_nu
omk            = 0
hubble         = $parameterHash->{'H_0'}->{'value'}
w              = -1
cs2_lam        = 1
temp_cmb           = $parameterHash->{'T_CMB'}->{'value'}
helium_fraction    = $parameterHash->{'Y_He'}->{'value'}
massless_neutrinos = 2.046
nu_mass_eigenstates = 1
massive_neutrinos  = 1
share_delta_neff = T
nu_mass_fractions = 1
nu_mass_degeneracies = 
initial_power_num         = 1
pivot_scalar              = 0.05
pivot_tensor              = 0.05
scalar_amp(1)             = 2.1e-9
scalar_spectral_index(1)  = 0.96
scalar_nrun(1)            = 0
tensor_spectral_index(1)  = 0
initial_ratio(1)          = 1
reionization         = T
re_use_optical_depth = T
re_optical_depth     = 0.09
re_redshift          = 11
re_delta_redshift    = 1.5
re_ionization_frac   = -1
RECFAST_fudge = 1.14
RECFAST_fudge_He = 0.86
RECFAST_Heswitch = 6
RECFAST_Hswitch  = T
initial_condition   = 1
initial_vector = -1 0 0 0 0
vector_mode = 0
COBE_normalize = F
CMB_outputscale = 7.42835025e12 
transfer_high_precision = F
transfer_kmax           = $kMax
transfer_k_per_logint   = 0
transfer_num_redshifts  = 1
transfer_interp_matterpower = T
transfer_redshift(1)    = 0
transfer_filename(1)    = transfer_out.dat
transfer_matterpower(1) = matterpower.dat
scalar_output_file = scalCls.dat
vector_output_file = vecCls.dat
tensor_output_file = tensCls.dat
total_output_file  = totCls.dat
lensed_output_file = lensedCls.dat
lensed_total_output_file  =lensedtotCls.dat
lens_potential_output_file = lenspotentialCls.dat
FITS_filename      = scalCls.fits
do_lensing_bispectrum = F
do_primordial_bispectrum = F
bispectrum_nfields = 1
bispectrum_slice_base_L = 0
bispectrum_ndelta=3
bispectrum_delta(1)=0
bispectrum_delta(2)=2
bispectrum_delta(3)=4
bispectrum_do_fisher= F
bispectrum_fisher_noise=0
bispectrum_fisher_noise_pol=0
bispectrum_fisher_fwhm_arcmin=7
bispectrum_full_output_file=
bispectrum_full_output_sparse=F
bispectrum_export_alpha_beta=F
feedback_level = 1
derived_parameters = T
lensing_method = 1
accurate_BB = F
massive_nu_approx = 1
accurate_polarization   = T
accurate_reionization   = T
do_tensor_neutrinos     = T
do_late_rad_truncation   = T
number_of_threads       = 0
high_accuracy_default=T
accuracy_boost          = 1
l_accuracy_boost        = 1
l_sample_boost          = 1
END
   close(cambInput);

   # Run CAMB.
   system($galacticusPath."aux/camb/camb ".$transferFunctionFile.".inp");

   # Read in the tabulated data and output as an XML file.
   my @transferFunctionData;
   my %transferFunction;
   open(inHndl,"camb_transfer_out.dat");
   while ( my $line = <inHndl> ) {
       $line =~ s/^\s*//;
       my @columns = split(/\s+/,$line);
       my $k = $columns[0]*($parameterHash->{'H_0'}->{'value'}/100.0);
       push(@transferFunctionData,$k." ".$columns[1]);
   }
   close(inHndl);
   unlink("camb_params.ini","camb_transfer_out.dat","camb_matterpower.dat",$transferFunctionFile.".inp");
   @{$transferFunction{'datum'}} = @transferFunctionData;
   @{$transferFunction{'column'}} = (
       "k [Mpc^{-1}] - wavenumber",
       "T(k) - transfer function"
       );
   @{$transferFunction{'description'}} = "Cold dark matter power spectrum created by CAMB.";
   foreach my $parameter ( keys(%{$parameterHash}) ) {
       push(@{$transferFunction{'parameter'}},	    
	    {
		"name"  => $parameter,
		"value" => $parameterHash->{$parameter}->{'value'}
	    }
	   );
   }
   # Add extrapolation data.
   ${$transferFunction{'extrapolation'}->{'wavenumber'}}[0]->{'limit' } = "low";
   ${$transferFunction{'extrapolation'}->{'wavenumber'}}[0]->{'method'} = "power law";
   ${$transferFunction{'extrapolation'}->{'wavenumber'}}[1]->{'limit' } = "high";
   ${$transferFunction{'extrapolation'}->{'wavenumber'}}[1]->{'method'} = "power law";
   # Add file format version.
   $transferFunction{'fileFormat'} = $fileFormatCurrent;
   # Add unique label.
   $transferFunction{'uniqueLabel'} = $data->{'uniqueLabel'};
   # Output the transfer function.
   my $transferFunctionStructure = \%transferFunction;
   my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"data");
   open(outHndl,">".$transferFunctionFile);
   print outHndl $xmlOutput->XMLout($transferFunctionStructure);
   close(outHndl);
}

exit;
