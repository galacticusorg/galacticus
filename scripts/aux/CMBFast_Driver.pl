#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;
use File::Copy;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");

# Driver script for CMBFast.
# Andrew Benson (27-Nov-2009)

# Get arguments.
if ( $#ARGV != 3 ) {die "Usage: CMBFast_Driver.pl <parameterFile> <transferFunctionFile> <kMax> <fileFormatVersion>"};
$parameterFile        = $ARGV[0];
$transferFunctionFile = $ARGV[1];
$kMax                 = $ARGV[2];
$fileFormat           = $ARGV[3];

# Ensure the requested file format version is compatible.
my $fileFormatCurrent = 1;
die('CMBFast_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Download the code.
unless ( -e $galacticusPath."aux/cmbfast.tar.gz" ) {
    print "CMBFast_Driver.pl: downloading CMBFast code.\n";
    system("wget http://web.archive.org/web/20070205144812/http://cfa-www.harvard.edu/~mzaldarr/CMBFAST/cmbfast.tar.gz -O ".$galacticusPath."aux/cmbfast.tar.gz");
    die("CMBFast_Driver.pl: FATAL - failed to download CMBFast code.") unless ( -e $galacticusPath."aux/cmbfast.tar.gz" );
}

# Unpack the code.
unless ( -e $galacticusPath."aux/cmbfast4.5.1" ) {
    print "CMBFast_Driver.pl: unpacking CMBFast code.\n";
    system("tar -x -v -z -C ".$galacticusPath."aux -f ".$galacticusPath."aux/cmbfast.tar.gz");
    die("CMBFast_Driver.pl: FATAL - failed to unpack CMBFast code.") unless ( -e $galacticusPath."aux/cmbfast4.5.1" );
}

# Patch the code.
unless ( -e $galacticusPath."aux/cmbfast4.5.1/configure.patch" ) {
    foreach $file ( "cmbflat.F.patch", "configure.patch", "recfast.f.patch" ) {
	copy($galacticusPath."aux/cmbfast4.5.1_Galacticus_Modifications/".$file,$galacticusPath."aux/cmbfast4.5.1/".$file);
	if ( $file =~ m/\.patch$/ ) {system("cd ".$galacticusPath."aux/cmbfast4.5.1; patch < $file")};
	print "$file\n";
    }
}

# Build the code.
unless ( -e $galacticusPath."aux/cmbfast4.5.1/cmb" ) {
    print "CMBFast_Driver.pl: compiling CMBFast code.\n";
    system("cd ".$galacticusPath."aux/cmbfast4.5.1/; configure --f77=gfortran --f77-flags=-ffixed-line-length-132 --f77-flags=-O4 --with-hispdhikmodes=yes; make");
    die("CMBFast_Driver.pl: FATAL - failed to build CMBFast code.") unless ( -e $galacticusPath."aux/cmbfast4.5.1/cmb" );
}

# Specify default parameters.
$Omega_nu          = 0.0;
$NNeutrino         = 0;
$NNeutrinoMassless = 3.04;
$gStarMassive      = 10.75;
$NPerLogk          = 10;

# Parse the parameter file.
$xml = new XML::Simple;
$data = $xml->XMLin($parameterFile);
$parameterHash = $data->{'parameter'};

# Check that required parameters exist.
@parameters = ( "Omega_b", "Omega_Matter", "Omega_DE", "H_0", "T_CMB", "Y_He" );
foreach $parameter ( @parameters ) {
    die("CMBFast_Driver.pl: FATAL - parameter ".$parameter." can not be found.") unless ( exists($data->{'parameter'}->{$parameter}) );
}

# Calculate derived parameters.
$Omega_c = $parameterHash->{'Omega_Matter'}->{'value'}-$parameterHash->{'Omega_b'}->{'value'};
$kMax = $kMax/($parameterHash->{'H_0'}->{'value'}/100.0);

my $makeFile = 0;
if ( -e $transferFunctionFile ) {
    my $xmlDoc = new XML::Simple;
    $transferFunction = $xmlDoc->XMLin($transferFunctionFile);
    if ( exists($transferFunction->{'fileFormat'}) ) { 
	$makeFile = 1 unless ( $transferFunction->{'fileFormat'} == $fileFormatCurrent );
    } else {
	$makeFile = 1;
    }
} else {
    $makeFile = 1;
}

# Create the file if necessary.
if ( $makeFile == 1 ) {
   # Run CMBFast.
   open(cmbPipe,"|".$galacticusPath."aux/cmbfast4.5.1/cmb");
   print cmbPipe "1\n";
   print cmbPipe "$kMax $NPerLogk\n";
   print cmbPipe "1 0\n";
   print cmbPipe $galacticusPath."data/transfer_function.tmp\n";
   print cmbPipe "1\n";
   print cmbPipe "-1\n";
   print cmbPipe $parameterHash->{'Omega_b'}->{'value'}." ";
   print cmbPipe $Omega_c." ";
   print cmbPipe $parameterHash->{'Omega_DE'}->{'value'}." ";
   print cmbPipe "$Omega_nu\n";
   print cmbPipe $parameterHash->{'H_0'}->{'value'}." ";
   print cmbPipe $parameterHash->{'T_CMB'}->{'value'}." ";
   print cmbPipe $parameterHash->{'Y_He'}->{'value'}." ";
   print cmbPipe "$NNeutrinoMassless $NNeutrino $gStarMassive\n";
   print cmbPipe "1\n";
   print cmbPipe "0\n";
   print cmbPipe "1\n";
   close(cmbPipe);
   
   # Read in the tabulated data and output as an XML file.
   open(inHndl,$galacticusPath."data/transfer_function.tmp");
   while ( $line = <inHndl> ) {
       $line =~ s/^\s*//;
       @columns = split(/\s+/,$line);
       $k = $columns[0]*($parameterHash->{'H_0'}->{'value'}/100.0);
       $transferFunctionData[++$#transferFunctionData] = $k." ".$columns[1];
   }
   close(inHndl);
   unlink($galacticusPath."data/transfer_function.tmp","nu1.dat",$parameterFile);
   @{$transferFunction{'datum'}} = @transferFunctionData;
   @{$transferFunction{'column'}} = (
       "k [Mpc^{-1}] - wavenumber",
       "T(k) - transfer function"
       );
   @{$transferFunction{'description'}} = "Cold dark matter power spectrum created by CMBFast.";
   foreach $parameter ( @parameters ) {
       %{${$transferFunction{'parameter'}}[++$#{$transferFunction{'parameter'}}]} = (
	   "name" => $parameter,
	   "value" => $parameterHash->{$parameter}->{'value'}
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
   $transferFunction = \%transferFunction;
   $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"data");
   open(outHndl,">".$transferFunctionFile);
   print outHndl $xmlOutput->XMLout($transferFunction);
   close(outHndl);
}

exit;
