#!/usr/bin/perl
use strict;
use Astro::FITS::CFITSIO qw( :constants );
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::IO::HDF5;
use Carp;
Astro::FITS::CFITSIO::PerlyUnpacking(0);

# Fit polynomials to the mass and redshift dependence of stellar mass errors in the PRIMUS survey fields.
# Andrew Benson (22-May-2014)

# Get field solid angles.
my $solidAngleFile = new PDL::IO::HDF5("constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/solidAngles.hdf5");
my $solidAngles    = $solidAngleFile->dataset("solidAngle")->get();
my $weights        = $solidAngles/sum($solidAngles);
# Define files for each field.
my @fileNames = 
    (
     "cosmos_benson.fits"    ,  
     "xmm_swire_benson.fits" ,
     "cfhtls_xmm_benson.fits",  
     "cdfs_benson.fits"      ,  
     "es1_benson.fits"
    );
# Iterate over files.
my $iField = 0;
open(my $fortranHndl,">constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massErrors.F90");
open(my $latexHndl  ,">constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massErrors.tex");
for my $fileName ( @fileNames ) {
    ++$iField;
    # Open the file and select the relevant HDU.
    my $status;
    my $filePath = "constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/".$fileName;
    die("massErrors.pl: this script requires the PRIMUS data tables be available")
	unless ( -e $filePath );
    my $file = Astro::FITS::CFITSIO::open_file($filePath,READONLY,$status);
    $file->movabs_hdu(2,ANY_HDU,$status);
    # Extract the number of rows in the table.
    $file->get_num_rows(my $rowCount,$status);
    # Extract the data from the required columns.
    my $cols;
    foreach ( "ZPRIMUS", "MASS_50", "MASS_ERR" ) {
	$file->get_colnum(0,$_,my $dataColumn,$status);
	$cols->{$_} = zeroes($rowCount)->double;
	$file->read_col(Astro::FITS::CFITSIO::TDOUBLE(),$dataColumn,1,1,$rowCount,0,${$cols->{$_}->get_dataref},undef,$status);
	$cols->{$_}->upd_data;
    }
    # Close the file.
    $file->close_file($status);
    # Construct matrices required to optimize the polynomial coefficients.
    my @sequence = # Sequence of exponents on mass and redshift.
	(
	 [0,0],
	 [1,0],
	 [2,0],
	 [0,1],
	 [0,2],
	 [1,1] 
	);
    my $b = pdl zeroes(scalar(@sequence)                  );
    my $a = pdl zeroes(scalar(@sequence),scalar(@sequence));
    for(my $i=0;$i<scalar(@sequence);++$i) {
	my $miExponent = ${$sequence[$i]}[0];
	my $ziExponent = ${$sequence[$i]}[1];
	$b->(($i)) .= sum($cols->{'MASS_ERR'}*$cols->{'MASS_50'}**$miExponent*$cols->{'ZPRIMUS'}**$ziExponent);
	for(my $j=0;$j<scalar(@sequence);++$j) {
	    my $mjExponent = ${$sequence[$j]}[0];
	    my $zjExponent = ${$sequence[$j]}[1];
	    $a->(($i),($j)) .= 
		sum(
		    $cols->{'MASS_50'}**$miExponent*$cols->{'ZPRIMUS'}**$ziExponent
		    *
		    $cols->{'MASS_50'}**$mjExponent*$cols->{'ZPRIMUS'}**$zjExponent
		);
	}
    }
    # Solve for polynomial coefficients.
    my $x = msolve($a,transpose($b));
    # Output results.
    (my $fieldName = uc($fileName)) =~ s/_BENSON\.FITS//;
    $fieldName =~ s/_/ /g;
    print $fortranHndl "   ! Error fit for field: ".$fieldName."\n";
    print $fortranHndl "   error(".$iField.")=                                             &\n";
    for(my $i=0;$i<scalar(@sequence);++$i) {
	my $miExponent   = ${$sequence[$i]}[0];
	my $ziExponent   = ${$sequence[$i]}[1];
	my $coefficient  = $x->((0),($i));
	(my $pCoefficient = sprintf("%16.10e",$coefficient)) =~ tr/e/d/;
	my $prefix      = "";
	$prefix         = "+"
	    if ( $coefficient >= 0.0 );
	my $massTerm = "*logMass**" .$miExponent;
	$massTerm    = "*logMass   "
	    if ( $miExponent == 1 );
	$massTerm    = "           "
	    if ( $miExponent == 0 );
	my $redshiftTerm = "*redshift**" .$ziExponent;
	$redshiftTerm    = "*redshift   "
	    if ( $ziExponent == 1 );
	$redshiftTerm    = "            "
	    if ( $ziExponent == 0 );
	print $fortranHndl "      &         ".$prefix.$pCoefficient.$massTerm.$redshiftTerm;
	print $fortranHndl " & "
	    unless ( $i == scalar(@sequence)-1 );
	print $fortranHndl "\n";
    }

    print $latexHndl $fieldName;
    for(my $i=0;$i<scalar(@sequence);++$i) {
	my $coefficient  = $x->((0),($i));
	my $pCoefficient = sprintf("%6.4f",$coefficient);
	my $prefix      = "";
	$prefix         = "+"
	    if ( $coefficient >= 0.0 );
	print $latexHndl " & ".$prefix.$pCoefficient;
    }
    print $latexHndl " \\\\\n";
}
close($fortranHndl);
close($latexHndl  );

exit;
