#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
require Mangle;

# Generate C_l's from a survey window using mangle (http://space.mit.edu/~molly/mangle/).
# Andrew Benson (18-April-2014)

# Get arguments.
die("Usage: mangleHarmonize.pl <file1> [<file2>]....... <outputFile> <lmax>")
    unless ( scalar(@ARGV) > 2 );
my @files          = @ARGV;
my $lMax           = pop(@files);
my $outputFileName = pop(@files);

# Build mangle.
&Mangle::Build;

# Count number of fields.
my $fieldCount = scalar(@files);

# Compute and sum spherical harmonic coefficients for each file.
my @Wlms;
foreach my $file ( @files ) {
    system("time ".$galacticusPath."../mangle/bin/harmonize -l".$lMax." ".$file." ".$file.".wlm")
	unless ( -e $file.".wlm" );
    die("mangleHarmonize.pl: failed to generate wlm file")
	unless ( -e $file.".wlm" );
    open(my $wlmFile,$file.".wlm");
    my $header = <$wlmFile>;
    my $Wlm    = pdl zeroes(2,($lMax+1)*($lMax+2)/2);
    my $i      = -1;
    for(my $l=0;$l<=$lMax;++$l) {
	for (my $m=0;$m<=$l;++$m) {
	    ++$i;
	    (my $line = <$wlmFile>) =~ s/^\s*(.*?)\s*$/$1/;
	    $Wlm->(:,($i)) .= pdl split(/\s+/,$line);
	}
    }
    close($wlmFile);
    push(@Wlms,$Wlm);
}

# Compute C_l's.
my $Cls  = pdl zeroes(                            $lMax+1);
my $Clfs = pdl zeroes($fieldCount*($fieldCount+1),$lMax+1);
my $i    = -1;
for(my $l=0;$l<=$lMax;++$l) {
    my $Cl  = pdl 0.0;
    my $Clf = pdl zeroes($fieldCount*($fieldCount+1));
    for (my $m=0;$m<=$l;++$m) {
	++$i;
	my $weight = pdl 2.0;
	$weight .= 1.0
	    if ( $m == 0 );
	my $Wlm = pdl zeroes(2);
	foreach ( @Wlms ) {
	    $Wlm += $_->(:,($i));
	}
	$Cl += $weight*sum($Wlm**2);
	my $j = -1;
	for(my $p=0;$p<$fieldCount;++$p) {
	    for(my $q=$p;$q<$fieldCount;++$q) {
		++$j;
		$Clf->(($j)) += $weight*sum($Wlms[$p]->(:,($i))*$Wlms[$q]->(:,($i)));
	    }
	}
    }
    $Cl             /= (2.0*$l+1.0);
    $Clf            /= (2.0*$l+1.0);
    $Cls ->(  ($l)) .= $Cl ;
    $Clfs->(:,($l)) .= $Clf;
}

# Store results to file.
my $outputFile = new PDL::IO::HDF5(">".$outputFileName);
my $ls         = pdl sequence(1+$lMax);
$outputFile->attrSet(files => join(" : ",@files));
$outputFile->dataset('l' )->set($ls );
$outputFile->dataset('Cl')->set($Cls);
my $j = -1;
for(my $p=0;$p<$fieldCount;++$p) {
    for(my $q=$p;$q<$fieldCount;++$q) {
	++$j;
	my $name = "Cl_".$p."_".$q;
	$outputFile->dataset($name)->set($Clfs->(($j),:));
    }
}

exit;
