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
unshift(@INC, $galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;

# Construct a mangle mask file for the UltraVISTA survey of Muzzin et al. (2013; http://adsabs.harvard.edu/abs/2013ApJ...777...18M)
# Andrew Benson (27-August-2014)

# Define work directory.
my $dataDirectoryName = $galacticusPath."constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/";
# Define work directory.
my $workDirectoryName = $dataDirectoryName."work/";
# Define mangle directory.
my $mangle            = $galacticusPath."aux/mangle/bin/";
# Generate raw field geometry.
unless ( -e $dataDirectoryName."field.ply" ) {
    open(my $rawFile,">".$dataDirectoryName."field.ply");
    print $rawFile "149.373 150.779 1.604 2.81\n";
    close($rawFile);    
}
# Specify radii for bright and medium stars.
my $brightRadius = 75.0/3600.0;
my $mediumRadius = 30.0/3600.0;
# Generate star mask holes.
unless ( -e $dataDirectoryName."stars.ply" ) {
    open(my $starFile,">".$dataDirectoryName."stars.ply");
    # Read USNO star list.
    (my $id_usno, my $ra_usno, my $dec_usno, my $b1, my $b1_sg, my $r1, my $r1_sg, my $b2, my $b2_sg, my $r2, my $r2_sg, my $i2, my $i2_sg) 
	= rcols($dataDirectoryName."USNO_star_list",0,1,2,3,4,5,6,7,8,9,10,11,12);
    # Find USNO bright and medium stars.
    my $usnoBright = which(($b1 >  0.0) & ($b1 <= 10.0));
    my $usnoMedium = which(($b1 > 10.0) & ($b1 <= 13.0));
    # Generate masks for USNO bright stars.
    for(my $i=0;$i<nelem($usnoBright);++$i) {
	print $starFile $ra_usno->($usnoBright)->(($i))."\t".$dec_usno->($usnoBright)->(($i))."\t".$brightRadius."\n";
    }
    # Generate masks for USNO medium stars.
    for(my $i=0;$i<nelem($usnoMedium);++$i) {
	print $starFile $ra_usno->($usnoMedium)->(($i))."\t".$dec_usno->($usnoMedium)->(($i))."\t".$mediumRadius."\n";
    }
    # Read 2MASS stars list.
    (my $ra_2mass, my $dec_2mass, my $_2massj, my $ej_2mass, my $jsig_2mass, my $h_2mass, my $eh_2mass, my $hsig_2mass, my $k_2mass, my $ek_2mass, my $ksig_2mass, my $phot_2mass) = rcols($dataDirectoryName."2mass_psc_list.dat",0,1,2,3,4,5,6,7,8,9,10,11,{EXCLUDE => qr/^[\\|]/});
    # Find 2MASS bright and medium stars.
    my $twomassBright = which(($k_2mass >  0.0) & ($k_2mass <=  8.0));
    my $twomassMedium = which(($k_2mass >  8.0) & ($k_2mass <= 10.5));
    # Generate masks for 2MASS bright stars.
    for(my $i=0;$i<nelem($twomassBright);++$i) {
	print $starFile $ra_2mass->($twomassBright)->(($i))."\t".$dec_2mass->($twomassBright)->(($i))."\t".$brightRadius."\n";
    }
    # Generate masks for 2MASS medium stars.
    for(my $i=0;$i<nelem($twomassMedium);++$i) {
	print $starFile $ra_2mass->($twomassMedium)->(($i))."\t".$dec_2mass->($twomassMedium)->(($i))."\t".$mediumRadius."\n";
    }
    close($starFile);
}
# Generate bad pixel map.
system("wget http://ultravista.org/release1/data/UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits -O ".$workDirectoryName."UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits")
    unless ( -e $workDirectoryName."UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits" );
system("python ".$dataDirectoryName."badPixels.py ".$workDirectoryName)
    unless ( -e $workDirectoryName."UVISTA_Ks_15_12_10_bad.ply" );
# Pixelize the survey field.
system($mangle."pixelize -Ps0,11 -ir ".$dataDirectoryName."field.ply ".$dataDirectoryName."fieldP.ply")
    unless ( -e $dataDirectoryName."fieldP.ply" );
# Pixelize the stars.
system($mangle."pixelize -Ps0,11 -ic ".$dataDirectoryName."stars.ply ".$dataDirectoryName."starsP.ply")
    unless ( -e $dataDirectoryName."starsP.ply" );
# Pixelize the bad pixels.
system($mangle."pixelize -Ps0,11 -ir /".$workDirectoryName."UVISTA_Ks_15_12_10_bad.ply /".$workDirectoryName."UVISTA_Ks_15_12_10_badP.ply")
    unless ( -e "/".$workDirectoryName."UVISTA_Ks_15_12_10_badP.ply" );
# Snap the survey field.
system($mangle."snap ".$dataDirectoryName."fieldP.ply ".$dataDirectoryName."fieldPS.ply")
    unless ( -e $dataDirectoryName."fieldPS.ply" );
# Snap the stars,
system($mangle."snap ".$dataDirectoryName."starsP.ply ".$dataDirectoryName."starsPS.ply")
    unless ( -e $dataDirectoryName."starsPS.ply" );
# Snap the bad pixels.
system($mangle."snap -a1.0e-7d -b1.0e-7d -t1.0e-7d /".$workDirectoryName."UVISTA_Ks_15_12_10_badP.ply /".$workDirectoryName."UVISTA_Ks_15_12_10_badPS.ply")
    unless ( -e "/".$workDirectoryName."UVISTA_Ks_15_12_10_badPS.ply" );
# Holeize the stars.
unless ( -e $dataDirectoryName."holesPS.ply" ) {
    open(my $starsIn ,    $dataDirectoryName."starsPS.ply");
    open(my $holesOut,">".$dataDirectoryName."holesPS.ply");
    while ( my $line = <$starsIn> ) {
	$line =~ s/1 weight/0 weight/;
	print $holesOut $line;
    }
    close($starsIn );
    close($holesOut);
}
# Holeize the bad pixels.
unless ( -e "/".$workDirectoryName."badHolesPS.ply" ) {
    open(my $badPixelsIn,    "/".$workDirectoryName."UVISTA_Ks_15_12_10_badPS.ply");
    open(my $holesOut   ,">/".$workDirectoryName."badHolesPS.ply" );
    while ( my $line = <$badPixelsIn> ) {
	$line =~ s/1 weight/0 weight/;
	print $holesOut $line;
    }
    close($badPixelsIn );
    close($holesOut);
}
# Balkanize the survey and holes.
system($mangle."balkanize ".$dataDirectoryName."fieldPS.ply ".$dataDirectoryName."holesPS.ply /".$workDirectoryName."badHolesPS.ply /".$workDirectoryName."mask.ply")
    unless ( -e "/".$workDirectoryName."mask.ply" );
# Unify the mask.
system($mangle."unify /".$workDirectoryName."mask.ply /".$dataDirectoryName."surveyMask.ply")
    unless ( -e "/".$dataDirectoryName."surveyMask.ply" );

exit;
