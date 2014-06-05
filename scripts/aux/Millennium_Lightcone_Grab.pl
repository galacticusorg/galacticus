#!/usr/bin/env perl
#use lib "./perlLocal";
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use Astro::Cosmology;
use Switch;
use XML::Simple;

# Extract a set of trees from the Millennium Simulation database that are required to populate a lightcone.
# Andrew Benson (20-January-2011)

# Get the data directory.
die("Usage: Millennium_Lightcone_Grab.pl <dataDirectory> <surveySize> <maxRedshift> [options....]") 
    unless ( scalar(@ARGV) >= 3 );
my $dataDirectory = $ARGV[0];
my $surveySize    = $ARGV[1];
my $maxRedshift   = $ARGV[2];

# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
while ( $iArg < scalar(@ARGV)-1 ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Make data directory.
system("mkdir -p ".$dataDirectory);

# Specify user and password.
my $sqlUser     = "";
my $sqlPassword = "";
$sqlUser     = $arguments{"user"    } if ( exists($arguments{"user"    }) );
$sqlPassword = $arguments{"password"} if ( exists($arguments{"password"}) );

# Number of trees to retrieve at one time.
my $treeBlockCount = 200;
$treeBlockCount = $arguments{"treesPerFile"} if ( exists($arguments{"treesPerFile"}) );

# Define constants.
my $degreesToRadians = pdl 3.1415927/180.0;

# Attempt to read saved geometry data.
my $xml;
my $geometryData;
my $gotGeometry;
if ( -e $dataDirectory."/geometry.xml" ) {
    $xml          = new XML::Simple;
    $geometryData = $xml->XMLin($dataDirectory."/geometry.xml");
    $gotGeometry  = 1;
} else {
    $gotGeometry  = 0;
}

# Define the size of the simulation box.
my $boxLength = pdl 500.0; # Mpc/h.
if ( $gotGeometry == 0 ) {
    @{$geometryData->{'boxLength'}} = $boxLength->list();
    $geometryData->{'units'}->{'length'}->{'hubbleExponent'} = -1;
    $geometryData->{'units'}->{'length'}->{'unitsInSI'}      = 3.08568025e22;
}

# Define cosmology.
my $omega0  = 0.25;
my $lambda0 = 0.75;
my $H0      = 100.0; # Use H_0=100km/s/Mpc since Millennium Simulation units are in "little h" system.
my $cosmology = Astro::Cosmology->new(
				      omega_matter => $omega0,
				      omega_lambda => $lambda0,
				      H0           => $H0
				      );
@{$geometryData->{'cosmology'}->{'parameter'}} = (
    {name => 'omega0' , value => $omega0 },
    {name => 'lambda0', value => $lambda0},
    {name => 'H0'     , value => $H0     }
    ) if ( $gotGeometry == 0);

# Define minimum size of survey field.
my $surveySizeRadians = $surveySize*$degreesToRadians;
my $invAngle          = 1.0/($surveySize*$degreesToRadians);
my $halfAngle         = 0.5*$surveySize*$degreesToRadians;
my $tanHalfAngle      = tan($halfAngle);
if ( $gotGeometry == 0 ) {
    @{$geometryData->{'fieldOfView'}->{'length'}} = $surveySizeRadians->list();
    $geometryData->{'fieldOfView'}->{'geometry'} = "square";
}

# Define maximum redshift of lightcone.
my $maxDistance =  $cosmology->comov_dist($maxRedshift);
@{$geometryData->{'maximumDistance'}} = $maxDistance->list() if ( $gotGeometry == 0);

# Redshifts corresponding to Millennium Simulation outputs.
my $millenniumSnapshots = pdl ( 
				0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
				32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 
				);
my $millenniumRedshifts = pdl (
			       127.0000000, 79.99789400, 49.99959000, 30.00006300, 19.91569000, 18.24372300, 16.72452500, 15.34307400, 
			       14.08591400, 12.94078000, 11.89657000, 10.94386400, 10.07346100, 9.277915000, 8.549912000, 7.883203500, 
			       7.272188000, 6.711586500, 6.196833600, 5.723864000, 5.288833600, 4.888449000, 4.519555600, 4.179468600, 
			       3.865682800, 3.575905000, 3.308097800, 3.060419000, 2.831182700, 2.618861400, 2.422044000, 2.239485500, 
			       2.070027400, 1.912632700, 1.766335800, 1.630270700, 1.503636500, 1.385718100, 1.275846200, 1.173416900, 
			       1.077874500, 0.988708140, 0.905462400, 0.827699100, 0.755035640, 0.687108800, 0.623590100, 0.564176600, 
			       0.508591400, 0.456577240, 0.407899440, 0.362340270, 0.319703430, 0.279801800, 0.242469090, 0.207548630, 
			       0.174897610, 0.144383420, 0.115883370, 0.089287830, 0.064493395, 0.041403063, 0.019932542, 0.000000000
			       );
my $millenniumDistances    = pdl zeroes(nelem($millenniumRedshifts));
my $millenniumDistancesMin = pdl zeroes(nelem($millenniumRedshifts));
my $millenniumDistancesMax = pdl zeroes(nelem($millenniumRedshifts));
for(my $i=0;$i<nelem($millenniumRedshifts);++$i) {
    my @z = $millenniumRedshifts->index($i)->list();
    $millenniumDistances->index($i) .= $cosmology->comov_dist($z[0]);
}
for(my $i=0;$i<nelem($millenniumRedshifts);++$i) {
    if ( $i == 0 ) {
 	$millenniumDistancesMax->index($i) .= $millenniumDistances->index($i);
    } else {
 	$millenniumDistancesMax->index($i) .= 0.5*($millenniumDistances->index($i)+$millenniumDistances->index($i-1));
    }
    if ( $i == nelem($millenniumRedshifts)-1 ) {
 	$millenniumDistancesMin->index($i) .= $millenniumDistances->index($i);
    } else {
 	$millenniumDistancesMin->index($i) .= 0.5*($millenniumDistances->index($i)+$millenniumDistances->index($i+1));
    }
}
@{$geometryData->{'outputs'}->{'redshift'       }} = $millenniumRedshifts   ->list() if ( $gotGeometry == 0);
@{$geometryData->{'outputs'}->{'minimumDistance'}} = $millenniumDistancesMin->list() if ( $gotGeometry == 0);
@{$geometryData->{'outputs'}->{'maximumDistance'}} = $millenniumDistancesMax->list() if ( $gotGeometry == 0);

# Compute comoving distances corresponding to Millennium Redshifts.

# Compute comoving distance/redshift relation.
my $redshiftCount     = 1000;
my $redshifts         = pdl (0..$redshiftCount)*$maxRedshift/$redshiftCount; 
my $comovingDistances = pdl zeroes(nelem($redshifts));
my @zs                = $redshifts->list();
for(my $i=0;$i<nelem($redshifts);++$i) {
    my $z =$zs[$i];
    if ( $z <= 0.0 ) {
 	$comovingDistances->index($i) .= 0.0;
    } else {
 	$comovingDistances->index($i) .= $cosmology->comov_dist($z);
    }
}

# Compute the values of m and n. Here we assume n=m+1 (to keep n and m as small as
# possible. According to Kitzbichler & White (2007; MNRAS; 376; 2) this will give a field of
# dimensions 1/m^2/ by 1/m/n^2 radians. Given the inverse anglular size of the survey "invAngle"
# we can solve the resulting cubic equation for m.
my $m = int((8.0+108.0*$invAngle+12.0*sqrt(12.0*$invAngle+81.0*$invAngle**2))**(1.0/3.0)/6.0+2.0/3.0*(8.0+108.0*$invAngle+12.0*sqrt(12.0*$invAngle+81.0*$invAngle**2))**(-1.0/3.0)-2.0/3.0);
my $n = $m+1;
my $angle1   = 1.0/$m/$n/$n/$degreesToRadians;
my $angle2   = 1.0/$m/$m/$n/$degreesToRadians;
my $maxDepth = sqrt($n**2+$m**2+($n*$m)**2)*$boxLength;
(my $maxRedshifts,my $error) = interpolate($maxDepth,$comovingDistances,$redshifts);
print "Integer divisors (m,n)         = (".$m.",".$n.")\n";
print "Request survey angle           = ".$surveySize." degrees\n";
print "Actual survey angles           = ".$angle1." x ".$angle2." degrees\n";
print "Maximum depth before repeats   = ".$maxDepth." Mpc/h\n";
print "Corresponding maximum redshift = ".$maxRedshifts."\n";

# Choose origin.
my $origin;
if ( $gotGeometry == 0 ) {
    $origin = pdl random(3)*$boxLength;
    @{$geometryData->{'origin'}->{'coordinate'}} = $origin->list();
} else {
    $origin = pdl $geometryData->{'origin'}->{'coordinate'};
}
print "Survey origin selected to be   = ".$origin."\n";

# Compute unit vectors.
my $unitVector1;
my $unitVector2;
my $unitVector3;
if ( $gotGeometry == 0 ) {
    my $xAxis       = pdl (1.0,0.0,0.0);
    $unitVector1 = pdl zeroes(3);
    $unitVector1->index(0) .= $n   /sqrt($n**2+$m**2+($n*$m)**2);
    $unitVector1->index(1) .= $m   /sqrt($n**2+$m**2+($n*$m)**2);
    $unitVector1->index(2) .= $m*$n/sqrt($n**2+$m**2+($n*$m)**2);
    $unitVector2 = crossp($unitVector1,$xAxis      );
    $unitVector3 = crossp($unitVector1,$unitVector2);
    $unitVector1 /= sqrt(sum($unitVector1**2));
    $unitVector2 /= sqrt(sum($unitVector2**2));
    $unitVector3 /= sqrt(sum($unitVector3**2));
    @{$geometryData->{'unitVector1'}->{'coordinate'}} = $unitVector1->list();
    @{$geometryData->{'unitVector2'}->{'coordinate'}} = $unitVector2->list();
    @{$geometryData->{'unitVector3'}->{'coordinate'}} = $unitVector3->list();
} else {
    $unitVector1 = pdl $geometryData->{'unitVector1'}->{'coordinate'};
    $unitVector2 = pdl $geometryData->{'unitVector2'}->{'coordinate'};
    $unitVector3 = pdl $geometryData->{'unitVector3'}->{'coordinate'};
}
print "Unit vectors:\n";
print "\t".$unitVector1."\n";
print "\t".$unitVector2."\n";
print "\t".$unitVector3."\n";

# Output the lightcone geometry data.
if ( $gotGeometry == 0 ) {
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"geometry");
    open(oHndl,">".$dataDirectory."/geometry.xml");
    print oHndl $xmlOutput->XMLout($geometryData,KeyAttr => {});
    close(oHndl);
}

# Open a file to output cone data to.
unless ( -e $dataDirectory."/treeIDs.data" ) {
    
    # Open a file to output tree data to.
    open(tHndl,">".$dataDirectory."/treeIDs.data");
    my %treeIDs;
        
    # Loop through redshifts.
    for(my $i=0;$i<nelem($millenniumRedshifts);++$i) {
	# Check that this snapshot contributes to the requested redshift range.
	if ( $millenniumDistancesMin->index($i) <= $maxDistance ) {
	    # Determine maximum width of the cone in this interval.
	    my $coneWidth = $millenniumDistancesMax->index($i)*$surveySize*$degreesToRadians*sqrt(2.0);
	    
	    # Determine the range of periodic boxes that must be explored.
	    my $periodicRange = pdl zeroes(3,2);
	    $periodicRange(:,(0)) .= floor(($origin+$millenniumDistancesMin->index($i)*$unitVector1)/$boxLength)-1;
	    $periodicRange(:,(1)) .= ceil (($origin+$millenniumDistancesMax->index($i)*$unitVector1)/$boxLength)+0;
	    
	    # Loop over periodic boxes and determine if the lightcone is in the box.
	    my $xReplicants = pdl $periodicRange((0),(0))..$periodicRange((0),(1));
	    my $yReplicants = pdl $periodicRange((1),(0))..$periodicRange((1),(1));
	    my $zReplicants = pdl $periodicRange((2),(0))..$periodicRange((2),(1));
	    
	    foreach my $ix ( $xReplicants->list() ) {
		foreach my $iy ( $yReplicants->list() ) {
		    foreach my $iz ( $zReplicants->list() ) {
			# Loop over each face of this cube.
			my $intersectsCone = 0;
			my $belowMin = 0;
			my $aboveMax = 0;
			for (my $iFace=0;$iFace<6;++$iFace) {
			    my $faceOrigin  = pdl (0,0,0);
			    my $faceVector1 = pdl (0,0,0);
			    my $faceVector2 = pdl (0,0,0);
			    if ( $iFace == 0 ) {
				$faceOrigin  .= pdl ($ix,$iy,$iz);
				$faceVector1 .= pdl (1,0,0);
				$faceVector2 .= pdl (0,1,0);
			    }
			    if ( $iFace == 1 ) {
				$faceOrigin  .= pdl ($ix,$iy,$iz);
				$faceVector1 .= pdl (1,0,0);
				$faceVector2 .= pdl (0,0,1);
			    }
			    if ( $iFace == 2 ) {
				$faceOrigin  .= pdl ($ix,$iy,$iz);
				$faceVector1 .= pdl (0,1,0);
				$faceVector2 .= pdl (0,0,1);
			    }
			    if ( $iFace == 3 ) {
				$faceOrigin  .= pdl ($ix+1,$iy,$iz);
				$faceVector1 .= pdl (0,1,0);
				$faceVector2 .= pdl (0,0,1);
			    }
			    if ( $iFace == 4 ) {
				$faceOrigin  .= pdl ($ix,$iy+1,$iz);
				$faceVector1 .= pdl (1,0,0);
				$faceVector2 .= pdl (0,0,1);
			    }
			    if ( $iFace == 5 ) {
				$faceOrigin  .= pdl ($ix,$iy,$iz+1);
				$faceVector1 .= pdl (1,0,0);
				$faceVector2 .= pdl (0,1,0);
			    }
			    # Compute coordinates of face mid-point.
			    my $faceMidpoint = ($faceOrigin+0.5*$faceVector1+0.5*$faceVector2)*$boxLength;
			    # Construct face normal vector.
			    my $faceNormal = crossp($faceVector1,$faceVector2);
			    # Solve linear equations to get solution for intersection point.
			    my $B  = pdl [$origin-$faceOrigin*$boxLength];
			    my $BT = t($B);
			    my $A  = pdl [$faceVector1,$faceVector2,-$unitVector1];
			    my $AT = t($A);
			    my $solution = msolve($AT,$BT);
			    my $intersectionDistance = $solution((0),(2));
			    my $intersectionPoint = $origin+$intersectionDistance*$unitVector1;
			    
			    # Compute distance from face midpoint along both axes.
			    my $faceDistance1 = (($intersectionPoint-$faceMidpoint)*$faceVector1)->sum();
			    my $faceDistance2 = (($intersectionPoint-$faceMidpoint)*$faceVector2)->sum();
			    
			    # Check if any part of the cone can intersect the cube.
			    my $cosineAngle = abs(inner($faceNormal,$unitVector1));
			    if (    abs($faceDistance1) < 0.5*$boxLength+$coneWidth/$cosineAngle
				    && abs($faceDistance2) < 0.5*$boxLength+$coneWidth/$cosineAngle ) {
				if (  $intersectionDistance >= $millenniumDistancesMin->index($i)
				      && $intersectionDistance <= $millenniumDistancesMax->index($i) ) {
				    $intersectsCone = 1;
				}
				if ( $intersectionDistance < $millenniumDistancesMin->index($i) ) {
				    $belowMin = 1;
				    if ( $aboveMax == 1 ) {$intersectsCone = 1};
				}
				if ( $intersectionDistance > $millenniumDistancesMax->index($i) ) {
				    $aboveMax = 1;
				    if ( $belowMin == 1 ) {$intersectsCone = 1};
				}
			    }
			}
			
			# If this cube intersects the cone, then process.
			if ( $intersectsCone == 1 ) {
			    print "Processing:\n";
			    print "\tSnapshot: ".$millenniumSnapshots->index($i)."\n";
			    print "\tReplicant: ".$ix.":".$iy.":".$iz."\n";
			    
			    # Construct equations giving distances along lightcone axes.
			    my $lightConeX = "x+".$ix."*".$boxLength."-".$origin->index(0);
			    my $lightConeY = "y+".$iy."*".$boxLength."-".$origin->index(1);
			    my $lightConeZ = "z+".$iz."*".$boxLength."-".$origin->index(2);
			    my $lightCone1 = 
				"(".$lightConeX.")*".$unitVector1->index(0)
				."+(".$lightConeY.")*".$unitVector1->index(1)
				."+(".$lightConeZ.")*".$unitVector1->index(2);
			    my $lightCone2 = 
				"(".$lightConeX.")*".$unitVector2->index(0)
				."+(".$lightConeY.")*".$unitVector2->index(1)
				."+(".$lightConeZ.")*".$unitVector2->index(2);
			    my $lightCone3 = 
				"(".$lightConeX.")*".$unitVector3->index(0)
				."+(".$lightConeY.")*".$unitVector3->index(1)
				."+(".$lightConeZ.")*".$unitVector3->index(2);
			    
			    # Construct geometrical and temporal condition on halo selection.
			    my $condition1 = $lightCone1." > ".$millenniumDistancesMin->index($i);
			    my $condition2 = $lightCone1." <= ".$millenniumDistancesMax->index($i);
			    my $condition3 = "abs(".$lightCone2.")/(".$lightCone1.") < ".$tanHalfAngle;
			    my $condition4 = "abs(".$lightCone3.")/(".$lightCone1.") < ".$tanHalfAngle;
			    my $condition5 = "snapNum = ";
			    unless ( exists($arguments{"fixedSnapshot"}) ) {
				$condition5 .= $millenniumSnapshots->index($i);
			    } else {
				$condition5 .= $arguments{"fixedSnapshot"};
			    }
			    
			    # Combine all conditions.
			    my $condition = 
				"where ".$condition1
				." and ".$condition2
				." and ".$condition3
				." and ".$condition4
				." and ".$condition5;
			    
			    # Construct SQL query.
			    $condition =~ s/\+/%2B/g;
			    $condition =~ s/\//%2F/g;
			    my $sqlQuery = "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=select treeId,snapNum,x,y,z from MPAHaloTrees..MHalo ".$condition;
			    my $outputFile = "tmp.data";
			    my $getCommand = "wget --http-user=".$sqlUser." --http-passwd=".$sqlPassword." \"".$sqlQuery."\" -O ".$outputFile;
			    system($getCommand);
			    my $lineCount = 0;
			    my $gotOK     = 0;
			    open(iHndl,$outputFile);
			    while ( my $line = <iHndl> ) {
				if ( $line =~ m/\#OK/ ) {$gotOK = 1};
				if ( $line =~ m/^\d/ ) {
				    chomp($line);
				    ++$lineCount;
				    my @columns = split(/,/,$line);
				    $treeIDs{$columns[0]} = 1;
				}
			    }
			    close(iHndl);
			    
			}
		    }
		}
	    }
	}
    }
    
    # Output tree IDs.
    foreach my $treeID ( sort(keys(%treeIDs)) ) {
	print tHndl $treeID."\n";
    }

    # Close file.
    close(tHndl);
}

# Download merger trees from the Millennium Database.

# Read in the list of tree IDs.
print "Reading tree IDs from file...\n";
my %treeIDs;
open(tHndl,$dataDirectory."/treeIDs.data");
while ( my $line = <tHndl> ) {
    chomp($line);
    $treeIDs{$line} = 1;
}
close(tHndl);
print "   ...done\n";
my @treeIDs        = sort(keys(%treeIDs));
my $temporaryFile  = "tmp.csv";
system("make Millennium_Merger_Tree_File_Maker.exe");
for(my $iTree=0;$iTree<scalar(@treeIDs);$iTree+=$treeBlockCount) {
    my $jTree = $iTree+$treeBlockCount-1;
    if ( $jTree > scalar(@treeIDs)-1 ) {$jTree = scalar(@treeIDs)-1};
    my $select = "root.treeId in (".join(",",@treeIDs[$iTree..$jTree]).")";
    
    # Define the name of the tree file.
    my $treeFile = $dataDirectory."/Lightcone_Trees_".$iTree.":".$jTree.".hdf5";

    # Check if we already have the file.
    unless ( -e $treeFile ) {
	print "Retrieving data for ".$treeFile."...\n";
	
	# Write the SQL selection to a file.
	open(oHndl,">sqlSelect.txt");
	print oHndl $select."\n";
	close(oHndl);

	# Run the script to grab the trees.
	system("scripts/aux/Millennium_Trees_Grab.pl --user ".$sqlUser." --password ".$sqlPassword." --table MPAHaloTrees..MHalo --traceParticles no --selectFile sqlSelect.txt --output ".$temporaryFile);
	
	# Process into Galacticus format.
	system("./Millennium_Merger_Tree_File_Maker.exe ".$temporaryFile." none ".$treeFile." galacticus 1");
	die("Millennium_Lightcone_Grab.pl: failed to build merger tree file")
	    unless ( $? == 0 );
	unlink($temporaryFile);

	print "   ...done\n";
    }

}

exit;
