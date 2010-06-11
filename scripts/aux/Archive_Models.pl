#!/usr/bin/env perl
use File::Find;
use File::Copy;
use XML::Simple;
use Data::Dumper;

# Archives Galacticus models:
#  Stores the full list of parameters used to run the model.
#  Stores the version and revision information.
#  Stores the run time.
#  Stores any fitting data.
#  Stores and comments associated with the model.
#  Places plots into a tar archive and stores them along with a link to them. 
# Andrew Benson (11-June-2010)

if ( $#ARGV != 1 ) {die("Usage: Archive_Models.pl <modelDirectory> <archiveDirectory>")};
$modelDirectory[0] = $ARGV[0];
$archiveDirectory  = $ARGV[1];
unless ( $archiveDirectory =~ m/^\// ) {
    $pwd = `pwd`;
    chomp($pwd);
    $archiveDirectory = $pwd."/".$archiveDirectory;
}

# See if we can find the bzr revision number.
if ( -e ".bzr/branch/last-revision" ) {
    open(bzrHndl,".bzr/branch/last-revision");
    $line = <bzrHndl>;
    if ( $line =~ m/^\s*(\d+)/ ) {
	$bzrRevision = $1;
    } else {
	$bzrRevision = "unknown";
    }
}

# Find and process all files in the source directory tree.
find(\&processFile, @modelDirectory);

exit;

sub processFile {

    # Process a file in the source directory tree. If it's a Galacticus model file, archive it.
    $fileName = $_;
    chomp($fileName);

    # Check if this is an HDF5 file.
    if ( $fileName =~ m/\.hdf5(\.bz2)??$/ ) {

	# Get current time to use as filename.
	$now = `date +%Y.%m.%d.%H.%M.%S`;
	chomp($now);

	# Unpack if necessary.
	if ( $fileName =~ m/\.hdf5\.bz2$/ ) {
	    system("bunzip2 ".$fileName);
	    $fileName =~ s/\.bz2$//;
	    $rePackFile = 1;
	} else {
	    $rePackFile = 0;
	}

	# Begin a data structure.
	$data->{'contents'} = "Archived Galacticus model.";

	# Extract all parameters from the file.
	open(pHndl,"h5ls -d ".$fileName."/Parameters|");
	$iParameter = -1;
	while ( $line = <pHndl> ) {
	    @columns = split(/\s+/,$line);
	    $name = $columns[0];
	    $line = <pHndl>;
	    $line = <pHndl>;
	    if ( $line =~ m/^\s*\(\d+\)\s*(.*)/ ) {
		$value = $1;
		$value =~ s/\s*\"\s*,\s*\"\s*/ /g;
		$value =~ s/\s*\,\s*/ /g;
		$value =~ s/\s*\"\s*//g;
		$value =~ s/&/&amp;/g;
		$value =~ s/\s*\'\s+\'\s+repeats\s+\d+\s+times\s*$//;
		$parameter->{'name'}  = $name;
		$parameter->{'value'} = $value;
		${$data->{'parameters'}->{'parameter'}}[++$iParameter] = $parameter;
		undef($parameter);
	    }
	}
	close(pHndl);

	# Extract version number and run time from the file.
	foreach $version ( "versionMajor", "versionMinor", "versionRevision", "runTime" ) {
	    open(pHndl,"h5ls -d -S ".$fileName."/Version/".$version."|");
	    $line = <pHndl>;
	    $line = <pHndl>;
	    $line = <pHndl>;
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    $line =~ s/\%20/ /g;
	    $line =~ s/\"//g;
	    close(pHndl);
	    $data->{'version'}->{$version} = $line;
	}
	$data->{'version'}->{'bzrRevision'} = $bzrRevision;

	# Looks for comments in the file and include them if necessary.
	if ( -e "comments.txt" ) {
	    open(cHndl,"comments.txt");
	    while ( $line = <cHndl> ) {
		$data->{'comments'} .= $line;
	    }
	    close(cHndl);
	    $data->{'comments'} =~ s/\n$//;
	}

	# Include any fitting results.
	if ( -e "galacticusFits.xml.bz2" ) {
	    system("bunzip2 galacticusFits.xml.bz2");
	    $reCompressFits = 1;
	} else {
	    $reCompressFits = 0;
	}
	if ( -e "galacticusFits.xml" ) {
	    $xml = new XML::Simple;
	    $fitData = $xml->XMLin("galacticusFits.xml");
	    $data->{'galacticusFits'} = $fitData;
	    system("bzip2 galacticusFits.xml") if ( $reCompressFits == 1 );
	}
	
	# Search for any plots in this directory.
	undef(@plotFiles);
	undef(@recompressFiles);
	opendir(mDir,".");
	while ( $plotFile = readdir(mDir) ) {
	    if ( $plotFile =~ m/\.pdf(\.bz2)??$/ ) {
		if ( $plotFile =~ m/\.pdf\.bz2$/ ) {
		    system("bunzip2 ".$plotFile);
		    $plotFile =~ s/\.bz2$//;
		    push(@recompressFiles,$plotFile);
		}
		push(@plotFiles,$plotFile);
	    }
	}
	closedir(mDir);
	if ( $#plotFiles >= 0 ) {
	    $tarFile = $now.".tar.bz2";
	    system("tar cvfj ".$tarFile." ".join(" ",@plotFiles));
	    system("mkdir -p ".$archiveDirectory."/plots");
	    move($tarFile,$archiveDirectory."/plots/".$tarFile);
	    $data->{'plotsArchive'} = $archiveDirectory."/plots/".$tarFile;
	    system("bzip2 ".join(" ",@recompressFiles)) if ( $#recompressFiles >= 0 );
	}

	# Repack the main file if necessary.
	system("bzip2 ".$fileName) if ( $rePackFile == 1 );
	
	# Output the data to an XML file.
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"archivedModel");
	system("mkdir -p ".$archiveDirectory."/models");
	open(outHndl,">".$archiveDirectory."/models/".$now.".xml");
	print outHndl $xmlOutput->XMLout($data);
	close(outHndl);
	undef($data);
	system("bzip2 ".$archiveDirectory."/models/".$now.".xml");
	
    }

}
