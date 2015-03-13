#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V094'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V094'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");

use 5.010;  # so filetest ops can stack
use Tie::File;
use Fcntl qw( SEEK_SET );

# Script to extract forests of merger trees from the Bolshoi simulation as input to Galacticus.
# Stephanie Dörschner (05-December-2014)

# Get arguments.
 die "Usage: Bolshoi_Trees_Grab.pl <workingDirectory> <outputFile> [forestFile=<forestFile>] [locationFile=<locationFile>] [ignoreFile=<ignoreFile>] [treeCountMaximum=<treeCountMaximum>] [startLine=<startLine>] \n"
     unless( scalar( @ARGV ) >= 2 && scalar( @ARGV ) <= 7 );
     my $workDir        = shift @ARGV	 ;
     my $outFile	= shift @ARGV	 ;
     my $forestFile	= "forests.list" ;
     my $locationFile	= "locations.dat";
     my $ignoreFile     = "ignore.dat"   ;
     my $append		= 0		 ;
     my $treeCountMaximum = 10000	 ;
     my $startLine	= 1		 ;

 my $lastChar = substr( $workDir,length($workDir)-1, 1 );
 if ( $lastChar eq '/' )
 {
    chop( $workDir );
 }

SWITCH: foreach my $arg ( @ARGV )
 {
    my @opt = split( "=", $arg );
    if( $opt[0] eq "forestFile"		) { $forestFile   = $opt[1]; 		  next SWITCH; }
    if( $opt[0] eq "locationFile"	) { $locationFile = $opt[1]; 		  next SWITCH; }
    if( $opt[0] eq "ignoreFile"   	) { $ignoreFile   = $opt[1]; $append = 1; next SWITCH; }
    if( $opt[0] eq "treeCountMaximum" 	) { $treeCountMaximum = $opt[1]; 	  next SWITCH; }
    if( $opt[0] eq "startLine"    	) { $startLine    = $opt[1]; 		  next SWITCH; }
 }

# Find forests to copy.
 
    my %forestHash;
    my %ignoreHash;
    my $treeCount  = 0;
    my $forestCount= 0;

    if ( -f $workDir."/".$ignoreFile )
    {
	# Make hash of forests to ignore.
	open( IFORESTS, "<", "$workDir/$ignoreFile" ) or die "Could not open $workDir/$ignoreFile: $!\n";
	while( my $ignoreLine = <IFORESTS> )
	{
	    $ignoreHash{ 'ignoreLine' } = 1;
	}
	close( IFORESTS ) or die "Could not close $workDir/$ignoreFile: $!\n";
    }


    # If needed download forests.list and locations.dat.
    if ( not  -f -s $workDir/$forestFile ) 
    {
	my $url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/forests.list.gz";
	system( "wget", "-e robots=off", "--quiet", "--directory-prefix=$workDir", "$url" ) == 0 or die "Could not download forests.list.gz from $url: $!\n";
	system( "gunzip", "$workDir/forests.list.gz") == 0  or die "Could not unzip forests.list.gz: $!\n";
	$forestFile = "forests.list";
    }
    if ( not  -f -s $workDir/$locationFile )
    {
	
        my $url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/locations.dat.gz";
        system( "wget", "-e robots=off", "--quiet", "--directory-prefix=$workDir", "$url" ) == 0 or die "Could not download locations.dat.gz from $url: $!\n";
	system( "gunzip", "$workDir/locations.dat.gz") == 0  or die "Could not unzip locations.dat.gz: $!\n";
	$locationFile = "locations.dat";
    }

    # Collect forests.
    open( FORESTS,  "<", "$workDir/$forestFile" ) or die "Could not open $workDir/$forestFile: $!\n";

    # Wrapper loop to find several forests.
    while( $treeCount < $treeCountMaximum )
    {
	# Reset Filepointer and line count.
	seek( FORESTS, 0, SEEK_SET);
	my $lineCount = 0;

	# Find forest to record.
FOREST:	while( my $forestLine = <FORESTS> )
	{
	  ++$lineCount;

	  # Skip header and offset.
	  if ( $lineCount == 1 || $lineCount < $startLine ) { next FOREST; }

	  # Extract forest ID.
	  my @forestID = split( " ", $forestLine );
          my $forest = $forestID[1];

	  # Skip forests that have already been recorded or should be ignored.
	  if ( exists( $forestHash{ $forest } ) || exists( $ignoreHash{ $forest } ) ) { next FOREST; }
		
	  # Record tree ID belonging to forest, build hash of hashes.
	  $forestHash{ $forest }{ $forestID[0] } = [];
	  ++$treeCount;
	
	  # Find all other trees belonging to $forest.
	  while( $forestLine = <FORESTS> )
	  {
	    my @forestID = split( " ", $forestLine );
	    if ( $forestID[1] == $forest )
	    {
		# Record tree ID belonging to forest.
		$forestHash{ $forest }{ $forestID[0] } = [];
		++$treeCount;
	    }
	  }

	  my $numberOfTrees = scalar( keys %{ $forestHash{$forest} } );
	
	  # Find locations of all trees belonging to $forest.
	  open( LOCATIONS, "<", "$workDir/$locationFile" ) or die "Could not open $workDir/$locationFile: $!\n";

	  # Skip head line.
	  my $locationLine = <LOCATIONS>;

	  my $counter = 0;
LOCATION: while ( $locationLine = <LOCATIONS> )
	  {
	    my @locID = split( " ", $locationLine );
	    if ( $forestHash{ $forest }{ $locID[0] } )
	    { 
		# Record file ID, offset and file name.
		$forestHash{ $forest }{ $locID[0] } = [ $locID[2], $locID[3] ];
		++$counter;
	    }
	    last LOCATION if ( $counter == $numberOfTrees );
	  }
	  close( LOCATIONS ) or die "Could not close $workDir/$locationFile: $!\n";

	  print STDOUT "Forest $forest has $numberOfTrees trees.\n";

	  # Count found forest and adjust start line for the next loop.
	  ++$forestCount;
	  $startLine = $lineCount;
	  last FOREST;
	}
    }
    close( FORESTS ) or die "Could not close $workDir/$forestFile: $!\n";

    # Append found forests to $ignoreFile.
    if ( $append == 1 )
    {
	open( IFORESTS, ">>", "$workDir/$ignoreFile" ) or die "Appending forests to ignoreFile: Could not open $workDir/$ignoreFile: $!\n";
	foreach my $forestID ( keys %forestHash )
	{
	    print IFOREST "$forestID\n";
	}
	close ( IFOREST ) or die "Appending forests to ignoreFile: Could not close $workDir/$ignoreFile: $!\n";
    }



########## Copy all found trees to output file in ascending order. ##########


 open( MYTREES, ">", "$workDir/$outFile" ) or die "Could not open $workDir/$outFile: $!\n";
 my $header = 0;

 # Copy forests in ascending order.
 foreach my $forestID ( sort( keys %forestHash ) )
 {
    # Copy trees belonging to $forestID.
    foreach my $treeID ( keys %{ $forestHash{ $forestID } } )
    {
	my $offset   = $forestHash{ $forestID }{ $treeID }[0];
	my $fileName = $forestHash{ $forestID }{ $treeID }[1];

	# Download and unzip treefile if necessary.
	if ( not ( -f -s  "$workDir/$fileName" ) )
	{
	    my $fileGZ = "$fileName.gz";
	    if ( -f -s "$workDir/$fileGZ" )
	    {
		system( "gunzip", "$workDir/$fileGZ" );
		if ( $? != 0 )
		{
		    print STDOUT "Warning: Could not unzip $fileGZ: $!. Downloading $fileGZ.\n";
		    system( "rm", "$workDir/$fileGZ" );
		    my $url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/$fileGZ";
		    system( "wget", "-e robots=off", "--quiet", "--directory-prefix=$workDir", "$url" ) == 0 or die "Could not download $fileGZ from $url: $!\n";
		    system( "gunzip", "$workDir/$fileGZ" ) 						== 0 or die "Could not unzip $fileGZ: $!\n";
		}
	    }
	    else
	    {
		my $url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/$fileGZ";
		system( "wget", "-e robots=off", "--quiet", "--directory-prefix=$workDir", "$url" ) 	== 0 or die "Could not download $fileGZ from $url: $!\n";
                system( "gunzip", "$workDir/$fileGZ" ) 							== 0 or die "Could not unzip $fileGZ: $!\n";
	    }
	}

	# Copy tree to file.
	open( NEWTREES, "<", "$workDir/$fileName" ) or die "Could not open $workDir/$fileName: $!\n";

	 # Copy header if necessary.
	 if ( $header == 0 )
	 {
	    my $counter = 0;
	    while( my $headerLine = <NEWTREES>)
	    {
		++$counter;
		if ( $counter == 1 ) { chomp $headerLine; $headerLine =  join( " ", "#tree_id", substr( $headerLine, 1 ), "\n" ); }
		if ( $counter == 4 ) 
		{ 
		    print MYTREES "#Tree_id: ID of the merger tree, the halo belongs to. Trees of the same forest have identical tree IDs.\n";
		}
		if ( $counter == 45) { print MYTREES "#Number of Trees:", $treeCount , "\n"; last; }
		print MYTREES $headerLine;
	    }
	    $header = 1;
	 }

	 # Find tree data at given offset.
	 seek( NEWTREES, $offset, SEEK_SET );
	
	 # Copy data.
TREE:	 while( my $output = <NEWTREES> )
	 {
	    if ( $output =~ m/#tree $treeID/ ){ next TREE; }
	    if ( $output =~ m/#tree/	     ){ last TREE; }
	    print MYTREES "$forestID $output";
	 }
	close( NEWTREES ) or die "Could not close $workDir/$fileName: $!\n";
    }
 }

 close( MYTREES ) or die "Could not close $workDir/$outFile: $!\n";
 print STDOUT "Copied $treeCount merger trees of $forestCount forests to file $workDir/$outFile, sorted in ascending order.\n";

 exit;
