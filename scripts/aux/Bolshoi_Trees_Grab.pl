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
$ENV{TZ} = 'Europe/London';

use Tie::File;
use Fcntl qw( SEEK_SET );
use POSIX::strftime::GNU;

# Get arguments.
 die "Usage: Bolshoi_Trees_Grab.pl <workingDirectory> <outputFile> [ignoreFile=<ignoreFile>] [treeCountMaximum=<treeCountMaximum>] [appendToIgnoreFile=<append>]\n"
     unless( scalar( @ARGV ) >= 2 && scalar( @ARGV ) <= 5 );
     my $workDir        = shift @ARGV	 ;	# Directory where all files are saved.
     my $outFile	= shift @ARGV	 ;	# Name of ouput file.
     my $ignoreFile     = "ignore.dat"   ;	# Name of file that contains the IDs of forests which should not be extracted.
     my $append		= 0		 ;	# Whether to append IDs of extracted forests to the ignore file.
     my $treeCountMaximum = 10000	 ;	# Approximate number of trees that are to be extracted. When the number is reached the current forest will be completed.

 my $lastChar = substr( $workDir,length($workDir)-1, 1 );
 if ( $lastChar eq '/' )
 {
    chop( $workDir );
 }

SWITCH: foreach my $arg ( @ARGV )
 {
    my @opt = split( "=", $arg );
    if( $opt[0] eq "ignoreFile"   	) { $ignoreFile   	= $opt[1]; 			  next SWITCH; }
    if( $opt[0] eq "treeCountMaximum" 	) { $treeCountMaximum 	= $opt[1]; 	  		  next SWITCH; }
    if( $opt[0] eq "appendToIgnoreFile"	) { if ( ($opt[1] eq "append") || ($opt[1] eq "1") ) { $append = 1; };  next SWITCH; }
 }




########################## Determine forests to copy. ##########################################


    my $forestFile	= "forests.list";
    my $locationFile	= "locations.dat";
    my $startLine   	= 1;
    my %forestHash;
    my %ignoreHash;
    my $treeCount  	= 0;
    my $forestCount	= 0;
    my $timestamp;
    my @timestamp;
    my $urlHeader;
    my $url;
    my $getFile;
    my $fileGZ;
    my $file;

    # Make hash of forests to ignore.
    $file = $ignoreFile; 
    if ( -f "$workDir/$file" )
    {
	open( IFORESTS, "<", "$workDir/$file" ) or die ("Could not open $workDir/$file : $!\n");
	while( my $ignoreLine = <IFORESTS> )
	{
	    $ignoreHash{ $ignoreLine } = 1;
	}
	close( IFORESTS ) or die "Could not close $workDir/$file: $!\n";
    }

    # If necessary download 'forests.list'.
    $file = $forestFile;
    $fileGZ = "$file.gz";
    $url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/$fileGZ";
    if ( (-f "$workDir/$file") && (-s "$workDir/$file") ) {
	@timestamp = localtime((stat "$workDir/$file")[9]);
    } else {
	@timestamp = localtime(0);
    }
    $timestamp = POSIX::strftime("%a, %d %b %Y %H:%M:%S %Z", @timestamp);
    $urlHeader = "If-Modified-Since: $timestamp";
    $getFile = system( "wget", "-e robots=off", "--quiet",  "--header=$urlHeader", "--output-document=$workDir/$fileGZ", "$url" );
    if ( $getFile == 0 ) { 
	system( "gunzip", "--force", "$workDir/$fileGZ") == 0  or die "Could not unzip $workDir/$fileGZ: $!\n"; 
    } else {
	system( "rm", "--force", "$workDir/$fileGZ" );
    }

    # If necessary download locations.dat.
    $file = $locationFile;
    $fileGZ = "$file.gz";
    $url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/$fileGZ";
    if ( (-f "$workDir/$file") && (-s "$workDir/$file") ) {
        @timestamp = localtime((stat "$workDir/$file")[9]);
    } else {
        @timestamp = localtime(0);
    }
    $timestamp = POSIX::strftime("%a, %d %b %Y %H:%M:%S %Z", @timestamp);
    $urlHeader = "If-Modified-Since: $timestamp";
    $getFile = system( "wget", "-e robots=off", "--quiet", "--header=$urlHeader", "--output-document=$workDir/$fileGZ", "$url" );
    if ( $getFile == 0 ) { 
        system( "gunzip", "--force", "$workDir/$fileGZ") == 0  or die "Could not unzip $workDir/$fileGZ: $!\n"; 
    } else {
        system( "rm", "--force", "$workDir/$fileGZ" );
    }

    # Collect forests.
    open( FORESTS,  "<", "$workDir/$forestFile" ) or die "Could not open $workDir/$forestFile: $!\n";

    # Wrapper loop to find more than one forest.
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
	open( IFOREST, ">>", "$workDir/$ignoreFile" ) or die "Appending forests to ignoreFile: Could not open $workDir/$ignoreFile: $!\n";
	foreach my $forestID ( keys %forestHash )
	{
	    print IFOREST "$forestID\n";
	}
	close ( IFOREST ) or die "Appending forests to ignoreFile: Could not close $workDir/$ignoreFile: $!\n";
    }




########## Copy all found forests to output file in ascending order. ##########


 open( MYTREES, ">", "$workDir/$outFile" ) or die "Could not open $workDir/$outFile: $!\n";
 my $header = 0;

 # Copy forests in ascending order.
 foreach my $forestID ( sort( keys %forestHash ) )
 {
    # Copy trees belonging to $forestID.
    foreach my $treeID ( keys %{ $forestHash{ $forestID } } )
    {
	my $offset   = $forestHash{ $forestID }{ $treeID }[0];
	   $file     = $forestHash{ $forestID }{ $treeID }[1];

	# Download and unzip treefile if necessary.
	$fileGZ = "$file.gz";
	$url = "http://www.slac.stanford.edu/~behroozi/Bolshoi_Trees/$fileGZ";
	if ( (-f "$workDir/$file") && (-s "$workDir/$file") ) {
		@timestamp = localtime((stat "$workDir/$file")[9]);
	} else {
		@timestamp = localtime(0);
	}
	$timestamp = POSIX::strftime("%a, %d %b %Y %H:%M:%S %Z", @timestamp);
	$urlHeader = "If-Modified-Since: $timestamp";
	$getFile = system( "wget", "-e robots=off", "--quiet", "--header=$urlHeader", "--output-document=$workDir/$fileGZ", "$url" );
	if ( $getFile == 0 ) {
		system( "gunzip", "--force", "$workDir/$fileGZ") == 0  or die "Could not unzip $workDir/$fileGZ: $!\n";
	} else {
		system( "rm", "--force", "$workDir/$fileGZ" );
		die "Could not download $file from $url. Programm aborted.\n" unless ( (-f "$workDir/$file") && (-s "$workDir/$file") )
        }

	# Copy tree to file.
	open( NEWTREES, "<", "$workDir/$file" ) or die "Could not open $workDir/$file: $!\n";

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
	close( NEWTREES ) or die "Could not close $workDir/$file: $!\n";
    }
 }

 close( MYTREES ) or die "Could not close $workDir/$outFile: $!\n";
 print STDOUT "Copied $treeCount merger trees of $forestCount forests to file $workDir/$outFile, sorted in ascending order.\n";

 exit;
