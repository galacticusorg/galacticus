#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;

# Grab treeand particle data from the Millennium databases.
# Andrew Benson (18-Mar-2010)

# Create a hash of named arguments.
$iArg = -1;
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Initialize database access username and password.
my $sqlUser;
my $sqlPassword;

# Specify user and password.
$sqlUser     = $arguments{"user"    }
    if ( exists($arguments{'user'    }) );
$sqlPassword = $arguments{"password"}
    if ( exists($arguments{'password'}) );

# Parse the Galacticus config file if it is present.
if ( -e $galacticusPath."/galacticusConfig.xml" ) {
    my $xml    = new XML::Simple;
    my $config = $xml->XMLin($galacticusPath."/galacticusConfig.xml");
    if ( exists($config->{'millenniumDB'}->{'host'}) ) {
	foreach ( keys(%{$config->{'millenniumDB'}->{'host'}}) ) {
	    if ( $_ eq $ENV{'HOSTNAME'} || $_ eq "default" ) {
		$sqlUser     = $config->{'millenniumDB'}->{'host'}->{$_}->{'user'    }
		    if ( exists($config->{'millenniumDB'}->{'host'}->{$_}->{'user'    }) );
		$sqlPassword = $config->{'millenniumDB'}->{'host'}->{$_}->{'password'}
		    if ( exists($config->{'millenniumDB'}->{'host'}->{$_}->{'password'}) );
		if ( exists($config->{'millenniumDB'}->{'host'}->{$_}->{'passwordFrom'}) ) {
		    if ( $config->{'millenniumDB'}->{'host'}->{$_}->{'passwordFrom'} eq "input" ) {
			$sqlPassword = <>;
			chomp($sqlPassword);
		    }
		}
	    }
	}
    }
}

# If no user name or password is supplied, exit.
die("Millennium_Trees_Grab.pl: SQL database access username and password must be supplied")
    unless ( defined($sqlUser) && defined($sqlPassword) );

# Specify any selection.
$selection   = $arguments{"select"};

# Specify the treeId property.
if ( exists($arguments{"treeId"}) ) {
    $treeId = $arguments{"treeId"};
} else {
    $treeId = "treeId";
}

# Specify the haloId property.
if ( exists($arguments{"haloId"}) ) {
    $haloId = $arguments{"haloId"};
} else {
    $haloId = "haloId";
}

# Specify the descendantId property.
if ( exists($arguments{"descendantId"}) ) {
    $descendantId = $arguments{"descendantId"};
} else {
    $descendantId = "descendantId";
}

# Specify the output file.
if ( exists($arguments{"output"}) ) {
    $outputFile = $arguments{"output"};
} else {
    $outputFile = "Millennium_Trees.csv";
}

# Specify the database table.
if ( exists($arguments{"table"}) ) {
    $table = $arguments{"table"};
} else {
    $table = "millimil..MPAHalo";
}

# Specify the index table.
if ( exists($arguments{"indexTable"}) ) {
    $indexNode  = "indexNode";
    $indexTable = $arguments{"indexTable"};
} else {
    $indexNode = "node";
    $indexTable = $table;
}

# Specify the snapshot table.
if ( exists($arguments{"snapshotTable"}) ) {
    $snapshotTable = $arguments{"snapshotTable"};
} else {
    $snapshotTable = "millimil..Snapshots";
}

# Specify the particle table.
if ( exists($arguments{"particleTable"}) ) {
    $particleTable = $arguments{"particleTable"};
} else {
    $particleTable = "MMSnapshots..MillimilSnapshots";
}

# Determine if particles should be traced.
if ( exists($arguments{"traceParticles"}) ) {
    $traceParticles = $arguments{"traceParticles"};
} else {
    $traceParticles = "yes";
}

# Determine mass to use.
if ( exists($arguments{"mass"}) ) {
    $mass = $arguments{"mass"};
} else {
    $mass = "m_tophat";
}

# Specify the database URL.
$databaseURL = "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=";

# Build the retrieve command base.
$getCommandBase = "wget";
$getCommandBase .= " --http-user="  .$sqlUser     unless ( $sqlUser     eq "" );
$getCommandBase .= " --http-passwd=".$sqlPassword unless ( $sqlPassword eq "" );

# Build the SQL query to retrieve basic node data.
$sqlQuery = $databaseURL."select ".$indexNode.".".$treeId.", ".$indexNode.".".$haloId.", ".$indexNode.".".$descendantId.", node.firstHaloInFOFgroupId, node.snapNum, node.redshift, node.".$mass.", node.np, node.x, node.y, node.z, node.velX, node.velY, node.velZ, node.spinX, node.spinY, node.spinZ, node.halfmassRadius, node.mostBoundID from ".$table." node, ".$table." root";
$sqlQuery .= ", ".$indexTable." indexNode"
    unless ( $indexNode eq "node" );
$sqlQuery .= " where node.haloId between root.haloId and root.haloId\%2B999999";
$sqlQuert .= " and node.haloId = indexNode.".$haloId
    unless ( $indexNode eq "node" );

# Append any required selection.
$sqlQuery .= " and ".$selection unless ( $selection eq "" );
# Add an order by statement.
$sqlQuery .= " order by ".$indexNode.".".$treeId;
# Retrieve the data.
$getCommand = $getCommandBase." \"".$sqlQuery."\" -O ".$outputFile;
system($getCommand);

# Trace particles if requested.
if ( $traceParticles eq "yes" ) {
    # Build the SQL query to retrieve particle data for lost subhalos.
    $sqlQuery = $databaseURL."select ".$indexNode.".".$descendantId.", count(*) as num into countTable from ".$table." node, ".$table." root";
    $sqlQuery .= ", ".$indexTable." indexNode"
    unless ( $indexNode eq "node" );
    $sqlQuery .= " where root.haloId = node.treeId";
    $sqlQuery .= " and node.haloId = indexNode.".$haloId
	unless ( $indexNode eq "node" );
    # Append any required selection.
    $sqlQuery .= " and ".$selection unless ( $selection eq "" );
    # Append grouping command.
    $sqlQuery .= " group by ".$indexNode.".".$descendantId;
    # Add command to find nodes which are not the only antescendents.
    $sqlQuery .= "; select node.mostBoundId, node.snapNum into boundTable from countTable, ".$table." node where countTable.num > 1 and node.descendantId = countTable.".$descendantId;
    # Add command to select most bound particles trajectories corresponding to these nodes.
    $sqlQuery .= "; select snap.id, times.redshift, snap.snapNum, snap.x, snap.y, snap.z, snap.vx, snap.vy, snap.vz from boundTable, ".$particleTable." snap, ".$snapshotTable." times where boundTable.mostBoundId = snap.id and boundTable.snapNum <= snap.snapnum and times.snapnum = snap.snapNum";
    # Add commands to drop the temporary tables;
    $sqlQuery .= "; drop table countTable; drop table boundTable";
    # Retrieve the data.
    $getCommand = $getCommandBase." \"".$sqlQuery."\" -O particles.tmp";
    system($getCommand);
    # Sort the output to get ordering by particle number and then by (decreasing) redshift.
    $outputFile =~ s/\.csv/_Particles.csv/;
    system("sort -g -k1,1 -k2,2r -t, -u particles.tmp > ".$outputFile);
    unlink("particles.tmp");
}

exit;
