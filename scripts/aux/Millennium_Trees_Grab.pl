#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;

# Grab tree and particle data from the Millennium databases.
# Andrew Benson (18-Mar-2010)

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

# Initialize database access username and password.
my $sqlUser;
my $sqlPassword;

# Specify user and password.
$sqlUser     = $arguments{"user"    }
    if ( exists($arguments{'user'    }) );
$sqlPassword = $arguments{"password"}
    if ( exists($arguments{'password'}) );

# Parse the Galacticus config file if it is present.
my $config = &Galacticus::Options::LoadConfig();
if ( defined($config) ) {
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
my $selection   = "";
$selection = $arguments{"select"}
   if ( exists($arguments{"select"}) );
if ( exists($arguments{"selectFile"}) ) {
    open(iHndl,$arguments{"selectFile"});
    $selection = <iHndl>;
    chomp($selection);
    close(iHndl);
}

# Specify the treeId property.
my $treeId;
if ( exists($arguments{"treeId"}) ) {
    $treeId = $arguments{"treeId"};
} else {
    $treeId = "treeId";
}

# Specify the haloId property.
my $haloId;
if ( exists($arguments{"haloId"}) ) {
    $haloId = $arguments{"haloId"};
} else {
    $haloId = "haloId";
}

# Specify the descendantId property.
my $descendantId;
if ( exists($arguments{"descendantId"}) ) {
    $descendantId = $arguments{"descendantId"};
} else {
    $descendantId = "descendantId";
}

# Specify the velocity dispersion property.
my $velocityDispersion;
if ( exists($arguments{"velocityDispersion"}) ) {
    $velocityDispersion = $arguments{"velocityDispersion"};
} else {
    $velocityDispersion = "velDisp";
}

# Specify the output file.
my $outputFile;
if ( exists($arguments{"output"}) ) {
    $outputFile = $arguments{"output"};
} else {
    $outputFile = "Millennium_Trees.csv";
}

# Specify the database table.
my $table;
if ( exists($arguments{"table"}) ) {
    $table = $arguments{"table"};
} else {
    $table = "millimil..MPAHalo";
}

# Specify the index table.
my $indexTable;
my $indexNode;
if ( exists($arguments{"indexTable"}) ) {
    $indexNode  = "indexNode";
    $indexTable = $arguments{"indexTable"};
} else {
    $indexNode = "node";
    $indexTable = $table;
}
$indexNode = "node"
    if ( $indexTable eq $table );

# Specify the snapshot table.
my $snapshotTable;
if ( exists($arguments{"snapshotTable"}) ) {
    $snapshotTable = $arguments{"snapshotTable"};
} else {
    $snapshotTable = "millimil..Snapshots";
}

# Specify the particle table.
my $particleTable;
if ( exists($arguments{"particleTable"}) ) {
    $particleTable = $arguments{"particleTable"};
} else {
    $particleTable = "MMSnapshots..MillimilSnapshots";
}

# Determine if particles should be traced.
my $traceParticles;
if ( exists($arguments{"traceParticles"}) ) {
    $traceParticles = $arguments{"traceParticles"};
} else {
    $traceParticles = "yes";
}

# Determine mass to use.
my $mass;
if ( exists($arguments{"mass"}) ) {
    $mass = $arguments{"mass"};
} else {
    $mass = "m_tophat";
}

# Determine half-mass radius to use.
my $halfMassRadius = "halfmassRadius";
$halfMassRadius = $arguments{"halfMassRadius"}
    if ( exists($arguments{"halfMassRadius"}) );

# Specify the database URL.
my $databaseURL = "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=";

# Build the retrieve command base.
my $getCommandBase = "wget";
$getCommandBase   .= " --no-proxy --retry-connrefused";
$getCommandBase   .= " --http-user="  .$sqlUser     unless ( $sqlUser     eq "" );
$getCommandBase   .= " --http-passwd=".$sqlPassword unless ( $sqlPassword eq "" );

# Build the SQL query to retrieve basic node data. We convert quantities containing a physical length into comoving units.
my $plus = "%2b";
my $sqlQuery = $databaseURL."select ".$indexNode.".".$treeId.", ".$indexNode.".".$haloId.", ".$indexNode.".".$descendantId.", node.firstHaloInFOFgroupId, node.snapNum, node.redshift, node.".$mass.", node.np, node.x, node.y, node.z, node.velX, node.velY, node.velZ, node.spinX*(1.0".$plus."node.redshift), node.spinY*(1.0".$plus."node.redshift), node.spinZ*(1.0".$plus."node.redshift), node.".$halfMassRadius."*(1.0".$plus."node.redshift), node.mostBoundID, node.vMax, node.".$velocityDispersion." from ".$table." node";
# Add the root node if it is used in the selection.
$sqlQuery .= ", ".$table." root"
    if ( $selection =~ m/root\./ );
# Add the index node if the index table differs from the node table.
$sqlQuery .= ", ".$indexTable." ".$indexNode.""
    if ( $indexTable ne $table );
# Create a where statement.
my $where = " where";
# Add implicit join on root node only if the selection uses the root node.
$where .= " root.haloId = node.treeId"
    if ( $selection =~ m/root\./ );
# Add implicit join on index node only if the index table differs from the node table.
if ( $indexTable ne $table ) {
    $where .= " and"
	unless ( $where eq " where" );
    $where .= " ".$indexNode.".".$haloId." = node.haloId";
}
# Append any required selection.
unless ( $selection eq "" ) {
    $where .= " and "
	unless ( $where eq " where" );
    $where .= $selection;
}
# Append the where statement.
$sqlQuery .= $where
    unless ( $where eq " where" );
# Add an order by statement.
$sqlQuery .= " order by ".$indexNode.".".$treeId;
# Retrieve the data.
my $getCommand = $getCommandBase." \"".$sqlQuery."\" -O ".$outputFile;
system($getCommand);

# Trace particles if requested.
if ( $traceParticles eq "yes" ) {
    # Build the SQL query to retrieve particle data for lost subhalos.
    my $sqlQuery = $databaseURL."select ".$indexNode.".".$descendantId.", count(*) as num into countTable from ".$indexTable." ".$indexNode.", ".$table." root, ".$indexTable." indexNode where root.haloId = node.treeId and node.haloId = indexNode.".$haloId;
    $sqlQuery   .= ", ".$indexTable." indexNode"
    unless ( $indexNode eq "node" );
    $sqlQuery   .= " and node.haloId = indexNode.".$haloId
	unless ( $indexNode eq "node" );
    # Append any required selection.
    $sqlQuery   .= " and ".$selection unless ( $selection eq "" );
    # Append grouping command.
    $sqlQuery .= " group by ".$indexNode.".".$descendantId;
    # Add command to find nodes which are not the only antescendents.
    $sqlQuery   .= "; select node.mostBoundId, node.snapNum into boundTable from countTable, ".$table." node where countTable.num > 1 and node.descendantId = countTable.".$descendantId;
    # Add command to select most bound particles trajectories corresponding to these nodes.
    $sqlQuery   .= "; select snap.id, times.redshift, snap.snapNum, snap.x, snap.y, snap.z, snap.vx, snap.vy, snap.vz from boundTable, ".$particleTable." snap, ".$snapshotTable." times where boundTable.mostBoundId = snap.id and boundTable.snapNum <= snap.snapnum and times.snapnum = snap.snapNum";
    # Add commands to drop the temporary tables;
    $sqlQuery   .= "; drop table countTable; drop table boundTable";
    # Retrieve the data.
    my $getCommand = $getCommandBase." \"".$sqlQuery."\" -O particles.tmp";
    system($getCommand);
    # Sort the output to get ordering by particle number and then by (decreasing) redshift.
    $outputFile =~ s/\.csv/_Particles.csv/;
    system("sort -g -k1,1 -k2,2r -t, -u particles.tmp > ".$outputFile);
    unlink("particles.tmp");
}

exit;
