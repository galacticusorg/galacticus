#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use Data::Dumper;

# Grab treeand particle data from the Millennium databases.
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

# Specify user and password.
my $sqlUser     = $arguments{"user"};
my $sqlPassword = $arguments{"password"};

# Specify any selection.
my $selection   = $arguments{"select"};

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
if ( exists($arguments{"indexTable"}) ) {
    $indexTable = $arguments{"indexTable"};
} else {
    $indexTable = $table;
}

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

# Specify the database URL.
my $databaseURL = "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=";

# Build the retrieve command base.
my $getCommandBase = "wget";
$getCommandBase   .= " --http-user="  .$sqlUser     unless ( $sqlUser     eq "" );
$getCommandBase   .= " --http-passwd=".$sqlPassword unless ( $sqlPassword eq "" );

# Build the SQL query to retrieve basic node data.
my $sqlQuery = $databaseURL."select indexNode.".$treeId.", indexNode.".$haloId.", indexNode.".$descendantId.", node.firstHaloInFOFgroupId, node.snapNum, node.redshift, node.".$mass.", node.np, node.x, node.y, node.z, node.velX, node.velY, node.velZ, node.spinX, node.spinY, node.spinZ, node.halfmassRadius, node.mostBoundID from ".$table." node, ".$table." root, ".$indexTable." indexNode where root.haloId = node.treeId and node.haloId = indexNode.".$haloId;
# Append any required selection.
$sqlQuery .= " and ".$selection unless ( $selection eq "" );
# Add an order by statement.
$sqlQuery .= " order by indexNode.".$treeId;
# Retrieve the data.
my $getCommand = $getCommandBase." \"".$sqlQuery."\" -O ".$outputFile;
system($getCommand);

# Trace particles if requested.
if ( $traceParticles eq "yes" ) {
    # Build the SQL query to retrieve particle data for lost subhalos.
    my $sqlQuery = $databaseURL."select indexNode.".$descendantId.", count(*) as num into countTable from ".$table." node, ".$table." root, ".$indexTable." indexNode where root.haloId = node.treeId and node.haloId = indexNode.".$haloId;
    # Append any required selection.
    $sqlQuery .= " and ".$selection unless ( $selection eq "" );
    # Append grouping command.
    $sqlQuery .= " group by indexNode.".$descendantId;
    # Add command to find nodes which are not the only antescendents.
    $sqlQuery .= "; select node.mostBoundId, node.snapNum into boundTable from countTable, ".$table." node where countTable.num > 1 and node.descendantId = countTable.".$descendantId;
    # Add command to select most bound particles trajectories corresponding to these nodes.
    $sqlQuery .= "; select snap.id, times.redshift, snap.snapNum, snap.x, snap.y, snap.z, snap.vx, snap.vy, snap.vz from boundTable, ".$particleTable." snap, ".$snapshotTable." times where boundTable.mostBoundId = snap.id and boundTable.snapNum <= snap.snapnum and times.snapnum = snap.snapNum";
    # Add commands to drop the temporary tables;
    $sqlQuery .= "; drop table countTable; drop table boundTable";
    # Retrieve the data.
    my $getCommand = $getCommandBase." \"".$sqlQuery."\" -O particles.tmp";
    system($getCommand);
    # Sort the output to get ordering by particle number and then by (decreasing) redshift.
    my $outputFile =~ s/\.csv/_Particles.csv/;
    system("sort -g -k1,1 -k2,2r -t, -u particles.tmp > ".$outputFile);
    unlink("particles.tmp");
}

exit;
