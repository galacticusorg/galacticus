#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use XML::Simple;
use File::Find;
use Term::ReadKey;
use Data::Dumper;
use Switch;

# Read in any configuration options.
my $config;
if ( -e "galacticusConfig.xml" ) {
    my $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Identify e-mail options for this host.
my $dbConfig;
if ( exists($ENV{'HOSTNAME'}) && exists($config->{'millenniumDB'}->{'host'}->{$ENV{'HOSTNAME'}}) ) {
    $dbConfig = $config->{'millenniumDB'}->{'host'}->{$ENV{'HOSTNAME'}};
} elsif ( exists($config->{'millenniumDB'}->{'host'}->{'default'}) ) {
    $dbConfig = $config->{'millenniumDB'}->{'host'}->{'default'};
} else {
    print "Please enter your Millennium database username:\n";
    while (1) {
	ReadMode "raw";
	my $c = ReadKey(0);
	ReadMode "normal";
	print $c;
	last if $c eq "\n";
	$dbConfig->{'user'} .= $c;
    }
    $dbConfig->{'passwordFrom'} = "input";
}

# Get any password now.
my $dbPassword;
switch ( $dbConfig->{'passwordFrom'} ) {
    case ( "input" ) {
	print "Please enter your Millennium database password:\n";
	$dbPassword = &getPassword;
    }
    require Net::DBus;
}

# Run the script to grab the trees.
system("scripts/aux/Millennium_Trees_Grab.pl --user ".$dbConfig->{'user'}." --password ".$dbPassword." --select \"root.m_tophat > 1.0e2 and root.m_tophat < 10.0e2\" --output tests/nBodyMergerTrees/Millennium_Trees.csv");

# Process into Galacticus format.
system("make Millennium_Merger_Tree_File_Maker.exe; ./Millennium_Merger_Tree_File_Maker.exe tests/nBodyMergerTrees/Millennium_Trees.csv tests/nBodyMergerTrees/Millennium_Trees_Particles.csv tests/nBodyMergerTrees/Millennium_Trees.hdf5 galacticus 1");

exit;

sub getPassword {
    # Read a password from standard input while echoing asterisks to the screen.
    ReadMode('noecho');
    ReadMode('raw');
    my $password = '';
    while (1) {
	my $c = ReadKey(0);
	last if $c eq "\n";
	print "*";
	$password .= $c;
    }
    ReadMode('restore');
    print "\n";
    return $password;
}
