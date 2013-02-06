#!/usr/bin/env perl
use lib "./perl";
use XML::Simple;
use File::Find;
use Term::ReadKey;
use Net::DBus;
use Switch;

# Read in any configuration options.
if ( -e "galacticusConfig.xml" ) {
    $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Identify e-mail options for this host.
if ( exists($config->{'millenniumDB'}->{'host'}->{$ENV{'HOST'}}) ) {
    $dbConfig = $config->{'millenniumDB'}->{'host'}->{$ENV{'HOST'}};
} elsif ( exists($config->{'millenniumDB'}->{'host'}->{'default'}) ) {
    $dbConfig = $config->{'millenniumDB'}->{'host'}->{'default'};
} else {
    print "Please enter your Millennium database username:\n";
    while (1) {
	$c = ReadKey(-1);
	last if $c eq "\n";
	$dbConfig->{'user'} .= $c;
    }
    $dbConfig->{'method'} = "input";
}

# Get any password now.
switch ( $dbConfig->{'passwordFrom'} ) {
    case ( "input" ) {
	print "Please enter your Millennium database password:\n";
	$dbPassword = &getPassword;
    }
    case ( "kdewallet" ) {
	$appName          = "Galacticus";
	$folderName       = "glc-millennium-db";
	my $bus           = Net::DBus->find;
	my $walletService = $bus->get_service("org.kde.kwalletd");
	my $walletObject  = $walletService->get_object("/modules/kwalletd");
	my $walletID      = $walletObject->open("kdewallet",0,$appName);
	if ( $walletObject->hasEntry($walletID,$folderName,"dbPassword",$appName) == 1 ) {
	    $dbPassword = $walletObject->readPassword($walletID,$folderName,"dbPassword",$appName); 
	} else {
	    print "Please enter your Millennium database password:\n";
	    $dbPassword = &getPassword;
	    $walletObject->writePassword($walletID,$folderName,"dbPassword",$dbPassword,$appName); 
	}
    }
}

# Run the script to grab the trees.
system("scripts/aux/Millennium_Trees_Grab.pl --user ".$dbConfig->{'user'}." --password ".$dbPassword." --select \"root.m_tophat > 1.0e2 and root.m_tophat < 10.0e2\" --output tests/nBodyMergerTrees/Millennium_Trees.csv");

# Process into Galacticus format.
system("make Millennium_Merger_Tree_File_Maker.exe; Millennium_Merger_Tree_File_Maker.exe tests/nBodyMergerTrees/Millennium_Trees.csv tests/nBodyMergerTrees/Millennium_Trees_Particles.csv tests/nBodyMergerTrees/Millennium_Trees.hdf5");

exit;

sub getPassword {
    # Read a password from standard input while echoing asterisks to the screen.
    ReadMode('noecho');
    ReadMode('raw');
    my $password = '';
    while (1) {
	my $c;
	1 until defined($c = ReadKey(-1));
	last if $c eq "\n";
	print "*";
	$password .= $c;
    }
    ReadMode('restore');
    print "\n";
    return $password;
}
