#!/usr/bin/env perl
use strict;
use warnings;
use JSON;
use File::Find;

# Retrieve and archive copies of all run-time downloaded data files for Galacticus.
# Andrew Benson (19-April-2024)

# Get arguments.
die("Usage: archive.pl <galacticusPath> <archivePath> <slackToken>")
    unless ( scalar(@ARGV) == 3 );
my $galacticusPath = $ARGV[0];
my $archivePath    = $ARGV[1];
my $slackToken     = $ARGV[2];

# Parse the dependencies file.
my $dependencies;
open(my $dependenciesFile,$galacticusPath."/aux/dependencies.yml");
while ( my $line = <$dependenciesFile> ) {
    if ( $line =~ m/^(.*):\s+([\d\.]+)/ ) {
	my $code    = $1;
	my $version = $2;
	my @versions = split(/\./,$version);
	$dependencies->{$code} =
	{
	    version      => $version    ,
	    versionMajor => $versions[0]
	};
    }
}
close($dependenciesFile);

# Scan the source directory for files
my @directories = ($galacticusPath."/source",$galacticusPath."/scripts");
my @links;
find(\&linkFinder,@directories);
$_ = $galacticusPath."/Makefile"; &linkFinder();

# Retrieve links.
my $report;
foreach my $link ( @links ) {
    if ( $link =~ m/^https??:\/\/(.+)\/(.+)/ ) {
	my $path = $1;
	my $file = $2;
	system("mkdir -p ".$archivePath."/".$path);
	my $fileName = $archivePath."/".$path."/".$file;
	unless ( -e $fileName ) {
	    $report->{'report'} .= "RETRIEVING: ".$link."\n";
	    system("wget ".$link." -O ".$fileName);
	    unless ( $? == 0 ) {
		$report->{'report'} .= "\tFAILED: ".$link."\n";
	    }
	} else {
	    $report->{'report'} .= "SKIPPING: (already archived) ".$link."\n";
	}
    }
}

# Report the results.
system("curl -X POST -H 'Content-type: application/json' --data '".encode_json($report)."' https://hooks.slack.com/triggers/".$slackToken);

exit;

sub linkFinder {
    # Find links that may be downloaded at run-time.
    my $fileName = $_;
    return
        unless ( $fileName =~ m/\.(F90|Inc|py)$/ || $fileName =~ m/\/Makefile$/ );
    open(my $file,$fileName);
    while ( my $line = <$file> ) {
	# Fortran source.
	if ( $fileName =~ m/\.(F90|Inc)$/ ) {
	    if ( $line =~ m/^\s*call\s+download\s*\(\s*(["'][^,]+)/ ) {
		my $link = $1;
		$link =~ s/["']//g;
		# Replace dependencies with the actual version number here. Handle the "cloudyVersion" case as a special instance as
		# we have to insert a "c" prefix.
		$link =~ s/\/\/char\(([a-zA-Z]+)VersionMajor\)\/\//$dependencies->{$1}->{'versionMajor'}/g;
		$link =~ s/\/\/char\(cloudyVersion\)\/\//c$dependencies->{'cloudy'}->{'version'}/g;
		$link =~ s/\/\/char\(([a-zA-Z]+)Version\)\/\//$dependencies->{$1}->{'version'}/g;
		# Add the link to the list to retrieve. We skip the "backup" ("old") Cloudy path here.
		push(@links,$link)
		    unless ( $link =~ m/cloudy_releases\/c\d+\/old\// );
	    }
	}
	# Python scripts.
	if ( $fileName =~ m/\.py$/ ) {
	    if ( $line =~ m/^\s*urllib\.request\.urlretrieve\s*\(\s*(["'][^,]+)/ ) {
		my $link = $1;
		$link =~ s/["']//g;
		# Add the link to the list to retrieve.
		push(@links,$link);
	    }
	}
	# Makefiles
	if ( $fileName =~ m/^Makefile*/ ) {
	    if ( $line =~ m/^\s*wget\s+(\-\-\S+\s+)*(\S+)/ ) {
		my $link = $2;
		# Add the link to the list to retrieve.
		push(@links,$link);
	    }
	}
    }
    close($file);
}
