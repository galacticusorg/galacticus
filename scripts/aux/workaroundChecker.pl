#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;

# Check for closed GCC PRs in workarounds.
# Andrew Benson (17-September-2021)

# Initialize list of workarounds.
our $workarounds;

# Scan for files.
my @directories = ("./source","./perl/Galacticus/Build");
find(\&fileMatcher,@directories);

# Get status of PRs.
&checkLinks();

# Generate and send report.
my $status;
my $resolved = grep {$workarounds->{$_}->{'status'} eq "RESOLVED"} keys(%{$workarounds});
if ( $resolved ) {
    $status = 1;
    print "!!! Resolved PRs with workarounds exist !!!\n\n";
} else {
    $status = 0;
    print "No resolved PRs with workarounds exist\n\n";
}
foreach my $PR ( sort(keys(%{$workarounds})) ) {
    print "!!! "
	if ( $workarounds->{$PR}->{'status'} eq "RESOLVED" );
    print "PR".$PR." (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=".$PR."):\n";
    foreach my $file ( sort(keys(%{$workarounds->{$PR}->{'files'}})) ) {
	print " -> ".$file."\n";
    }
    print "\n\n";
}

exit $status;

sub fileMatcher {
    my $fileName = $_;
    my $fullName = $File::Find::name;
    return
	unless ( $fullName =~ m/\.(F90|Inc|pm)$/ );
    open(my $file,$fullName);
    while ( my $line = <$file> ) {
	if ( $line =~ m/^\s*#??\s*<workaround\s.*PR="(\d+)"/ ) {
	    my $PR = $1;
	    ++$workarounds->{$PR}->{'files'}->{$fullName};
	}
    }
    close($file);
}

sub checkLinks {
    foreach my $PR ( keys(%{$workarounds}) ) {
	$workarounds->{$PR}->{'status'} = "UNKNOWN";
	my $url = "https://gcc.gnu.org/bugzilla/show_bug.cgi?id=".$PR;
	system("sleep 1; curl --insecure --location --output pr.html --fail \"".$url."\"");
	open(my $report,"pr.html");
	while ( my $line = <$report> ) {
	    if ( $line =~ m/<span id="static_bug_status">([A-Z]+)/ ) {
		$workarounds->{$PR}->{'status'} = $1;
	    }
	}
	close($report);
	unlink("pr.html");
    }
}
