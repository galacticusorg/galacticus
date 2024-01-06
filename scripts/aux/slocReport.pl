#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use JSON;

# Count source lines of code in Galacticus source files, accounting for embedded XML and LaTeX.
# Andrew Benson (04-January-2024)

# Initialize counts.
my $count;

# Get a list of all source files to process.
my @sourceFiles;
find(sub {push(@sourceFiles,$_ =~ m/\.F90$/ || $_ =~ m/\.Inc/ ? $File::Find::name : ()) },"source");
my @docFiles;
find(sub {push(@docFiles   ,$_ =~ m/\.tex$/                   ? $File::Find::name : ()) },"doc"   );

# Iterate over source files.
foreach my $fileName ( @sourceFiles ) {
    # Initialize state of embedded code.
    my $inXML        = 0;
    my $inLaTeX      = 0;
    # Read the file.
    open(my $file,$fileName);	
    while ( my $line = <$file> ) {
	# Detect entry and exit of embedded regions.
	if      ( $line =~ m/^\s*!!\[\s*$/ ) {
	    $inXML   = 1;
	} elsif ( $line =~ m/^\s*!!\]\s*$/ ) {
	    $inXML   = 0;
	} elsif ( $line =~ m/^\s*!!\{\s*$/ ) {
	    $inLaTeX = 1;
	} elsif ( $line =~ m/^\s*!!\}\s*$/ ) {
	    $inLaTeX = 0;
	}
	# Skip empty lines, continuations of lines, and comments (unless they are an OpenMP directive).
	next
	    if (
		$line =~ m/^\s*$/
		||
		$line =~ m/^\s*&/
		||
		(
		 $line =~ m/^\s*!/
		 &&
		 $line !~ m/^\s*!\$/
		)
	    );
	# Count the line.
	if      ( $inXML   ) {
	    ++$count->{'xml'    };
	} elsif ( $inLaTeX ) {
	    ++$count->{'latex'  };
	} else {
	    ++$count->{'fortran'};
	}
    }
    close($file);
}

# Iterate over documentation files.
foreach my $fileName ( @docFiles ) {
    # Read the file.
    open(my $file,$fileName);	
    while ( my $line = <$file> ) {
	# Skip empty lines, and comments.
	next
	    if (
		$line =~ m/^\s*$/
		||
		$line =~ m/^\s*%/
	    );
	# Count the line.
	++$count->{'latex'};
    }
    close($file);
}

# Count everything else using sloccount.
system("sloccount aux constraints parameters parameters.xml perl plots schema scripts testSuite source > sloccount.out");
open(my $sloccount,"sloccount.out");
while ( my $line = <$sloccount> ) {
    if ( $line =~ m/^([a-z0-9]+):\s*(\d+)/ ) {
	my $language   = $1;
	my $countLines = $2;
	$count->{$language} += $countLines
	    unless ( $language eq "f90" );
    }
}
close($sloccount);
unlink("sloccount.out");

# Report the results.
system("curl -X POST -H 'Content-type: application/json' --data '".encode_json($count)."' ".$ENV{'SLACK_WEBHOOK_SLOCREPORT_URL'});

exit 0;
