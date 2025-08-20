#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Galacticus::Build::Directives;
use List::ExtraUtils;

# Monitor for updates to certain files, outputting warnings if so.
# Andrew Benson (09-April-2025)

# Gets list of changed files.
my @filesChanged = @ARGV;

# Get list of changes to watch.
my $xml     = new XML::Simple();
my $watches = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/scripts/aux/watches.xml", KeyAttr => 0);

# Initialize a list of warnings.
my $warnings;

# Iterate over files, extracting directives.
foreach my $fileName ( @filesChanged ) {
    # Iterate over all watches.
    foreach my $watch ( &List::ExtraUtils::as_array($watches->{'watch'}) ) {
	next
	    unless ( exists($watch->{'file'}) );
	if ( $watch->{'file'} eq $fileName ) {
	    $warnings .= ":warning: File `".$fileName."` has changed. ".$watch->{'message'}."\n";
	}
    }
    # Iterate over all directives in this file.
    foreach my $directive ( &Galacticus::Build::Directives::Extract_Directives($fileName,"*",setRootElementType => 1) ) {
	# Iterate over all watches.
	foreach my $watch ( &List::ExtraUtils::as_array($watches->{'watch'}) ) {
	    next
		unless ( exists($watch->{'type'}) );
	    if (
		$directive->{'rootElementType'} eq $watch->{'type'}
		&&
		$directive->{'name'} eq $watch->{'name'}
		) {
		$warnings .= ":warning: File `".$fileName."` has changed. ".$watch->{'message'}."\n";
	    }
	}
    }
}

# Output warnings if any exist.
if ( defined($warnings) ) {
    open(my $warningFile,">",$ENV{'GALACTICUS_EXEC_PATH'}."/warnings.md");
    print $warningFile "# Watched files have changed\n";
    print $warningFile $warnings;
    close($warningFile);
}

exit;
