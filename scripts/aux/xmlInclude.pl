#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use XML::LibXML;
use XML::Twig;

# Parse an XML file, include any files via XInclude elements, and write out the results.
# Andrew Benson (8-January-2014)

die("Usage: xmlInclude.pl <inFile> <outFile>")
    unless ( scalar(@ARGV) == 2 );
my $inFile  = $ARGV[0];
my $outFile = $ARGV[1];

# Parse the content.
my $parser = XML::LibXML->new();
my $dom    = $parser->load_xml(location => $inFile);
$parser->process_xincludes($dom);
my $twig = XML::Twig->new (comments => 'drop', pretty_print => 'indented', escape_gt => 1);
$twig->parse($dom->serialize());

# Dump the output to file.
open(my $outHndl,">".$outFile);
print $outHndl $twig->sprint();
close($outHndl);

exit;
