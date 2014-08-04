#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use XML::SAX;
use XML::SAX::Writer;
use XML::Filter::XInclude;
use XML::Twig;

# Parse an XML file, include any files via XInclude elements, and write out the results.
# Andrew Benson (8-January-2014)

die("Usage: xmlInclude.pl <inFile> <outFile>")
    unless ( scalar(@ARGV) == 2 );
my $inFile  = $ARGV[0];
my $outFile = $ARGV[1];

# Define output buffer.
my $content;
# Construct an XML writer.
my $writer = XML::SAX::Writer->new(Output => \$content);
# Construct an XInclude filter.
my $filter = XML::Filter::XInclude->new(Handler => $writer);
# Construct an XML parser.
my $parser = XML::SAX::ParserFactory->parser(Handler => $filter);
# Parse the file/
$parser->parse_uri($inFile);
# Strip comments.
my $twig = XML::Twig->new (comments => 'drop', pretty_print => 'indented');
$twig->parse($content);
# Dump the output to file.
open(my $outHndl,">".$outFile);
print $outHndl $twig->sprint();
close($outHndl);

exit;
