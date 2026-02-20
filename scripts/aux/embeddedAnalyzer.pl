#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use XML::Validator::Schema;
use XML::SAX::ParserFactory;
use File::Slurp qw(slurp);
use LaTeX::SpellCheck;

# Perform checks on embedded XML and LaTeX
# Andrew Benson (28-February-2023)

# Get the file to process.
die('Usage: embeddedAnalyzer.pl <fileName> <warningFile>')
    unless ( scalar(@ARGV) == 2 );
my $fileName        = $ARGV[0];
my $warningFileName = $ARGV[1];

# Begin parsing the file.
my $status            = 0              ;
my $inDirective       = 0              ;
my $inXML             = 0              ;
my $inLaTeX           = 0              ;
my $inComment         = 0              ;
my $lineNumber        = 0              ;
my $xml               = new XML::Simple;
my $warnings                           ;
my $directiveRoot                      ;
my $rawCode                            ;
my $rawDirective                       ;
my $rawLaTeX                           ;
my $rawComment                         ;
my $strippedDirective                  ;
open(my $code,$fileName);
while ( my $line = <$code> ) {
    # Detect the end of a LaTeX section and change state.
    $inLaTeX = 0
	if ( $line =~ m/^\s*!!\}/ );
    # Detect the end of an XML section and change state.
    $inXML = 0
	if ( $line =~ m/^\s*!!\]/ );
    # Detect the end of a LaTeX section and change state.
    $inComment = 0
	if ( $line !~ m/^\s*!\s/ );
    # Process LaTeX blocks.
    $rawLaTeX .= $line
	if ( $inLaTeX );
    if ( defined($rawLaTeX) && ! $inLaTeX ) {
	my $LaTeXLog = &testLaTeX($rawLaTeX);
	if ( defined($LaTeXLog) ) {
	    print "LaTeX fragment compilation failed (".$fileName.":".$lineNumber."):\n".$LaTeXLog;
	    $status = 1;
	}
	$warnings .= &LaTeX::SpellCheck::spellCheck($rawLaTeX,"latex",$fileName);
	undef($rawLaTeX);
    }
    # Process comment blocks.
    $rawComment .= $line
	if ( $inComment );
    if ( defined($rawComment) && ! $inComment ) {
	$warnings .= &LaTeX::SpellCheck::spellCheck($rawComment,"text",$fileName);
	undef($rawComment);
    }
    # Process XML blocks.
    my $isDirective  = 0;
    my $endDirective = 0;
    (my $strippedLine = $line) =~ s/^\s*\!<\s*//;
    $strippedLine =~ s/&nbsp;/ /g;
    if ( $inXML ) {
	# Determine if line is a directive line.
	$isDirective    = 1
	    if ( $strippedLine =~ m/^\s*\<([^\s\>\/]+)/ || $inDirective == 1 );
	$directiveRoot = $1
	    if ( $isDirective == 1 && $inDirective == 0 );		
	# Catch the end of directives.
	$endDirective = 1
	    if ( $isDirective == 1 && $strippedLine =~ m/\s*\<\/$directiveRoot\>/ );
	$endDirective = 1
	    if ( $isDirective == 1 && $inDirective == 0 && ( $strippedLine =~ m/\s*\<$directiveRoot\s.*\/\>/ || $strippedLine =~ m/\s*\<$directiveRoot\/\>/ ) );
	# Record whether we are currently in or out of a directive.
	$inDirective = 1
	    if ( $isDirective == 1 );
    }
    # Accumulate raw text.
    if ( $inDirective ) {
	$rawDirective      .= $line;
	$strippedDirective .= $strippedLine;
    } elsif ( $line !~ m/^\s*!!(\[|\])/ ) {
	$rawCode           .= $line;
    }
    # Process code and directive blocks as necessary.
    if ( ( $inDirective == 0 || eof($code) || $endDirective ) && $rawDirective ) {
	# Attempt to parse the directive XML.
	my $directive = eval{$xml->XMLin($strippedDirective, keepRoot => 1)};
	if ( $@ ) {
	    print "Parsing XML fragment failed (".$fileName.":".$lineNumber.")\n".$@."\n";
	    print $strippedDirective;
	    $status = 1;
	}
	my $directiveName = (keys %{$directive})[0];
	# Validate the directive if possible.
	if ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/schema/".$directiveName.".xsd" ) {
	    my $validator = XML::Validator::Schema->new(file => $ENV{'GALACTICUS_EXEC_PATH'}."/schema/".$directiveName.".xsd");
	    my $parser    = XML::SAX::ParserFactory->parser(Handler => $validator); 
	    eval { $parser->parse_string($strippedDirective) };
	    if ( $@ ) {
		print "XML fragment validation failed (".$fileName.":".$lineNumber."):\n".$@."\n";
		$status = 1;
	    }
	}
	# Look for a "description" element in the directive - this should be LaTeXable.
	if ( exists($directive->{$directiveName}->{'description'}) ) {
	    my $LaTeXLog = &testLaTeX($directive->{$directiveName}->{'description'});
	    if ( defined($LaTeXLog) ) {
		print "XML LaTeX description compilation failed (".$fileName.":".$lineNumber."):\n".$LaTeXLog;
		$status = 1;
	    }
	    $warnings .= &LaTeX::SpellCheck::spellCheck($directive->{$directiveName}->{'description'},"latex",$fileName);
	}
	# Reset the raw directive text.
	$inDirective = 0;
	undef($rawDirective     );
	undef($strippedDirective);
    }
    # Detect the start of an XML section and change state.
    $inXML = 1
	if ( $line =~ m/^\s*!!\[/ );
    # Detect the start of a LaTeX section and change state.
    $inLaTeX = 1
	if ( $line =~ m/^\s*!!\{/ );
    # Detect the start of a comment section and change state.
    if ( $line =~ m/^\s*!\s/ ) {
	$inComment   = 1;
	$rawComment .= $line;
    }
    # Increment line number count.
    ++$lineNumber;
}
close($code);

# Write accumulated warnings to file (appending).
unless ( $warnings eq "" ) {
    open(my $warningFile,">>",$warningFileName);
    print $warningFile $warnings;
    close($warningFile);
}

exit $status;

sub testLaTeX {
    # Test that a LaTeX fragment compiles.
    my $rawLaTeX = shift();
    $rawLaTeX =~ s/&amp;/&/g;
    $rawLaTeX =~ s/&lt;/>/g;
    open(my $LaTeX,">doc/frag.tex");
    print $LaTeX "\\documentclass[letterpaper,10pt,headsepline]{scrbook}\n";
    print $LaTeX "\\usepackage{natbib}\n";
    print $LaTeX "\\usepackage{epsfig}\n";
    print $LaTeX "\\usepackage[acronym]{glossaries}\n";
    print $LaTeX "\\usepackage[backref,colorlinks]{hyperref}\n";
    print $LaTeX "\\usepackage{amssymb}\n";
    print $LaTeX "\\usepackage{amsmath}\n";
    print $LaTeX "\\usepackage{color}\n";
    print $LaTeX "\\usepackage{listings}\n";
    print $LaTeX "\\usepackage{tensor}\n";
    print $LaTeX "\\input{".$ENV{'GALACTICUS_EXEC_PATH'}."/doc/commands}\n";
    print $LaTeX "\\input{".$ENV{'GALACTICUS_EXEC_PATH'}."/doc/Glossary}\n";
    print $LaTeX "\\newcommand{\\docname}{tmp}\n";
    print $LaTeX "\\begin{document}\n";
    print $LaTeX $rawLaTeX;
    print $LaTeX "\\end{document}\n";
    close($LaTeX);
    system("cd doc; pdflatex -halt-on-error frag > frag.tmp");
    my $log = $? ? slurp("doc/frag.log") : undef();
    foreach my $fileName ( "frag.tex", "frag.pdf", "frag.log", "frag.aux", "frag.tmp", "frag.glo" ) {
	unlink("doc/".$fileName)
	    if ( -e "doc/".$fileName );
    }
    return $log;
}
