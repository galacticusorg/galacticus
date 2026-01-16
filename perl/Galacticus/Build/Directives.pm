# Contains a Perl module which implements processing of directives.

package Galacticus::Build::Directives;
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;
use File::Names;

sub Extract_Directive {
    # Extract a named directive from a given file handle.
    my $fileHandle    = shift();
    my $directiveName = shift();
    my $state         = shift();
    (my %options)     = @_
	if ( scalar(@_) > 1 );
    $state->{'inXML'} = 0
	unless ( exists($state->{'inXML'}) );
    my $directive;
    my $xmlText;
    my $depth = 0;
    while ( my $line = <$fileHandle> ) {
	# Detect the end of an XML section and change state.
	$state->{'inXML'} = 0
	    if ( $line =~ m/^\s*!!\]/ );
	# Skip instrumentation lines.
	next
	    if ( $line =~ m/^\!\-\->/ );
	# If we're actively processing XML content, accumulate the text.	
	if ( $state->{'inXML'} || $depth > 0 ) {
	    $line =~ s/^(\!\<)?\s*//
		if ( $state->{'inXML'} );
	    (my $processedLine = $line) =~ s/&nbsp;/ /g;
	    $xmlText .= $processedLine;
	    $depth += () = ( $line =~ /<([a-zA-Z0-9]+)[^\/>]*>/g ); # Increment depth by count of any opening elements.
	    $depth -= () = ( $line =~ /<\/([a-zA-Z0-9]+)>/g      ); # Decrement depth by count of any closing elements.
	    if ( defined($xmlText) && $xmlText !~ m/^\s*$/ && $depth == 0 ) {
		# Parse the XML.
		my $xml    = new XML::Simple(KeepRoot => 1);
		$directive = eval{$xml->XMLin($xmlText)};
		if ( $@ ) {
		    print "Extract_Directive: while extracting directive '".$directiveName."' from line ".$fileHandle->input_line_number()." of file '".&File::Names::Get_Name($fileHandle)."' failed parsing with message:\n".$@."\n";
		    print " XML content was:\n";
		    print $xmlText;
		    die();
		}
		if ( $directiveName eq "*" || exists($directive->{$directiveName}) ) {
		    my $matchedDirectiveName = (keys(%{$directive}))[0];
		    $directive = $directive->{$matchedDirectiveName};
		    # Include the root element name if requested.
		    $directive->{'rootElementType'} = $matchedDirectiveName
			if ( exists($options{'setRootElementType'}) && $options{'setRootElementType'} == 1 );
		    # Check any conditions.
		    if ( %options && exists($options{'conditions'}) ) {
			my $matched = 1;
			foreach ( keys(%{$options{'conditions'}}) ) {
			    unless ( exists($directive->{$_}) && $directive->{$_} eq ${options{'conditions'}}{$_} ) {
				$matched = 0;
				last;
			    }
			}
			next
			    unless ( $matched );
		    }
		    last;
		}
		undef($directive);
		undef($xmlText  );
	    }
	}
	# Detect the start of an XML section and change state.
	$state->{'inXML'} = 1
	    if ( $line =~ m/^\s*!!\[/ );
    }
    # Return the directive.
    return $directive;
}

sub Extract_Directives {
    # Extract all named directives from a file.
    my $fileName      = shift();
    my $directiveName = shift();
    (my %options)     = @_
	if ( scalar(@_) > 1 );
    # Open the file and iterate until all directives are obtained.
    my @directives;
    return
	unless ( -e $fileName );
    my %state;
    open(my $fileHandle,$fileName);
    while ( my $directive = &Extract_Directive($fileHandle,$directiveName,\%state,%options) ) {
	push(@directives,$directive);
    }
    close($fileHandle);
    # Return the list of directives.
    return @directives;
}

1;
