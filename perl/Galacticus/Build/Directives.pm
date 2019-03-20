# Contains a Perl module which implements processing of directives.

package Galacticus::Build::Directives;
use strict;
use warnings;

sub Extract_Directive {
    # Extract a named directive from a given file handle.
    my $fileHandle    = shift();
    my $directiveName = shift();
    (my %options)     = @_
	if ( scalar(@_) > 1 );
    my $directive;
    my $xmlText;
    my $commentRegEx       = exists($options{'comment'}) ? $options{'comment'} : qr/^\s*\!#/;
    my $xmlTagRegEx        = $directiveName eq "*" ? qr/^\s*<([a-zA-Z0-9]+).*>/ : qr/^\s*<($directiveName)(>|\s.*>)/; 
    my $depth              = 0;
    while ( my $line = <$fileHandle> ) {
	if ( $line =~ s/$commentRegEx// || $depth > 0 ) {
	    $xmlText .= $line;
	    $depth += () = ( $line =~ /<([a-zA-Z0-9]+)[^\/>]*>/g ); # Increment depth by count of any opening elements.
	    $depth -= () = ( $line =~ /<\/([a-zA-Z0-9]+)>/g      ); # Decrement depth by count of any closing elements.
	    if ( defined($xmlText) && $depth == 0 ) {
		# Parse the XML.
		my $xml    = new XML::Simple(KeepRoot => 1);
		$directive = $xml->XMLin($xmlText);
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
    open(my $fileHandle,$fileName);
    while ( my $directive = &Extract_Directive($fileHandle,$directiveName,%options) ) {
	push(@directives,$directive);
    }
    close($fileHandle);
    # Return the list of directives.
    return @directives;
}

1;
