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
    my $matchedDirectiveName;
    my $commentRegEx = exists($options{'comment'}) ? $options{'comment'} : qr/^\s*\!#/;
    my $xmlTagRegEx  = $directiveName eq "*" ? qr/^\s*<([a-zA-Z0-9]+).*>/ : qr/^\s*<($directiveName)(>|\s.*>)/;    
    while ( my $line = <$fileHandle> ) {
	if ( $line =~ s/$commentRegEx// ) {
	    $xmlText .= $line;
	    if ( ! defined($matchedDirectiveName) && $line =~ $xmlTagRegEx ) {
		$matchedDirectiveName = $1;
		$xmlText              = $line;
	    }
	    if ( defined($matchedDirectiveName) && ( $line =~ m/^\s*<\/$matchedDirectiveName[\s>]/ || $line =~ m/^\s*<$matchedDirectiveName\s.*\/>/ ) ) {
		# Parse the XML.
		my $xml    = new XML::Simple();
		$directive = $xml->XMLin($xmlText);
		if ( %options && exists($options{'conditions'}) ) {
		    my $matched = 1;
		    foreach ( keys(%{$options{'conditions'}}) ) {
			unless ( exists($directive->{$_}) && $directive->{$_} eq ${options{'conditions'}}{$_} ) {
			    $matched = 0;
			    last;
			}
		    }
		    unless ( $matched ) {
			undef($directive           );
			undef($xmlText             );
			undef($matchedDirectiveName);
			next;
		    }
		}
		last;
	    }
	} elsif ( defined($matchedDirectiveName) ) {
	    $xmlText .= $line
	}
    }
    # Include the root element name if requested.
    $directive->{'rootElementType'} = $matchedDirectiveName
	if ( defined($matchedDirectiveName) && exists($options{'setRootElementType'}) && $options{'setRootElementType'} == 1 );
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
