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
    my $xmlTagRegEx  = $directiveName eq "*" ? qr/^\s*<([a-zA-Z]+)[\s>]/ : qr/^\s*<($directiveName)[\s>]/;    
    while ( my $line = <$fileHandle> ) {
	if ( $line =~ s/$commentRegEx// ) {
	    if ( $line =~ $xmlTagRegEx ) {
		$matchedDirectiveName = $1;
		$xmlText              = "";
	    }
	    $xmlText .= $line;
	    if ( $line =~ m/^\s*<\/$matchedDirectiveName[\s>]/ || $line =~ m/^\s*<$matchedDirectiveName\s.*\/>/ ) {
		# Parse the XML.
		my $xml    = new XML::Simple();
		$directive = $xml->XMLin($xmlText);
		if ( %options && exists($options{'conditions'}) ) {
		    foreach ( keys(%{$options{'conditions'}}) ) {
			unless ( exists($directive->{$_}) && $directive->{$_} eq ${options{'conditions'}}{$_} ) {
			    undef($directive);
			    next;
			}
		    }
		}
		last;
	    }
	}
    }
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
