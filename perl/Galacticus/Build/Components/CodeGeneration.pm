# Contains a Perl module which implements functions required for code generation.

package Galacticus::Build::Components::CodeGeneration;
use strict;
use warnings;
use utf8;

sub Function_Arguments {
    # Return a list of all variables which are arguments.
    my @data = @{shift()};
    my @arguments;
    foreach my $datum ( @data ) {
	push(@arguments,@{$datum->{'variables'}})
	    if ( exists($datum->{'attributes'}) && grep {$_ =~ m/intent/} @{$datum->{'attributes'}} );
    }
    return @arguments;
}

sub Importables {
    # Return a list of all types which must be imported into an interface.
    my @data = @{shift()};
    my @importables;
    foreach my $datum ( @data ) {
	push(@importables,$datum->{'type'})
	    if ( $datum->{'intrinsic'} eq "class" || $datum->{'intrinsic'} eq "type" );
    }
    return @importables;
}

1;
