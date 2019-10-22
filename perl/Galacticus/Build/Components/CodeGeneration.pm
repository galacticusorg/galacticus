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
	    if ( 
		(
		 $datum->{'intrinsic'} eq "class"
		 &&
		 $datum->{'type'     } ne "*"
		)
		||
		$datum ->{'intrinsic'} eq "type"
		)
	    {
		push(@importables,$datum->{'type'});
	    } elsif (
		$datum->{'intrinsic'} ne "class" 
		&&
		exists($datum->{'type'})
		&&
		$datum->{'type'} !~ m/\s*\s*len\s*=/
		&&
		$datum->{'type'} =~ m/\s*(kind\s*=\s*)??([a-zA-Z0-9_]+)/
		) {
		push(@importables,$2);
	    }	
    }
    return @importables;
}

1;
