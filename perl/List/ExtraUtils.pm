# Contains a Perl module which implements extra list utilities.

package ExtraUtils;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Scalar::Util 'reftype';

sub smart_push {
    # Intelligently push an object onto an array.
    my $array = shift;
    my $list  = shift;
    # Proceed only if the list is defined.
    if ( defined($list) ) {
	# Is the list a reference?
	if ( reftype($list) ) {
	    # It is, check if it's an array?
	    if ( reftype($list) eq "ARRAY" ) {
		# It is an array, simply push as an array.
		push(@{$array},@{$list});
	    } else {
		# It's not an array, push the object.
		push(@{$array},  $list );
	    }
	} else {
	    # List is not a reference, push the object.
	    push    (@{$array},  $list );
	}
    }
}

sub as_array {
    # Return an object as an array.
    my @array;
    my $list  = shift;
    &smart_push(\@array,$list);
    return @array;
}

sub hashList {
    # Return a list containing all elements of the given hash.
    my $hashRef = shift();
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );
    my @list    = map {
	$hashRef->{$_}->{$options{'keyAs'}} = $_ if ( exists($options{'keyAs'}) );
	$hashRef->{$_}
    } sort(keys(%{$hashRef}));
}

sub sortedKeys {
    # Return a sorted list of hash keys.
    my $hashRef = shift();
    return ()
	unless ( $hashRef );	
    return sort(keys(%{$hashRef}));
}

1;
