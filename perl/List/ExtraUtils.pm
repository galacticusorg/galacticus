# Contains a Perl module which implements extra list utilities.

package ExtraUtils;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
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

1;
