# Contains a Perl module which implements hooks for the component build system.

package Utils;
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
use utf8;

# Define a hash into which modules can insert their hooks.
our %componentHooks;

# Global verbosity level.
our $verbosityLevel = 1;

# Boolean labels.
our @booleanLabel = ( "false", "true" );

# Intrinsic types.
our %intrinsicTypes =
    (
     "integer"     => "integer"                ,
     "longInteger" => "integer(kind=kind_int8)",
     "logical"     => "logical"                ,
     "double"      => "double precision"       ,
     "void"        => "void"
    );

sub isIntrinsic {
    # Return true if the given type matches an intrinsic type.
    my $type = shift();
    return (grep {$_ eq $type} keys(%intrinsicTypes)) == 1 ? 1 : 0;
}

sub isOutputIntrinsic {
    # Return true if the given type matches an outputtable intrinsic type.
    my $type = shift();
    return (grep {$_ eq $type} ("double","integer","longInteger")) == 1 ? 1 : 0;
}

1;
