# Contains a Perl module which implements hooks for the tree preprocessor system.

package Hooks;
use strict;
use warnings;
use utf8;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 

# Define a hash into which modules can insert their hooks.
our %parseHooks;
our %processHooks;

1;
