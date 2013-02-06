# Contains a Perl module which implements hooks for the build system.

package Hooks;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;

# Define a hash into which modules can insert their hooks.
our %moduleHooks;

1;
