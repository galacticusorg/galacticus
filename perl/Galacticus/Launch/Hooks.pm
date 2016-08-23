# Contains a Perl module which implements hooks for the launch system.

package Galacticus::Launch::Hooks;
use strict;
use warnings;
use utf8;

# Define a hash into which modules can insert their hooks.
our %moduleHooks;

1;
