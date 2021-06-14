# Contains a Perl module which implements hooks for the tree preprocessor system.

package Galacticus::Build::SourceTree::Hooks;
use strict;
use warnings;
use utf8;

# Define a hash into which modules can insert their hooks.
our %parseHooks;
our %processHooks;
our %processDependencies;
our %analyzeHooks;

1;
