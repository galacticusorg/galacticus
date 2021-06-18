#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Fortran::Utils;

# Test the Fortran::Utils modules.
# Andrew Benson (17-June-2021)

# Regression test for issue #169 (https://github.com/galacticusorg/galacticus/issues/169), wrong code produced from valid variable declarations.
my $variables =
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponentTest",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "label" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "name" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "isEvolvable" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 attributes => [ "allocatable", "dimension(:)" ],
		 variables  => [ "labelsTmp", "namesTmp" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "allocatable", "dimension(:)" ],
		 variables  => [ "evolvableTmp" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "found" ]
	     }
	    ];
my $code = &Fortran::Utils::Format_Variable_Definitions($variables);
if ( $code =~ m/,\s*intent\(in   \)\s*,\s*optional\s*::/ ) {
    print "success: issue #169 (https://github.com/galacticusorg/galacticus/issues/169) - wrong code produced from valid variable declarations\n";
} else {
    print "FAILED: issue #169 (https://github.com/galacticusorg/galacticus/issues/169) - wrong code produced from valid variable declarations\n";
}

exit;
