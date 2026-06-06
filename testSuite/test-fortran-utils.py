#!/usr/bin/env python3
import subprocess
import sys
import os

# Test the Fortran::Utils modules.
# Andrew Benson (ported to Python)

execPath = os.environ.get("GALACTICUS_EXEC_PATH", "..")

# Run the Perl test via subprocess (testing Fortran::Utils::Format_Variable_Definitions).
perlCode = r"""
use lib "$execPath/perl";
use Fortran::Utils;
my $variables = [
    { intrinsic => "class",     type => "nodeComponentTest", attributes => ["intent(inout)"], variables => ["self"] },
    { intrinsic => "type",      type => "varying_string",    attributes => ["intent(in   )"], variables => ["label"] },
    { intrinsic => "character", type => "len=*",             attributes => ["intent(in   )"], variables => ["name"] },
    { intrinsic => "logical",                                 attributes => ["intent(in   )", "optional"], variables => ["isEvolvable"] },
    { intrinsic => "type",      type => "varying_string",    attributes => ["allocatable", "dimension(:)"], variables => ["labelsTmp", "namesTmp"] },
    { intrinsic => "logical",                                attributes => ["allocatable", "dimension(:)"], variables => ["evolvableTmp"] },
    { intrinsic => "logical",                                variables => ["found"] }
];
my $code = &Fortran::Utils::Format_Variable_Definitions($variables);
if ($code =~ m/,\s*intent\(in   \)\s*,\s*optional\s*::/) {
    print "success: issue #169 (https://github.com/galacticusorg/galacticus/issues/169) - wrong code produced from valid variable declarations\n";
} else {
    print "FAILED: issue #169 (https://github.com/galacticusorg/galacticus/issues/169) - wrong code produced from valid variable declarations\n";
}
"""

perlCode = perlCode.replace("$execPath", execPath)
result = subprocess.run(["perl", "-e", perlCode])
