#!/usr/bin/env python3
import subprocess
import sys
import os

# Test extraction of directives from source files.
# Andrew Benson (ported to Python)

execPath = os.environ.get("GALACTICUS_EXEC_PATH", "..")

# Run the Perl module via subprocess.
result = subprocess.run(
    f'perl -e \'use lib "{execPath}/perl"; use Galacticus::Build::Directives; '
    f'my @directives = &Galacticus::Build::Directives::Extract_Directives("{execPath}/testSuite/data/sourceWithDirectives.F90","*",setRootElementType => 1); '
    f'my $status = scalar(@directives) == 10 ? "success" : "FAILED"; '
    f'if ($status eq "success") {{ print $status.": all directives found\\n"; }} '
    f'else {{ print $status.": all directives NOT found\\n"; }}\'',
    shell=True
)
