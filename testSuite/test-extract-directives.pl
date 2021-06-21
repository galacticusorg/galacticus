#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::Directives;
use Data::Dumper;

# Test extraction of directives from source files.
# Andrew Benson (10-June-2021)

my @directives = &Galacticus::Build::Directives::Extract_Directives($ENV{'GALACTICUS_EXEC_PATH'}."/testSuite/data/sourceWithDirectives.F90","*",setRootElementType => 1);
my $status = scalar(@directives) == 6 ? "success" : "FAILED";
print $status.": all directives found\n";

exit;
