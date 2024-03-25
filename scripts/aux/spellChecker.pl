#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use LaTeX::SpellCheck;

# Perform spell checking of files
# Andrew Benson (28-February-2023)

# Get the file to process.
die('Usage: spellChecker.pl <fileName> <warningFile>')
    unless ( scalar(@ARGV) == 2 );
my $fileName        = $ARGV[0];
my $warningFileName = $ARGV[1];

# Check the file.
my $warnings = &LaTeX::SpellCheck::spellCheckFile($fileName,$fileName);

# Write accumulated warnings to file (appending).
unless ( $warnings eq "" ) {
    open(my $warningFile,">>",$warningFileName);
    print $warningFile $warnings;
    close($warningFile);
}

exit 0;
