#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::SourceTree;
use Galacticus::Options;
use File::Copy;

# Format module use sections of source files.
# Andrew Benson (02-February-2021)

# Get arguments.
die("Usage: formatModuleUses.pl <infile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $inputFileName  = $ARGV[0];
my %options =
    (
     suffix => "~"
     );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse the file to build a tree.
my $tree = &Galacticus::Build::SourceTree::ParseFile($inputFileName, instrument => 0);

# Walk the tree and regenerate any module use blocks.
my $node  = $tree;
my $depth = 0;
while ( $node ) {
    if ( $node->{'type'} eq "moduleUse" ) {
	&Galacticus::Build::SourceTree::Parse::ModuleUses::UpdateUses($node);
    }
    $node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
}    


# Serialize back to source code.
move($inputFileName,$inputFileName.$options{'suffix'});
open(my $outputFile,">",$inputFileName);
print $outputFile &Galacticus::Build::SourceTree::Serialize($tree, annotate => 0);
close($outputFile);

exit;
