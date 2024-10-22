# Contains a Perl module which implements processing of dependency directives.

package Galacticus::Build::SourceTree::Process::Dependencies;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'dependencies'} = \&Process_Dependencies;

sub Process_Dependencies {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for hook directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Handle dependenciesInitialize directives by building initialization calls.
	if ( $node->{'type'} eq "dependenciesInitialize" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Parse the dependency file.
	    my %dependencies;
	    open(my $dependencyFile,$ENV{'GALACTICUS_EXEC_PATH'}."/aux/dependencies.yml");
	    while ( my $line = <$dependencyFile> ) {
		if ( $line =~ m/^(.+):\s+([0-9\.]+)/ ) {
		    $dependencies{$1} = $2;
		} else {
		    die("can not parse dependency file line:\n".$line);
		}
	    }
	    close($dependencyFile);
	    # Generate code.
	    my $code = join("\n",map {"call dependencies_%set(var_str('".$_."'),var_str('".$dependencies{$_}."'))"} sort(keys(%dependencies)))."\n";
	    # Insert our code.
	    my $initializationNode =
	    {
		type       => "code",
		content    => $code,
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_Dependencies()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$initializationNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
