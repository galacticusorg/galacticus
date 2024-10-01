# Contains a Perl module which implements processing of parameter migration directives.

package Galacticus::Build::SourceTree::Process::ParameterMigration;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'parameterMigration'} = \&Process_ParameterMigration;

sub Process_ParameterMigration {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get code directive locations.
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
    # Walk the tree, looking for hook directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Handle parameterMigration directives by building eventHook objects for all events.
	if ( $node->{'type'} eq "parameterMigration" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Find all known migrations.
	    my $xml          = new XML::Simple();
	    my $migrations   = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/scripts/aux/migrations.xml");
	    my @commitHashes = map {$_->{'commit'}} &List::ExtraUtils::as_array($migrations->{'migration'});
	    # Generate code to initialize an array of commit hashes.
	    my $code;
	    for(my $i=0;$i<scalar(@commitHashes);++$i) {
		$code .= "commitHash(".($i+1).")=\"".$commitHashes[$i]."\"//c_null_char\n";
	    }
	    my $newNode   =
	    {
		type       => "code",
		content    => $code,
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_ParameterMigration()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Insert a variable for the commit hash array.
	    my @declarations =
		(
		 {
		     intrinsic    => "character",
		     type         => "len=41",
		     attributes   => [ "dimension(".scalar(@commitHashes).")" ],
		     variables    => [ "commitHash" ],
		     preprocessor => "GIT2AVAIL"
		 }
		);
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
