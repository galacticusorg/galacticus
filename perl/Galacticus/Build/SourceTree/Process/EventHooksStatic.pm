# Contains a Perl module which implements processing of static event directives.

package Galacticus::Build::SourceTree::Process::EventHooksStatic;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use Storable;
use List::ExtraUtils;
use Galacticus::Build::Directives;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'eventHooksStatic'} = \&Process_EventHooksStatic;

sub Process_EventHooksStatic {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get code directive locations if we do not have them.
    our $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml")
	unless ( $directiveLocations );
    # Initialize event hook names.
    our @eventHookNames;
    unless ( @eventHookNames ) {
	if ( -e $ENV{'BUILDPATH'}."/eventHooksStatic.blob" ) {
	    @eventHookNames = @{retrieve($ENV{'BUILDPATH'}."/eventHooksStatic.blob")};
	} else {	  
	    @eventHookNames = map {$_->{'name'}} map {&Galacticus::Build::Directives::Extract_Directives($_,'eventHookStatic')} &List::ExtraUtils::as_array($directiveLocations->{'eventHookStatic'}->{'file'});
	    store(\@eventHookNames,$ENV{'BUILDPATH'}."/eventHooksStatic.blob");
	}
    }
    # Walk the tree, looking for hook directives.
    my $node  = $tree;
    my $depth = 0;
    my $moduleNode;
    my @functions;
    while ( $node ) {
	# Locate the containing module.
	$moduleNode = $node
	    if ( $node->{'type'} eq "module" );	     
	# Handle eventHookStatic directives by creating code to call any hooked functions.
	if ( $node->{'type'} eq "eventHookStatic" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Find locations of all matching directives.
	    my @eventHookedLocations = &List::ExtraUtils::as_array($directiveLocations->{$node->{'directive'}->{'name'}}->{'file'})
		if ( exists($directiveLocations->{$node->{'directive'}->{'name'}}) );
	    # Iterate over event hook locations, finding names of hooked functions and the modules in which they are contained.
	    my @hookedFunctions;
	    foreach my $eventHookedLocation ( @eventHookedLocations ) {
		my $eventHookedTree  = &Galacticus::Build::SourceTree::ParseFile($eventHookedLocation);
		my $eventHookedNode  = $eventHookedTree;
		my $eventHookedDepth = 0;
		my $moduleName;
		my $functionName;
		while ( $eventHookedNode ) {
		    $moduleName = $eventHookedNode->{'name'}
		        if ( $eventHookedNode->{'type'} eq "module" );
		    $functionName = $eventHookedNode->{'directive'}->{'function'}
		        if ( $eventHookedNode->{'type'} eq $node->{'directive'}->{'name'} );
		    $eventHookedNode = &Galacticus::Build::SourceTree::Walk_Tree($eventHookedNode,\$eventHookedDepth);
		}
		die("Galacticus::Build::SourceTree::Process::EventHooksStatic: unable to find module containing hooked function")
		    unless ( defined($moduleName  ) );
		die("Galacticus::Build::SourceTree::Process::EventHooksStatic: unable to find name of hooked function"          )
		    unless ( defined($functionName) );
		push(@hookedFunctions,{function => $functionName, module => $moduleName});
	    }
	    # Insert the code.
	    my $newNode =
	    {
		type       => "code",
		content    => join("\n",map {"call ".$_->{'function'}."()"} @hookedFunctions)."\n",
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooksStatic::Process_EventHooksStatic()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse => {},
		source     => "Galacticus::Build::SourceTree::Process::EventHooksStatic::Process_EventHooksStatic()",
		line       => 1
	    };
	    $usesNode->{'moduleUse'}->{$_->{'module'}} = {intrinsic => 0, only => {$_->{'function'} => 1}}
	        foreach ( @hookedFunctions );
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	}
	# Look for event hooked functions.
	if ( (grep {$node->{'type'} eq $_} @eventHookNames) && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} = 1;
	    push(@functions,$node->{'directive'}->{'function'});
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Set visibilities for hooked functions.
    if ( scalar(@functions) > 0 ) {
	if ( defined($moduleNode) ) {
	    &Galacticus::Build::SourceTree::SetVisibility($moduleNode,$_,"public")
		foreach ( @functions );
	} else {
	    die("Galacticus::Build::SourceTree::Process::EventHooksStatic::Process_EventHooksStatic(): hooked function is not in a module");
	}
    }
}

1;
