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
use Galacticus::Build::Dependencies;

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
    # Get state storable information if we do not have it.
    our $stateStorables    = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml"     )
	unless ( $stateStorables     );
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
    my $node            = $tree;
    my $depth           = 0;
    my $isFunctionClass = 0;
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
		my $functionClass;
		my $after;
		my $before;
		my $useGlobal;
		while ( $eventHookedNode ) {
		    # Capture the module name directly.
		    $moduleName = $eventHookedNode->{'name'}
		        if ( $eventHookedNode->{'type'} eq "module" );
		    # Capture the details of the hooked function.
		    if ( $eventHookedNode->{'type'} eq $node->{'directive'}->{'name'} ) {
			$functionName = $eventHookedNode->{'directive'}->{'function' };
			$after        = $eventHookedNode->{'directive'}->{'after'    }
			    if ( exists($eventHookedNode->{'directive'}->{'after'    }) );
			$before       = $eventHookedNode->{'directive'}->{'before'   }
			    if ( exists($eventHookedNode->{'directive'}->{'before'   }) );
			$useGlobal    = $eventHookedNode->{'directive'}->{'useGlobal'}
			    if ( exists($eventHookedNode->{'directive'}->{'useGlobal'}) );
		    }
		    # Capture any functionClass instances.
		    $functionClass = $eventHookedNode->{'type'}."Class"
		        if ( exists(${$stateStorables->{'functionClasses'}}{$eventHookedNode->{'type'}."Class"}) );
		    # Move to the next node.
		    $eventHookedNode = &Galacticus::Build::SourceTree::Walk_Tree($eventHookedNode,\$eventHookedDepth);
		}
		# If we have a file containing a functionClass instance. Any function it provides for static hooking will be
		# available through the associated functionClass module.
		$moduleName = ${$stateStorables->{'functionClasses'}}{$functionClass}->{'module'}
		   if ( ! defined($moduleName) && defined($functionClass) );
		die("Galacticus::Build::SourceTree::Process::EventHooksStatic: unable to find module containing hooked function")
		    unless ( defined($moduleName  ) );
		die("Galacticus::Build::SourceTree::Process::EventHooksStatic: unable to find name of hooked function"          )
		    unless ( defined($functionName) );
		# Handle useGlobal: append "_" to function name and import from Functions_Global module.
		if ( defined($useGlobal) && $useGlobal eq "yes" ) {
		    $functionName .= "_";
		    $moduleName    = "Functions_Global";
		}
		push(@hookedFunctions,{function => $functionName, module => $moduleName, after => $after, before => $before});
	    }
	    # Sort hooked functions by dependency order if any have after/before constraints.
	    if ( grep { defined($_->{'after'}) || defined($_->{'before'}) } @hookedFunctions ) {
		my %tasks;
		foreach my $hf ( @hookedFunctions ) {
		    $tasks{$hf->{'function'}} = {};
		    $tasks{$hf->{'function'}}->{'after' } = $hf->{'after' }
		        if ( defined($hf->{'after' }) );
		    $tasks{$hf->{'function'}}->{'before'} = $hf->{'before'}
		        if ( defined($hf->{'before'}) );
		}
		my %sortData;
		&Galacticus::Build::Dependencies::Dependency_Sort(\%tasks, \%sortData);
		my %functionByName = map { $_->{'function'} => $_ } @hookedFunctions;
		@hookedFunctions = map { exists($functionByName{$_}) ? $functionByName{$_} : () } @{$sortData{'unitNames'}};
	    } else {
		# Simply sort by function name to ensure consistent build.
		@hookedFunctions = sort {$a->{'function'} cmp $b->{'function'}} @hookedFunctions;
	    }
	    # Get callWith and onReturn from the directive.
	    my $callWith = exists($node->{'directive'}->{'callWith'}) ? $node->{'directive'}->{'callWith'} : '';
	    my $onReturn = exists($node->{'directive'}->{'onReturn'}) ? $node->{'directive'}->{'onReturn'} : '';
	    # Insert the code.
	    my @callLines;
	    foreach my $hf ( @hookedFunctions ) {
		push(@callLines, "call ".$hf->{'function'}."(".$callWith.")");
		push(@callLines, $onReturn)
		    if ( $onReturn ne '' );
	    }
	    my $newNode =
	    {
		type       => "code",
		content    => join("\n",@callLines)."\n",
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
	# Record if the file contains a functionClass instance.
	$isFunctionClass = 1
	    if ( exists(${$stateStorables->{'functionClasses'}}{$node->{'type'}."Class"}) );
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Set visibilities for hooked functions.
    if ( scalar(@functions) > 0 ) {
	if ( defined($moduleNode) ) {
	    &Galacticus::Build::SourceTree::SetVisibility($moduleNode,$_,"public")
		foreach ( @functions );
	} elsif ( ! $isFunctionClass ) {
	    # If no module was found, and this file does not contain a functionClass instance, we have no way to set function
	    # visibility - this is an error. (For files containing functionClass instances, visibilities must have been set
	    # directly in the file in order to cause the functions to be available through the associated module, so we do not
	    # need to act on them here.)
	    die("Galacticus::Build::SourceTree::Process::EventHooksStatic::Process_EventHooksStatic(): hooked ".(scalar(@functions) > 1 ? "functions" : "function")." ".join(", ",map {"'".$_."'"} @functions)." ".(scalar(@functions) > 1 ? "are" : "is")." not in a module");
	}
    }
}

1;
