# Contains a Perl module which implements processing of global function directives.

package Galacticus::Build::SourceTree::Process::FunctionsGlobal;
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
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'functionsGlobal'} = \&Process_FunctionsGlobal;

sub Process_FunctionsGlobal {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get code directive locations if we do not have them.
    our $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml")
	unless ( $directiveLocations );
    # Get state storables if we do not have them.
    our $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml")
	unless ( $stateStorables );
    # Walk the tree, looking for functionsGlobal directives.
    my $node  = $tree;
    my $depth = 0;
    my $moduleNode;
    while ( $node ) {
	# Locate the containing module.
	$moduleNode = $node
	    if ( $node->{'type'} eq "module" );
	# Handle functionsGlobal directives.
	if ( $node->{'type'} eq "functionsGlobal" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Discover all global functions.
	    my @files = &List::ExtraUtils::as_array($directiveLocations->{'functionGlobal'}->{'file'});
	    my @functionsGlobal;
	    foreach my $file ( @files ) {
		my @directives = &Galacticus::Build::Directives::Extract_Directives($file,'functionGlobal');
		push(@functionsGlobal,map {$_->{'file'} = $file; $_} @directives);
	    }
	    # Handle the different types of directive.
	    if ( $node->{'directive'}->{'type'} eq "pointers" ) {
		# Generate pointer definitions.
		my $pointers = join("\n",map {"    procedure(".$_->{'unitName'}."_Null), pointer :: ".$_->{'unitName'}."_ => ".$_->{'unitName'}."_Null"} @functionsGlobal)."\n";
		# Generate and insert a node.
		my $pointerNode =
		{
		    type       => "code",
		    content    => $pointers,
		    firstChild => undef(),
		    source     => "Galacticus::Build::SourceTree::Process::FunctionsGlobal::Process_FunctionsGlobal()",
		    line       => 1
		};
		&Galacticus::Build::SourceTree::InsertPreContains($node,[$pointerNode]);
		# Generate pointer definitions.
		my $definitions;
		foreach my $functionGlobal ( @functionsGlobal ) {
		    my $opener;
		    my $closer;
		    if ( $functionGlobal->{'type'} eq "void" ) {
			$opener = "subroutine";
			$closer = "subroutine";
		    } else {
			$opener = "function";
			$closer = "function";
		    }
		    my @names;
		    foreach ( &List::ExtraUtils::as_array($functionGlobal->{'arguments'}) ) {
			if ( $_ =~ m/::(.*)$/ ) {
			    push(@names,&Fortran::Utils::Extract_Variables($1));
			}
		    }
		    $definitions .= " ".$opener." ".$functionGlobal->{'unitName'}."_Null(".join(",",@names).")\n";
		    $definitions .= "  use :: Error, only : Error_Report\n";
		    if ( exists($functionGlobal->{'module'}) ) {
			foreach my $module ( &List::ExtraUtils::as_array($functionGlobal->{'module'}) ) {
			    (my $moduleName = $module) =~ s/\s*,.*//;
			    $definitions .= "use".($moduleName eq "ISO_C_Binding" ? ", intrinsic" : "")." :: ".$module."\n";
			}
		    }
		    $definitions .= "  ".$functionGlobal->{'type'}." :: ".$functionGlobal->{'unitName'}."_Null\n"
			if ( $functionGlobal->{'type'} ne "void" );
		    foreach ( &List::ExtraUtils::as_array($functionGlobal->{'arguments'}) ) {
			$definitions .= "  ".$_."\n";
		    }
		    $definitions .= "  !\$GLC attributes unused :: ".join(",",@names)."\n"
			unless ( scalar(@names) == 0 );
		    if ( $functionGlobal->{'type'} eq "double precision" ) {
			$definitions .= "  ".$functionGlobal->{'unitName'}."_Null=0.0d0\n";
		    } elsif ( $functionGlobal->{'type'} =~ m/,\s*pointer/ ) {
			$definitions .= "  ".$functionGlobal->{'unitName'}."_Null => null()\n";
		    }
		    $definitions .= "  call Error_Report('global functions have not been initialized'//{introspection:location})\n";
		    $definitions .= " end ".$closer." ".$functionGlobal->{'unitName'}."_Null\n";
		}
		# To allow processing of directives by our preprocessor, we parse and process our generated code here.
		my $treeTmp = &Galacticus::Build::SourceTree::ParseCode($definitions,'Galacticus::Build::SourceTree::Process::FunctionsGlobal::Process_FunctionsGlobal()');
		&Galacticus::Build::SourceTree::ProcessTree($treeTmp);
		# Insert the preprocessed tree.	      
		&Galacticus::Build::SourceTree::InsertPostContains($node,[$treeTmp]);
	    } elsif ( $node->{'directive'}->{'type'} eq "establish" ) {
		# Generate pointer assignments.
		my $code = join("\n",map {"    ".$_->{'unitName'}."_ => ".$_->{'unitName'}} @functionsGlobal)."\n";
		# Generate module requirements.
		my %moduleUses;
		foreach my $functionGlobal ( @functionsGlobal ) {
		    my $treeTarget = &Galacticus::Build::SourceTree::ParseFile($functionGlobal->{'file'});
		    my $nodeTarget = $treeTarget->{'firstChild'};
		    while ( $nodeTarget ) {
			$moduleUses{$nodeTarget->{'name'}}->{'only'}->{$functionGlobal->{'unitName'}} = 1
			    if ( $nodeTarget->{'type'} eq "module" );
			$moduleUses{${$stateStorables->{'functionClasses'}}{$nodeTarget->{'type'}."Class"}->{'module'}}->{'only'}->{$functionGlobal->{'unitName'}} = 1
			    if ( exists(${$stateStorables->{'functionClasses'}}{$nodeTarget->{'type'}."Class"}) );
			$nodeTarget = $nodeTarget->{'sibling'};
		    }
		}
		my $moduleUseNode =
		{
		    type       => "moduleUse",
		    sibling    => undef()    ,
		    parent     => undef()    ,
		    firstChild => undef()    ,
		    source     => "Galacticus::Build::SourceTree::Process::FunctionsGlobal::Process_FunctionsGlobal()",
		    line       => 1          ,
		    moduleUse  => \%moduleUses
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$moduleUseNode);
		# Generate and insert a node.
		my $newNode =
		{
		    type       => "code",
		    content    => $code,
		    firstChild => undef(),
		    source     => "Galacticus::Build::SourceTree::Process::FunctionsGlobal::Process_FunctionsGlobal()",
		    line       => 1
		};
		&Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    } else {
		die("Galacticus::Build::SourceTree::Process::FunctionsGlobal: unknown type '".$node->{'directive'}->{'type'}."'");
	    }
	}
	# Walk the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
