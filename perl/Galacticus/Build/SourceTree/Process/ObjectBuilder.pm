# Contains a Perl module which implements processing of object builder directives.

package Galacticus::Build::SourceTree::Process::ObjectBuilder;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use List::ExtraUtils;
use XML::Simple;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'objectBuilder'} = \&Process_ObjectBuilder;

sub Process_ObjectBuilder {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $xml   = new XML::Simple();
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "objectBuilder" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the object builder. The logic here is that we search for
	    # a matching parameter in the given parameter set. If none is found, we step up
	    # through parent parameters until we do find one or we reach the top-level parameter
	    # set. If a match is found, we use that definition to build our object, unless this
	    # is the global parameter set and we're at the top level (in which case we use the
	    # default object of the relevant class). If no object is found, we again use the
	    # default object of the class if this is the global parameter set, otherwise we
	    # abort. If using a specific, given definition to build the object, we first check
	    # if it has already been built, reusing if it has, and building and storing if it
	    # has not. This prevents creating instances more than once when not necessary.
	    my $parameterName = exists($node->{'directive'}->{'parameterName'}) ? $node->{'directive'}->{'parameterName'} : $node->{'directive'}->{'class'}."Method";
	    my $defaultName   = $parameterName eq $node->{'directive'}->{'class'}."Method";
	    die("objects with defaults must have explicit parameter names")
		if ( exists($node->{'directive'}->{'default'}) && ! exists($node->{'directive'}->{'parameterName'}) );
	    my $parametersDefaultRequired = 0;
	    my $builderCode;
	    $builderCode .= "   ! Determine where to build+store or point to the required object....\n";
	    $builderCode .= "   parametersObjectBuildIsPrivate=.".($node->{'directive'}->{'threadPrivate'} eq "yes" ? "true" : "false").".\n"
		if ( exists($node->{'directive'}->{'threadPrivate'}) );
	    $builderCode .= "   parametersCurrent => ".$node->{'directive'}->{'source'}."\n";
	    if ( exists($node->{'directive'}->{'parameterName'}) ) {
		if ( exists($node->{'directive'}->{'default'}) ) {
		    $parametersDefaultRequired  =  1;
		    my $defaultXML              =  $xml->XMLout($node->{'directive'}->{'default'}, RootName => "parameters");
		    $defaultXML                 =~ s/\s*\n\s*//g;
		    $defaultXML                 =~ s/\s{2,}/ /g;
		    $builderCode               .=  "   if (.not.parametersCurrent%isPresent('".$parameterName."')) then\n";
		    $builderCode               .=  "    parametersDefault=inputParameters(var_str('".$defaultXML."'),noOutput=.true.)\n";
		    $builderCode               .= "     call parametersDefault%parametersGroupCopy(parametersCurrent)\n";
		    $builderCode               .=  "    parametersCurrent => parametersDefault\n";
		    $builderCode               .=  "  end if\n";
		} else {
		    $builderCode               .=  "   if (.not.parametersCurrent%isPresent('".$parameterName."')) call Galacticus_Error_Report('[".$parameterName."] object is undefined'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		}
	    } else {	    
		$builderCode .= "   do while (.not.parametersCurrent%isPresent('".$parameterName."').and.associated(parametersCurrent%parent))\n";
		$builderCode .= "      parametersCurrent => parametersCurrent%parent\n";
		$builderCode .= "   end do\n";
	    }
	    $builderCode .= "   if (parametersCurrent%isPresent('".$parameterName."').and.(.not.".$node->{'directive'}->{'source'}."%isGlobal().or.associated(parametersCurrent%parent))) then\n"
		if ( $defaultName );
	    $builderCode .= "      ! Object should belong to the parameter node. Get the node and test whether the object has already been created in it.\n";
	    $builderCode .= "      parameterNode => parametersCurrent%node('".$parameterName."'".(exists($node->{'directive'}->{'copy'}) ? ",copyInstance=".$node->{'directive'}->{'copy'} : "").")\n";
	    $builderCode .= "      if (parameterNode%objectCreated()) then\n";
	    $builderCode .= "         ! Object already exists - simply get a pointer to it.\n";
	    $builderCode .= "         genericObject => parameterNode%objectGet()\n";
	    $builderCode .= "         select type (genericObject)\n";
	    $builderCode .= "         class is (".$node->{'directive'}->{'class'}."Class)\n";
	    $builderCode .= "            ".$node->{'directive'}->{'name'}." => genericObject\n";
	    $builderCode .= "         class default\n";
	    $builderCode .= "            call Galacticus_Error_Report('parameter-stored object is not of [$node->{'directive'}->{'class'}] class'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
	    $builderCode .= "         end select\n";
	    $builderCode .= "      else\n";
	    $builderCode .= "         ! Object does not yet exist - build it and store in the parameter node.\n";
	    $builderCode .= "         ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."(parametersCurrent".(exists($node->{'directive'}->{'copy'}) ? ",copyInstance=".$node->{'directive'}->{'copy'} : "").(exists($node->{'directive'}->{'parameterName'}) ? ",parameterName='".$parameterName."'" : "").")\n";
	    $builderCode .= "         call parameterNode%objectSet(".$node->{'directive'}->{'name'}.")\n";
	    $builderCode .= "      end if\n";
	    if ( $defaultName ) {
		$builderCode .= "   else if (".$node->{'directive'}->{'source'}."%isGlobal()) then\n";
		$builderCode .= "      ! This is the global parameter set - so we can use the default object of this class.\n";
		$builderCode .= "      ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."()\n";
		$builderCode .= "   else\n";
		$builderCode .= "      ! No means to define the object.\n";
		$builderCode .= "      call Galacticus_Error_Report('[".$parameterName."] object is undefined'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		$builderCode .= "   end if\n";
	    }
	    # Build a code node.
	    my $newNode =
	    {
		type       => "code"      ,
		content    => $builderCode,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Add new variables and attributes.
	    unless ( exists($node->{'parent'}->{'objectBuilderDeclarations'}) ) {
		my @declarations =
		    (
		     {
			 intrinsic  => "type"                 ,
			 type       => "inputParameters"      ,
			 variables  => [ "parametersCurrent" ],
			 attributes => [ "pointer"           ]
		     },
		     {
			 intrinsic  => "type"                 ,
			 type       => "inputParameter"       ,
			 variables  => [ "parameterNode"     ],
			 attributes => [ "pointer"           ]
		     },
		     {
			 intrinsic  => "class"                ,
			 type       => "*"                    ,
			 variables  => [ "genericObject"     ],
			 attributes => [ "pointer"           ]
		     }
		    );
		&Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
		# Ensure error reporting module is used.
		my $usesNode =
		{
		    type      => "moduleUse",
		    moduleUse =>
		    {
			Galacticus_Error =>
			{
			    intrinsic => 0,
			    all       => 1
			}
		    }
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
		# Record that we have added the necessary declarations to the parent.
		$node->{'parent'}->{'objectBuilderDeclarations'} = 1;
	    }
	    if ( $parametersDefaultRequired && ! exists($node->{'parent'}->{'objectBuilderDefaultDeclarations'}) ) {
		my @declarations =
		    (
		     {
			 intrinsic  => "type"                 ,
			 type       => "inputParameters"      ,
			 variables  => [ "parametersDefault" ],
			 attributes => [ "target"            ]
		     }
		    );
		&Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
		# Record that we have added the necessary declarations to the parent.
		$node->{'parent'}->{'objectBuilderDefaultDeclarations'} = 1;
	    }
	    unless ( exists($node->{'parent'}->{'objectBuilderAttributes'}->{$node->{'directive'}->{'source'}}) ) {
		&Galacticus::Build::SourceTree::Parse::Declarations::AddAttributes($node->{'parent'},$node->{'directive'}->{'source'},["target"])
		    if ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node->{'parent'},$node->{'directive'}->{'source'}) );
		$node->{'parent'}->{'objectBuilderAttributes'}->{$node->{'directive'}->{'source'}} = 1;
	    }
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	if ( $node->{'type'} eq "objectDestructor" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the object destructor. The pointer should only be deallocated if associated and
	    # finalizable. Otherwise simply nullify it.
	    my $destructorCode;
	    $destructorCode .= "if (associated(".$node->{'directive'}->{'name'}.").and.".$node->{'directive'}->{'name'}."%isFinalizable()) then\n";
	    $destructorCode .= "   ! Deallocate the pointer.\n";
	    $destructorCode .= "   deallocate(".$node->{'directive'}->{'name'}.")\n";
	    $destructorCode .= "else\n";
	    $destructorCode .= "   ! Nullify the pointer.\n";
	    $destructorCode .= "   nullify(".$node->{'directive'}->{'name'}.")\n";
	    $destructorCode .= "end if\n";
	    # Build a code node.
	    my $newNode =
	    {
		type       => "code"      ,
		content    => $destructorCode,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
