# Contains a Perl module which implements processing of directives related to construction/destruction/referencing of
# functionClass objects. Specifically five separate directives are supported:
#
#    objectBuilder: This directive will get a functionClass object from the given parameter set. The object will be constructed if
#    necessary (and stored in the inputParameter node for re-use), or will be retrieved from the inputParameter node and
#    re-used. In either case the reference count of the object is incremented.
#
#    objectDestructor: This directive will remove the reference to the given object. The reference count is decremented. If the
#    reference count reaches zero as a result then the object is deallocated, otherwise the pointer to the object is simply
#    nullified.
#
#    referenceCountIncrement: This directive will simply increment the reference count to an object. It is intended for use where
#    an object is passed to a functionClass which wants to retain a reference to that object. The reference count of the object
#    will be incremented.
#
#    referenceAcquire: This directive will acquire a pointer to an object (typically as a result of a function). It is intended
#    for use where some other functionClass object provides a method that returns a pointer to a functionClass object. The
#    reference count to the object will be incremented.
#
#    referenceConstruct: This directive is intended for use when a functionClass object is being constructed directly. The
#    reference count to the object is incremented.
#
# All of these directives supported output of debugging information (detailing the location of the object in question along with
# that of its "owner") which can be analyzed by ./scripts/aux/functionClassReferencesDebug.pl.

package Galacticus::Build::SourceTree::Process::ObjectBuilder;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use XML::Simple;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'objectBuilder'} = \&Process_ObjectBuilder;

sub Process_ObjectBuilder {
    # Get the tree.
    my $tree = shift();
    # Determine if debugging output is required.
    my $debugging = exists($ENV{'GALACTICUS_OBJECTS_DEBUG'}) && $ENV{'GALACTICUS_OBJECTS_DEBUG'} eq "yes";
    # Walk the tree, looking for code blocks.
    my $node               = $tree;
    my $xml                = new XML::Simple();
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
    my $depth              = 0;
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
	    # Determine function return value name.
	    my $returnValueLabel;
	    if ( $node->{'parent'}->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/ ) {
		$returnValueLabel = $1;
	    } else {
		$returnValueLabel = $node->{'parent'}->{'name'}; 
	    }
	    # Determine the parameter name.
	    my $parameterName = exists($node->{'directive'}->{'parameterName'}) ? $node->{'directive'}->{'parameterName'} : $node->{'directive'}->{'class'}."Method";
	    my $defaultName   = $parameterName eq $node->{'directive'}->{'class'}."Method";
	    die("objects with defaults must have explicit parameter names")
		if ( exists($node->{'directive'}->{'default'}) && ! exists($node->{'directive'}->{'parameterName'}) );
	    # If including debugging information determine what will own the built object.
	    my $debugMessage;
	    if ( $debugging ) {
		my $ownerName;
		my $ownerLoc;
		if ( $node->{'directive'}->{'name'} =~ m/^(.+)\%[a-zA-Z0-9_]+$/ ) {
		    $ownerName = $1;
		    if ( lc($ownerName) eq lc($returnValueLabel) ) {
			$ownerLoc = "debugStackGet()";
		    } else {
			$ownerLoc = "loc(".$ownerName.")";
		    }
		} else {
		    $ownerName = $node->{'parent'}->{'name'};
		    $ownerLoc  = "'code:unknown'";
		    my $nodeParent = $node->{'parent'};
		    while ( $nodeParent ) {
			my $nodeChild = $nodeParent->{'firstChild'};
			while ( $nodeChild ) {
			    if ( $nodeChild->{'type'} eq "declaration" ) {
				foreach my $declaration ( @{$nodeChild->{'declarations'}} ) {
				    if ( grep {lc($_) eq lc($node->{'directive'}->{'name'})} @{$declaration->{'variables'}} ) {
					$ownerLoc = "'".$nodeParent->{'type'}.":".$nodeParent->{'name'}."'";
				    }
				}
			    }
			    $nodeChild = $nodeChild->{'sibling'};
			}
			$nodeParent = $nodeParent->{'parent'};
		    }
		}
		$debugMessage = "call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$node->{'directive'}->{'class'}." : ".$ownerName." : ')//".$ownerLoc."//' : '//loc(".$node->{'directive'}->{'name'}.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).")\n";
	    } else {
		$debugMessage = "";
	    }
	    # Begin code construction.
	    my $parametersDefaultRequired = 0;
	    my $builderCode;
	    $builderCode .= "   ! Determine where to build+store or point to the required object....\n";
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
		    $builderCode               .= "   do while (.not.parametersCurrent%isPresent('".$parameterName."').and.associated(parametersCurrent%parent))\n";
		    $builderCode               .= "      parametersCurrent => parametersCurrent%parent\n";
		    $builderCode               .= "   end do\n";
		    $builderCode               .=  "   if (.not.parametersCurrent%isPresent('".$parameterName."')) call Galacticus_Error_Report('[".$parameterName."] object is undefined'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		}
	    } else {	    
		$builderCode .= "   do while (.not.parametersCurrent%isPresent('".$parameterName."').and.associated(parametersCurrent%parent))\n";
		$builderCode .= "      parametersCurrent => parametersCurrent%parent\n";
		$builderCode .= "   end do\n";
	    }
	    $builderCode .= "   if (parametersCurrent%isPresent('".$parameterName."').and.(.not.".$node->{'directive'}->{'source'}."%isGlobal().or.associated(parametersCurrent%parent))) then\n"
		if ( $defaultName );
	    # Handle multiple copies.
	    my $copyInstance  = "";
	    my $copyLoopOpen  = "";
	    my $copyLoopClose = "";
	    if ( exists($node->{'directive'}->{'copy'}) ) {
		if ( $node->{'directive'}->{'copy'} =~ m/^([a-zA-Z0-9_]+)=/ ) {
		    $copyInstance  = ",copyInstance=".$1;
		    $copyLoopOpen  = "      do ".$node->{'directive'}->{'copy'}."\n";
		    $copyLoopClose = "      end do\n";
		} else {
		    $copyInstance = ",copyInstance=".$node->{'directive'}->{'copy'};
		}
	    }
	    $builderCode .= $copyLoopOpen;
	    $builderCode .= "      ! Object should belong to the parameter node. Get the node and test whether the object has already been created in it.\n";
	    $builderCode .= "      parameterNode => parametersCurrent%node('".$parameterName."'".$copyInstance.")\n";
	    $builderCode .= "      if (parameterNode%objectCreated()) then\n";
	    $builderCode .= "         ! Object already exists - simply get a pointer to it. Increment the reference counter as this is a new reference to an existing object.\n";
	    $builderCode .= "         genericObject => parameterNode%objectGet()\n";
	    $builderCode .= "         select type (genericObject)\n";
	    $builderCode .= "         class is (".$node->{'directive'}->{'class'}."Class)\n";
	    $builderCode .= "            ".$node->{'directive'}->{'name'}." => genericObject\n";
	    $builderCode .= "            call ".$node->{'directive'}->{'name'}."%referenceCountIncrement()\n";
	    $builderCode .= $debugMessage;
	    $builderCode .= "         class default\n";
	    $builderCode .= "            call Galacticus_Error_Report('parameter-stored object is not of [$node->{'directive'}->{'class'}] class'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
	    $builderCode .= "         end select\n";
	    $builderCode .= "      else\n";
	    $builderCode .= "         ! Object does not yet exist - build it and store in the parameter node. Increment reference counter here as this is a newly constructed object.\n";
	    $builderCode .= "         ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."(parametersCurrent".$copyInstance.(exists($node->{'directive'}->{'parameterName'}) ? ",parameterName='".$parameterName."'" : "").")\n";
	    $builderCode .= "            call ".$node->{'directive'}->{'name'}."%referenceCountIncrement()\n";
	    $builderCode .= $debugMessage;
	    $builderCode .= "         call parameterNode%objectSet(".$node->{'directive'}->{'name'}.")\n";
	    $builderCode .= "         call ".$node->{'directive'}->{'name'}."%autoHook()\n"
		if ( 
		    grep 
		    {exists($_->{'autoHook'}) && $_->{'autoHook'} eq "yes"} 
		    map
		    {&Galacticus::Build::Directives::Extract_Directives($_,$node->{'directive'}->{'class'})}
		    &List::ExtraUtils::as_array($directiveLocations->{$node->{'directive'}->{'class'}}->{'file'})
		);
	    $builderCode .= "      end if\n";
	    $builderCode .= $copyLoopClose;
	    if ( $defaultName ) {
		$builderCode .= "   else if (".$node->{'directive'}->{'source'}."%isGlobal()) then\n";
		$builderCode .= "      ! This is the global parameter set - so we can use the default object of this class. Increment the reference counter as this is a new reference to an existing object.\n";
		$builderCode .= $copyLoopOpen;
		$builderCode .= "      ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."()\n";
		$builderCode .= $debugMessage;
		$builderCode .= "      call ".$node->{'directive'}->{'name'}."%referenceCountIncrement()\n";
		$builderCode .= $copyLoopClose;
		$builderCode .= "   else\n";
		$builderCode .= "      ! Object is not explicitly defined, and this is not the global parameters object. Cause a default object of the class to be added to the parameters. Increment the reference count here as this is a new object.\n";
		$builderCode .= $copyLoopOpen;
		$builderCode .= "      ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."(parametersCurrent)\n";
		$builderCode .= "      call ".$node->{'directive'}->{'name'}."%referenceCountIncrement()\n";
		$builderCode .= $debugMessage;
		$builderCode .= $copyLoopClose;
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
		if ( $debugging ) {
		    $usesNode->{'moduleUse'}->{'Galacticus_Display'} =
		    {
			intrinsic => 0,
			all       => 1
		    };
		    $usesNode->{'moduleUse'}->{'String_Handling'   } =
		    {
			intrinsic => 0,
			all       => 1
		    };
		    $usesNode->{'moduleUse'}->{'ISO_Varying_String'} =
		    {
			intrinsic => 0,
			all       => 1
		    };
		}
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
	    # If including debugging information determine what will own the built object.
	    my $debugMessage;
	    if ( $debugging ) {
		my $ownerName;
		my $ownerLoc;
		if ( exists($node->{'directive'}->{'owner'}) ) {
		    $ownerName =     $node->{'directive'}->{'owner'}    ;
		    $ownerLoc  = "'".$node->{'directive'}->{'owner'}."'";
		} elsif ( $node->{'directive'}->{'name'} =~ m/^(.+)\%[a-zA-Z0-9_]+$/ ) {
		    $ownerName = $1;
		    $ownerLoc  = "loc(".$ownerName.")";
		} else {
		    $ownerName = $node->{'parent'}->{'name'};
		    $ownerLoc  = "'code:unknown'";
		    my $nodeParent = $node->{'parent'};
		    while ( $nodeParent ) {
			my $nodeChild = $nodeParent->{'firstChild'};
			while ( $nodeChild ) {
			    if ( $nodeChild->{'type'} eq "declaration" ) {
				foreach my $declaration ( @{$nodeChild->{'declarations'}} ) {
				    if ( grep {lc($_) eq lc($node->{'directive'}->{'name'})} @{$declaration->{'variables'}} ) {
					$ownerLoc = "'".$nodeParent->{'type'}.":".$nodeParent->{'name'}."'";
				    }
				}
			    }
			    $nodeChild = $nodeChild->{'sibling'};
			}
			$nodeParent = $nodeParent->{'parent'};
		    }
		}
		$debugMessage = "call Galacticus_Display_Message(var_str('functionClass[disown] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".($node->{'directive'}->{'name'} =~ m/([a-zA-Z0-9]+)_$/ ? $1 : "unknown")."] : ".$ownerName." : ')//".$ownerLoc."//' : '//loc(".$node->{'directive'}->{'name'}.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).")\n";
	    } else {
		$debugMessage = "";
	    }
	    # Generate source code for the object destructor. The pointer should only be deallocated if associated and
	    # finalizable. Otherwise simply nullify it.
	    my $destructorCode;
	    $destructorCode .= "if (associated(".$node->{'directive'}->{'name'}.")) then\n";
	    $destructorCode .= "   ! Decrement the reference count, and decide if this object can be destroyed.\n";
	    $destructorCode .= "   referenceCount_=".$node->{'directive'}->{'name'}."%referenceCountDecrement()\n";
            $destructorCode .= "   if (referenceCount_ == 0) then\n";
            $destructorCode .= "      ! Deallocate the pointer.\n";
 	    $destructorCode .= $debugMessage;
	    $destructorCode .= "      deallocate(".$node->{'directive'}->{'name'}.")\n";
            $destructorCode .= "   else if (referenceCount_ < 0) then\n";
            $destructorCode .= "      ! Negative counter - should not happen.\n";
	    if ( $debugging ) {
		$destructorCode .= "      call Galacticus_Display_Message(var_str('objectDestructor: negative reference counter (should abort, but will nullify for debugging) ".$node->{'directive'}->{'name'}." {loc: ')//loc(".$node->{'directive'}->{'name'}.")//'} at'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		$destructorCode .= "      nullify(".$node->{'directive'}->{'name'}.")\n";		
	    } else {
		$destructorCode .= "      call Galacticus_Error_Report('negative reference counter in object \"".$node->{'directive'}->{'name'}."\"'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
	    }
            $destructorCode .= "   else\n";
            $destructorCode .= "      ! Nullify the pointer.\n";
	    $destructorCode .= $debugMessage;
            $destructorCode .= "      nullify(".$node->{'directive'}->{'name'}.")\n";
            $destructorCode .= "   end if\n";
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
	    # Add variables needed by destructor.
	    if ( ! exists($node->{'parent'}->{'objectDestructorDeclarations'}) ) {
		my @declarations =
		    (
		     {
			 intrinsic  => "integer"            ,
			 variables  => [ "referenceCount_" ],
			 attributes => [                   ]
		     }
		    );
		&Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
		# Record that we have added the necessary declarations to the parent.
		$node->{'parent'}->{'objectDestructorDeclarations'} = 1;
	    }
	    # Ensure error reporting modules are used.
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
	    if ( $debugging ) {
		$usesNode->{'moduleUse'}->{'Galacticus_Display'} = {intrinsic => 0, all => 1};
		$usesNode->{'moduleUse'}->{'String_Handling'   } = {intrinsic => 0, all => 1};
		$usesNode->{'moduleUse'}->{'ISO_Varying_String'} = {intrinsic => 0, all => 1};
	    }
	    # Insert the modules node.
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	if ( $node->{'type'} eq "referenceCountIncrement" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the reference count increment.
	    my $incrementCode = "call ".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'object'}."%referenceCountIncrement()\n";
	    # If including debugging information generate required code.
	    my $debugMessage = "";
	    if ( $debugging ) {
		my $ownerName;
		my $ownerLoc;
		my $isResult = 0;
		if ( exists($node->{'directive'}->{'owner'}) ) {
		    $ownerName = $node->{'directive'}->{'owner'};
		    if ( exists($node->{'directive'}->{'isResult'}) && $node->{'directive'}->{'isResult'} eq "yes" ) {
			$ownerLoc = "debugStackGet()";
			$isResult = 1;
		    } else {
			$ownerLoc = "loc(".$ownerName.")";
		    }
		} else {
		    $ownerName = $node->{'parent'}->{'name'};
		    $ownerLoc  = "'code:unknown'";
		    my $nodeParent = $node->{'parent'};
		    while ( $nodeParent ) {
			my $nodeChild = $nodeParent->{'firstChild'};
			while ( $nodeChild ) {
			    if ( $nodeChild->{'type'} eq "declaration" ) {
				foreach my $declaration ( @{$nodeChild->{'declarations'}} ) {
				    if ( grep {lc($_) eq lc($node->{'directive'}->{'object'})} @{$declaration->{'variables'}} ) {
					$ownerLoc = "'".$nodeParent->{'type'}.":".$nodeParent->{'name'}."'";
				    }
				}
			    }
			    $nodeChild = $nodeChild->{'sibling'};
			}
			$nodeParent = $nodeParent->{'parent'};
		    }
		}
		$incrementCode .= "call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".($node->{'directive'}->{'object'} =~ m/([a-zA-Z0-9_]+)$/ ? $1 : "unknown")."] : ".$ownerName." : ')//".$ownerLoc."//' : '//loc(".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'object'}.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).")\n";
		my $usesNode =
		{
		    type      => "moduleUse",
		    moduleUse =>
		    {
			Galacticus_Display =>
			{
			    intrinsic => 0,
			    all       => 1
			},
			String_Handling    =>
		        {
			    intrinsic => 0,
			    all       => 1
			},
			ISO_Varying_String =>
		        {
			    intrinsic => 0,
			    all       => 1
			}
		    }
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    }
	    # Build a code node.
	    my $newNode =
	    {
		type       => "code"      ,
		content    => $incrementCode,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	if ( $node->{'type'} eq "referenceAcquire" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the reference aquisition.
	    my $acquireCode = (exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'target'}." => ".$node->{'directive'}->{'source'}."\n";
	    # If including debugging information push the target location to the debug stack.
	    if ( $debugging ) {
		my $ownerName;
		my $ownerLoc;
		if ( exists($node->{'directive'}->{'owner'}) ) {
		    $ownerName = $node->{'directive'}->{'owner'};
		    if ( exists($node->{'directive'}->{'isResult'}) && $node->{'directive'}->{'isResult'} eq "yes" ) {
			$ownerLoc = "debugStackGet()";
		    } else {
			$ownerLoc = "loc(".$node->{'directive'}->{'owner'}.")";
		    }
		} else {
		    $ownerName = $node->{'parent'}->{'name'};
		    $ownerLoc = "'code:unknown'";
		    my $nodeParent = $node->{'parent'};
		    while ( $nodeParent ) {
			my $nodeChild = $nodeParent->{'firstChild'};
			while ( $nodeChild ) {
			    if ( $nodeChild->{'type'} eq "declaration" ) {
				foreach my $declaration ( @{$nodeChild->{'declarations'}} ) {
				    if ( grep {lc($_) eq lc($node->{'directive'}->{'target'})} @{$declaration->{'variables'}} ) {
					$ownerLoc = "'".$nodeParent->{'type'}.":".$nodeParent->{'name'}."'";
				    }
				}
			    }
			    $nodeChild = $nodeChild->{'sibling'};
			}
			$nodeParent = $nodeParent->{'parent'};
		    }
		}		
		$acquireCode = "call debugStackPush(".$ownerLoc.")\n".$acquireCode."call debugStackPop()\n"
		    unless ( $ownerLoc eq "debugStackGet()" );
		$acquireCode .= "call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".($node->{'directive'}->{'target'} =~ m/([a-zA-Z0-9_]+)$/ ? $1 : "unknown")."] : ".$ownerName." : ')//".$ownerLoc."//' : '//loc(".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'target'}.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).")\n";
	    }
	    $acquireCode .= "call ".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'target'}."%referenceCountIncrement()\n";
	    # Build a code node.
	    my $newNode =
	    {
		type       => "code"      ,
		content    => $acquireCode,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	if ( $node->{'type'} eq "referenceConstruct" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the reference aquisition.
	    my $constructCode  = (exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'object'}."=".$node->{'directive'}->{'constructor'}."\n";
	    $constructCode    .= "call ".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'object'}."\%referenceCountIncrement()\n";
 	    # If including debugging information push the target location to the debug stack.
	    if ( $debugging ) {
		my $objectLoc = "loc(".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'object'}.")";	
		$constructCode  = "call debugStackPush(".$objectLoc.")\n".$constructCode."call debugStackPop()\n";
		my $ownerName;
		my $ownerLoc;
		my $isResult = 0;
		if ( exists($node->{'directive'}->{'isResult'}) ) {
		    if ( exists($node->{'directive'}->{'owner'}) ) {
			$ownerName = $node->{'directive'}->{'owner'};
		    } else {
			if ( $node->{'parent'}->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/ ) {
			    $ownerName = $1;
			} else {
			    $ownerName = $node->{'parent'}->{'name'}; 
			}
		    }
		    $ownerLoc = "debugStackGet()";
		    $isResult = 1;
		} elsif ( exists($node->{'directive'}->{'owner'}) ) {
		    $ownerName = $node->{'directive'}->{'owner'};
		    $ownerLoc = "loc(".$ownerName.")";
		} else {
		    $ownerName = $node->{'parent'}->{'name'};
		    if ( exists($node->{'directive'}->{'ownerLoc'}) ) {
			$ownerLoc = "'".$node->{'directive'}->{'ownerLoc'}."'";
		    } else {
			$ownerLoc  = "'code:unknown'";
			my $nodeParent = $node->{'parent'};
			while ( $nodeParent ) {
			    my $nodeChild = $nodeParent->{'firstChild'};
			    while ( $nodeChild ) {
				if ( $nodeChild->{'type'} eq "declaration" ) {
				    foreach my $declaration ( @{$nodeChild->{'declarations'}} ) {
					if ( grep {lc($_) eq lc($node->{'directive'}->{'object'})} @{$declaration->{'variables'}} ) {
					    $ownerLoc = "'".$nodeParent->{'type'}.":".$nodeParent->{'name'}."'";
					}
				    }
				}
				$nodeChild = $nodeChild->{'sibling'};
			    }
			    $nodeParent = $nodeParent->{'parent'};
			}
		    }
		}
		$constructCode .= "call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".($node->{'directive'}->{'object'} =~ m/([a-zA-Z0-9_]+)$/ ? $1 : "unknown")."] : ".$ownerName." : ')//".$ownerLoc."//' : '//loc(".(exists($node->{'directive'}->{'owner'}) ? $node->{'directive'}->{'owner'}."%" : "").$node->{'directive'}->{'object'}.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).")\n";
		my $usesNode =
		{
		    type      => "moduleUse",
		    moduleUse =>
		    {
			Galacticus_Display =>
			{
			    intrinsic => 0,
			    all       => 1
			},
			String_Handling    =>
		        {
			    intrinsic => 0,
			    all       => 1
			},
			ISO_Varying_String =>
		        {
			    intrinsic => 0,
			    all       => 1
			},
			Function_Classes   =>
			{
			    intrinsic => 0,
			    all       => 1
			}
		    }
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
		
	    }
	    # Build a code node.
	    my $newNode =
	    {
		type       => "code"      ,
		content    => $constructCode,
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
