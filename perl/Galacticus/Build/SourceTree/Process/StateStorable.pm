# Contains a Perl module which implements processing of stateStorable directives.

package Galacticus::Build::SourceTree::Process::StateStorable;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Template 'fill_in_string';
use Galacticus::Build::SourceTree::Process::Utils qw(performIO);

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'stateStorable'} = \&Process_StateStorable;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'stateStorable'} = [ "generics" ];

sub Process_StateStorable {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Initialize the directives and classes.
    my @directiveNodes;
    my %classes;
    # Initialize state storables database.
    my $stateStorables;
    # Walk the tree.
    my $node       = $tree;
    my $moduleNode        ;
    my $depth      = 0    ;
    while ( $node ) {
	# Capture stateStorable directives.
	if ( $node->{'type'} eq "stateStorable" && ! $node->{'directive'}->{'processed'} ) {
	    # Get the module node if contained within a module.
	    $moduleNode = $node->{'parent'}
	        if ( $node->{'parent'}->{'type'} eq "module" );	    
	    # Extract the directive.
	    push(@directiveNodes,$node);	    
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} = 1;
	    # Get state storables database if we do not have it.
	    $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml")
		unless ( $stateStorables );	    
	}
	# Capture derived type definitions.
	if ( $node->{'type'} eq "type" ) {
	    # Parse class openers to find dependencies.
	    if ( $node->{'opener'} =~ m/^\s*type\s*(,\s*(abstract)\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*([a-z0-9_]+)\s*$/i ) {
		my $type     = $5;
		my $extends  = $3;
		my $abstract = defined($2);
		$classes{$type} =
		{
		     node     => $node   ,
		     extends  => $extends,
		     abstract => $abstract
		};
	    } else {
		die("Galacticus::Build::SourceTree::Process::StateStorable::Process_StateStorable: unable to parse type definition opener");
	    }
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }

    # Remove any type definitions which do not match the storable class.
    if ( scalar(@directiveNodes) > 0 ) {
	foreach my $className ( keys(%classes) ) {
	    my $matched = 0;
	    my $parentClassName = $className;
	    while ( ! $matched ) {
		$matched = 1
		    if ( grep {$_->{'directive'}->{'class'} eq $parentClassName} @directiveNodes );
		last
		    if ( $matched );
		last
		    unless ( exists($classes{$parentClassName}) && defined($classes{$parentClassName}->{'extends'}) );
		$parentClassName = $classes{$parentClassName}->{'extends'};
	    }
	    delete($classes{$className})
		unless ( $matched );
	}
    }
    
    # Record whether we need class functions.
    my $classFunctionsRequired = defined($moduleNode);
    # Process any stateStorable directives.
    foreach my $directiveNode ( @directiveNodes ) {
	my $directive = $directiveNode->{'directive'};
	# Assert that class must exist.
	die("Galacticus::Build::SourceTree::Process::StateStorable::Process_StateStorable: class '".$directive->{'class'}."' not found")
	    unless ( exists($classes{$directive->{'class'}}) );
	# Begin building state store function code.
	$code::className    = $directive->{'class'};
	my $classIdentifier = -1;
	my %classIdentifiers;
	my $outputCodeOpener = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine {$className}StateStore(self,stateFile,gslStateFile,storeIdentifier)
 !!\{
 Store the state of this object to file.
 !!\}
 use, intrinsic :: ISO_C_Binding     , only : c_size_t, c_ptr
 use            :: Display
 implicit none
 class    ({$className}), intent(inout)              :: self
 integer                , intent(in   )              :: stateFile
 type     (c_ptr       ), intent(in   )              :: gslStateFile
 logical                , intent(in   ), optional    :: storeIdentifier
 character(len=16      )                             :: label
CODE
	my $inputCodeOpener  = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine {$className}StateRestore(self,stateFile,gslStateFile)
 !!\{
 Store the state of this object to file.
 !!\}
 use, intrinsic :: ISO_C_Binding     , only : c_size_t, c_ptr
 use            :: Display
 implicit none
 class  ({$className}), intent(inout)               :: self
 integer              , intent(in   )               :: stateFile
 type   (c_ptr       ), intent(in   )               :: gslStateFile
CODE
	my $storedShapeRequired  = 0;
	my $wasAllocatedRequired = 0;
	# Close function.
	my $outputCodeCloser = fill_in_string(<<'CODE', PACKAGE => 'code');
 end select
 call displayUnindent('done',verbosity=verbosityLevelWorking)
 return
end subroutine {$className}StateStore
CODE
	my $inputCodeCloser  = fill_in_string(<<'CODE', PACKAGE => 'code');
 end select
 call displayUnindent('done',verbosity=verbosityLevelWorking)
 return
end subroutine {$className}StateRestore
CODE
	my $outputCode = fill_in_string(<<'CODE', PACKAGE => 'code');
 call displayIndent('storing state for "{$className}"',verbosity=verbosityLevelWorking)
 select type (self)
CODE
	my $inputCode  = fill_in_string(<<'CODE', PACKAGE => 'code');
 call displayIndent('restoring state for "{$className}"',verbosity=verbosityLevelWorking)
 select type (self)
CODE
	my @outputUnusedVariables;
	my @inputUnusedVariables;
	my $gslStateFileUsed = 0;
	my $labelUsed        = 0;
	my $transferUsed     = 0;
	my $rankMaximum      = 0;
	# Scan all known classes, finding all which derive from the base class.
	foreach my $className ( sort(keys(%classes)) ) {
	    my $matches = 0;
	    my $parentClassName = $className;
	    while ( defined($parentClassName) ) {
		if ( $parentClassName eq $directive->{'class'} ) {
		    $matches = 1;
		    last;
		}
		$parentClassName = $classes{$parentClassName}->{'extends'};
	    }
	    if ( $matches && ! $classes{$className}->{'abstract'} ) {
		# This class is derived from the class of interest, and is not abstract, so should be stored.
		++$classIdentifier;
		$classIdentifiers{$className} = $classIdentifier;
		$outputCode .= " type is (".$className.")\n";
		$inputCode  .= " type is (".$className.")\n";
		$outputCode .= &performIO("  if (present(storeIdentifier).and.storeIdentifier) write (stateFile) ".$classIdentifier."\n");
		# Search the class node for declarations.
		my @staticVariables;
		my @methodCalls;
		my $parentClassName = $className;
		while ( defined($parentClassName) ) {
		    # Find any variables to be excluded from state store/restore.
		    my @excludes = exists($directive->{$parentClassName}->{'exclude'}->{'variables'}) ? split(/\s*,\s*/,$directive->{$parentClassName}->{'exclude'}->{'variables'}) : ();
		    my $classNode = $classes{$parentClassName}->{'node'}->{'firstChild'};
		    while ( $classNode ) {
			if ( $classNode->{'type'} eq "declaration" ) {
			    foreach my $declaration ( @{$classNode->{'declarations'}} ) {
				# Identify variable type.
				if ( $declaration->{'intrinsic'} eq "procedure" ) {
				    # Type-bound procedure - nothing to do.
				} elsif ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) {
				    (my $type = $declaration->{'type'}) =~ s/\s//g;
				    if (
					(  grep {$_->{'type'} eq $type    } @{$stateStorables->{'stateStorables'}})
					&&
			                (! grep {$_           eq "pointer"} @{$declaration   ->{'attributes'     }})				       
					){
					# This is a non-pointer object which is explicitly stateStorable.
					my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
					my $rank                = ($dimensionDeclarator =~ tr/://);
					$rankMaximum            = $rank
					    if ( $rank > $rankMaximum );
					my $allocatable         =  grep {$_ eq "allocatable"} @{$declaration->{'attributes'}};
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $labelUsed   = 1;
					    if ( $allocatable ) {
						$outputCode .= "  if (allocated(self%".$variableName.")) then\n";
						$outputCode .= &performIO("   write (stateFile) .true.\n"
							    .             "   write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n");
					    }
					    $outputCode .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
					    if ( $declaration->{'intrinsic'} eq "class" ) {
						$outputCode .= "  select type (c__ => self%".$variableName.")\n";
						$outputCode .= "  class is (".$declaration->{'type'}.")\n";
						$outputCode .= &performIO("   write (label,'(i16)') sizeof(c__)\n");
						$outputCode .= "  end select\n";
					    } else {
						$outputCode .= &performIO("   write (label,'(i16)') sizeof(self%".$variableName.")\n");
					    }
					    $outputCode .= "  call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= " end if\n";
					    if ( $rank > 0 ) {
						for(my $i=1;$i<=$rank;++$i) {
						    $outputCode .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
						}
					    }					    
					    $outputCode .= " call self%".$variableName.($rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "")."%stateStore  (stateFile,gslStateFile,storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.").")\n";
					    if ( $rank > 0 ) {
						for(my $i=$rank;$i>0;--$i) {
						    $outputCode .= (" " x $i)."end do\n";
						}
					    }					    
					    if ( $allocatable ) {
						$outputCode .= "  else\n";
						$outputCode .= &performIO("   write (stateFile) .false.\n");
						$outputCode .= "  end if\n";
					    }
					    if ( $allocatable ) {
						$wasAllocatedRequired =1;
						$inputCode  .= &performIO(" read (stateFile) wasAllocated\n");
						$inputCode  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
						$inputCode  .= " if (wasAllocated) then\n";
					    }
					    $inputCode  .= "  call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
					    if ( $allocatable ) {
						$storedShapeRequired = 1;
						$inputCode  .= "  allocate(storedShape(".$rank."))\n";
						$inputCode  .= &performIO("  read (stateFile) storedShape\n");
						$inputCode  .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
						$inputCode  .= "  deallocate(storedShape)\n";
					    }
					    if ( $declaration->{'intrinsic'} eq "class" ) {
						$inputCode  .= " call ".$type."ClassRestore(self%".$variableName.",stateFile)\n";
						$classFunctionsRequired = 1;
					    }
					    if ( $rank > 0 ) {
						for(my $i=1;$i<=$rank;++$i) {
						    $inputCode .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
						}
					    }					    
					    $inputCode  .= " call self%".$variableName.($rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "")."%stateRestore(stateFile,gslStateFile)\n";
					    if ( $rank > 0 ) {
						for(my $i=$rank;$i>0;--$i) {
						    $inputCode .= (" " x $i)."end do\n";
						}
					    }					    
					    if ( $allocatable ) {
						$inputCode  .= " end if\n";
					    }
					}
				    }
				} else {
				    # Intrinsic type.
				    if ( grep {$_ eq "pointer"} @{$declaration->{'attributes'}} ) {
					# Pointers are currently not handled.
				    } elsif ( exists($declaration->{'type'}) && defined($declaration->{'type'}) && $declaration->{'type'} =~ m/^\s*omp_lock_kind\s*/ ) {
					# Do not store OpenMP lock variables.
				    } elsif ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
					# For allocatable variables we must first store the shape so that they can be reallocated on restore.
					my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
					my $rank = ($dimensionDeclarator =~ tr/://);
					foreach my $variableName ( @{$declaration->{'variables'}} ) {
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $storedShapeRequired  = 1;
					    $wasAllocatedRequired = 1;
					    $labelUsed            = 1;
					    $outputCode .= "  if (allocated(self%".$variableName.")) then\n";
					    $outputCode .= "   if (displayVerbosity() >= verbosityLevelWorking) then\n";
					    $outputCode .= &performIO("    write (label,'(i16)') sizeof(self%".$variableName.")\n");
					    $outputCode .= "    call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= "   end if\n";
					    $outputCode .= &performIO("   write (stateFile) .true.\n"
					                .  "   write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
					                .  "   write (stateFile) self%".$variableName."\n");
					    $outputCode .= "  else\n";
					    $outputCode .= &performIO("   write (stateFile) .false.\n");
					    $outputCode .= "  end if\n";
		    			    $inputCode  .= &performIO(" read (stateFile) wasAllocated\n");
					    $inputCode  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
					    $inputCode  .= " if (wasAllocated) then\n";
					    $inputCode  .= "  call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
					    $inputCode  .= "  allocate(storedShape(".$rank."))\n";
		    			    $inputCode  .= &performIO("  read (stateFile) storedShape\n");
					    $inputCode  .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
		    			    $inputCode  .= "  deallocate(storedShape)\n";
		    			    $inputCode  .= &performIO("  read (stateFile) self%".$variableName."\n");
		    			    $inputCode  .= " end if\n";
					}
				    } else {
					# Statically-sized variable.
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    my $store = 1;
					    if ( exists($directive->{$parentClassName}) && exists($directive->{$parentClassName}->{'restoreTo'}) ) {
					        foreach ( &List::ExtraUtils::as_array($directive->{$parentClassName}->{'restoreTo'}) ) {
						    my @variables = split(/\s*,\s*/,$_->{'variables'});
						    if ( grep {lc($_) eq lc($variableName)} @variables ) {
							$store = 0;
							$inputCode .= " self%".$variableName."=".$_->{'state'}."\n";
						    }
					        }
					    }
					    push(@staticVariables,$variableName)
						if ( $store );
					}
				    }
				}
			    }
			}
			$classNode = $classNode->{'sibling'};
		    }
		    # Add any explicit method calls.
		    if ( exists($directive->{$parentClassName}) && exists($directive->{$parentClassName}->{'methodCall'}) ) {
			foreach ( &List::ExtraUtils::as_array($directive->{$parentClassName}->{'methodCall'}) ) {
			    push(@methodCalls,"  call self%".$_->{'method'}."(".(exists($_->{'arguments'}) ? $_->{'arguments'} : "").")");
			}
		    }
		    $parentClassName = $classes{$parentClassName}->{'extends'};
		}
		foreach ( @staticVariables ) {
		    $labelUsed   = 1;
		    $outputCode .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    $outputCode .= &performIO("  write (label,'(i16)') sizeof(self%".$_.")\n");
		    $outputCode .= "  call displayMessage('storing \"".$_."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $outputCode .= " end if\n";
		}
		foreach ( @staticVariables ) {
		    $inputCode .= " call displayMessage('restoring \"".$_."\"',verbosity=verbosityLevelWorking)\n";
		}
		$outputCode .= &performIO("  write (stateFile) ".join(", &\n  & ",map {"self%".$_} @staticVariables)."\n")
		    if ( scalar(@staticVariables) > 0 );
		$inputCode  .= &performIO("  read  (stateFile) ".join(", &\n  & ",map {"self%".$_} @staticVariables)."\n")
		    if ( scalar(@staticVariables) > 0 );
		$inputCode  .= join("\n",@methodCalls)."\n"
		    if ( @methodCalls );
	    }
	}
	# Insert type-bindings.
	my $content = fill_in_string(<<'CODE', PACKAGE => 'code');
    !![
    <methods>
     <method method="stateStore"   description="Store the state of this object to file."    />
     <method method="stateRestore" description="Restore the state of this object from file."/>
    </methods>
    !!]
    procedure :: stateStore   => {$className}StateStore
    procedure :: stateRestore => {$className}StateRestore
CODE
	my $bindingTree = &Galacticus::Build::SourceTree::ParseCode($content,"Galacticus::Build::SourceTree::Process::StateStorable()");
	my @bindingNodes = &Galacticus::Build::SourceTree::Children($bindingTree);
	&Galacticus::Build::SourceTree::InsertPostContains($classes{$directive->{'class'}}->{'node'},\@bindingNodes);
	# Record unused variables.
	push(@outputUnusedVariables,"gslStateFile"                    )
	    unless ( $gslStateFileUsed );
	push(@inputUnusedVariables ,"gslStateFile"                    )
	    unless ( $gslStateFileUsed );
	push(@outputUnusedVariables,"label"                           )
	    unless ( $labelUsed        );

	if ( $storedShapeRequired ) {
	    $inputCodeOpener .= fill_in_string(<<'CODE', PACKAGE => 'code');
 integer(c_size_t    ), allocatable  , dimension(:) :: storedShape
CODE
	}
	if ( $wasAllocatedRequired ) {
	    $inputCodeOpener .= fill_in_string(<<'CODE', PACKAGE => 'code');
 logical                                            :: wasAllocated
CODE
	}
	if ( $transferUsed     ) {
	    $outputCodeOpener .= fill_in_string(<<'CODE', PACKAGE => 'code');
 integer  (kind=1      ), dimension(:) , allocatable :: transferred
 integer  (c_size_t    )                             :: transferredSize
CODE
	    $inputCodeOpener .= fill_in_string(<<'CODE', PACKAGE => 'code');
 integer  (kind=1      ), dimension(:) , allocatable :: transferred
 integer  (c_size_t    )                             :: transferredSize
CODE
	}
	if ( $rankMaximum > 0  ) {
	    $outputCodeOpener .= " integer  (c_size_t    )                             :: ".join(", ",map {"i".$_} 1..$rankMaximum)."\n";
	    $inputCodeOpener  .= " integer  (c_size_t    )                             :: ".join(", ",map {"i".$_} 1..$rankMaximum)."\n";
	}
	# Join code fragments.
	my $functionCode = 
	    $outputCodeOpener.
	    (@outputUnusedVariables ? " !\$GLC attributes unused :: ".join(", ",@outputUnusedVariables)."\n" : "").
	    $outputCode      .
	    $outputCodeCloser.
	    "\n"             .
	    $inputCodeOpener .
	    (@inputUnusedVariables ? " !\$GLC attributes unused :: ".join(", ",@inputUnusedVariables )."\n" : "").
	    $inputCode       .
	    $inputCodeCloser .
	    "\n"             ;
	# Build the class restore functions.
	if ( $classFunctionsRequired ) {
	    foreach $code::parentClassName ( sort(keys(%classes)) ) {
		for(my $rank=0;$rank<=1;++$rank) {
		    $code::rankSuffix  = $rank > 0 ? $rank."D"      : "";
		    $code::storedShape = $rank > 0 ? ",storedShape" : "";
		    $code::dimensions  = $rank > 0 ? ", dimension(".join(",",map {":"} 1..$rank).")" : "";
		    my $classRestoreCode = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine {$parentClassName}ClassRestore{$rankSuffix}(self,stateFile{$storedShape})
 !!\{
 Restore the class of this object from file.
 !!\}
 use            :: Error        , only : Error_Report
 use, intrinsic :: ISO_C_Binding, only : c_size_t
 implicit none
 class  ({$parentClassName}), intent(inout), allocatable{$dimensions} :: self
 integer                    , intent(in   )               :: stateFile
 integer                                                  :: classIdentifier
CODE
		    $classRestoreCode .= "integer(c_size_t), intent(in   ), dimension(".join(",",map {":"} 1..$rank).") :: storedShape\n"
			if ( $rank > 0 );
		    $classRestoreCode .= &performIO(" read (stateFile) classIdentifier\n");
		    $classRestoreCode .= " select case (classIdentifier)\n";
		    foreach my $childClassName( sort(keys(%classIdentifiers)) ) {
			my $parentClassName = $childClassName;
			while ( defined($parentClassName) ) {
			    if ( $parentClassName eq $code::parentClassName ) {
				# This class extends our parent class, so include it in this function.
				$classRestoreCode .= " case (".$classIdentifiers{$childClassName}.")\n";
				$classRestoreCode .= "  if (allocated(self)) deallocate(self)\n";
				$classRestoreCode .= "  allocate(".$childClassName." :: self".($rank > 0 ? "(".join(",",map {"storedShape(".$_.")"} 1..$rank).")" : "").")\n";
			    }
			    $parentClassName = $classes{$parentClassName}->{'extends'};
			}
		    }
		    $classRestoreCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
 case default
  call Error_Report('serialized object is of incorrect class')
 end select
 return
end subroutine {$parentClassName}ClassRestore{$rankSuffix}
CODE
		    $functionCode .= $classRestoreCode."\n";
		    # Set visibility of class restore function.
		    &Galacticus::Build::SourceTree::SetVisibility($moduleNode,$code::parentClassName."ClassRestore".$code::rankSuffix,"public")
			if ( defined($moduleNode) );
		}
	    }
	}
	# Insert code.
	my $treeTmp = &Galacticus::Build::SourceTree::ParseCode($functionCode,"Galacticus::Build::SourceTree::Process::StateStorable()");
	&Galacticus::Build::SourceTree::ProcessTree($treeTmp);
	my @nodesTmp = &Galacticus::Build::SourceTree::Children($treeTmp);
	&Galacticus::Build::SourceTree::InsertPostContains($directiveNode->{'parent'},\@nodesTmp);
    }
}

1;
