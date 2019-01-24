# Contains a Perl module which implements processing of stateStorable directives.

package Galacticus::Build::SourceTree::Process::StateStorable;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'stateStorable'} = \&Process_StateStorable;

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
	if ( $node->{'type'} eq "stateStorable" ) {
	    # Assert that our parent is a module (for now).
	    die("Process_StateStorale: parent node must be a module")
		unless ( $node->{'parent'}->{'type'} eq "module" );
	    $moduleNode = $node->{'parent'};
	    # Extract the directive.
	    push(@directiveNodes,$node);	    
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
    # Specify types that are storable via direct transfer.
    my @transferrables =
	(
	 "fgsl_interp_type"
	);
    # Process any stateSotrable directives.
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
subroutine {$className}StateStore(self,stateFile,fgslStateFile,storeIdentifier)
 !% Store the state of this object to file.
 use, intrinsic :: ISO_C_Binding
 use            :: FGSL              , only : fgsl_file
 use            :: Galacticus_Display
 implicit none
 class    ({$className}), intent(inout)              :: self
 integer                , intent(in   )              :: stateFile
 type     (fgsl_file   ), intent(in   )              :: fgslStateFile
 logical                , intent(in   ), optional    :: storeIdentifier
 character(len=16      )                             :: label
 integer  (kind=1      ), dimension(:) , allocatable :: transferred
 integer  (c_size_t    )                             :: transferredSize
CODE
	my $inputCodeOpener  = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine {$className}StateRestore(self,stateFile,fgslStateFile)
 !% Store the state of this object to file.
 use, intrinsic :: ISO_C_Binding
 use            :: FGSL              , only : fgsl_file
 use            :: Galacticus_Display
 implicit none
 class  ({$className}), intent(inout)               :: self
 integer              , intent(in   )               :: stateFile
 type   (fgsl_file   ), intent(in   )               :: fgslStateFile
 integer(c_size_t    ), allocatable  , dimension(:) :: storedShape
 logical                                            :: wasAllocated
 integer(kind=1      ), dimension(:) , allocatable  :: transferred
 integer(c_size_t    )                              :: transferredSize
CODE
	# Close function.
	my $outputCodeCloser = fill_in_string(<<'CODE', PACKAGE => 'code');
 end select
 call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
 return
end subroutine {$className}StateStore
CODE
	my $inputCodeCloser  = fill_in_string(<<'CODE', PACKAGE => 'code');
 end select
 call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
 return
end subroutine {$className}StateRestore
CODE
	my $outputCode = fill_in_string(<<'CODE', PACKAGE => 'code');
 call Galacticus_Display_Indent('storing state for "{$className}"',verbosity=verbosityWorking)
 select type (self)
CODE
	my $inputCode  = fill_in_string(<<'CODE', PACKAGE => 'code');
 call Galacticus_Display_Indent('restoring state for "{$className}"',verbosity=verbosityWorking)
 select type (self)
CODE
	my @outputUnusedVariables;
	my @inputUnusedVariables;
	my $fgslStateFileUsed = 0;
	my $labelUsed         = 0;
	my $transferUsed      = 0;
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
		$outputCode .= "  if (present(storeIdentifier).and.storeIdentifier) write (stateFile) ".$classIdentifier."\n";
		# Search the class node for declarations.
		my @staticVariables;
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
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $labelUsed   = 1;
					    $outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    if ( $declaration->{'intrinsic'} eq "class" ) {
						$outputCode .= "  select type (c__ => self%".$variableName.")\n";
						$outputCode .= "  class is (".$declaration->{'type'}.")\n";
						$outputCode .= "   write (label,'(i16)') sizeof(c__)\n";
						$outputCode .= "  end select\n";
					    } else {
						$outputCode .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
					    }
					    $outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= " end if\n";
					    $inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $inputCode  .= " call ".$type."ClassRestore(self%".$variableName.")\n"
						if ( $declaration->{'intrinsic'} eq "class" );
					    $inputCode  .= " call self%".$variableName."%stateRestore(stateFile,fgslStateFile)\n";
					    $outputCode .= " call self%".$variableName."%stateStore  (stateFile,fgslStateFile,storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.").")\n";
					}
				    } elsif (
					$declaration->{'intrinsic'} eq "type"
					&&
					(! grep {$_ eq "dimension"} @{$declaration->{'attributes'}})	
					&&
					(  grep {$_ eq $type      } @transferrables                )
					) {
					# Type storable via direct transfer.
					$transferUsed = 1;
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $outputCode .= "transferredSize=sizeof(self%".$variableName.")\n";
					    $outputCode .= "if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    $outputCode .= " write (label,'(i16)') transferredSize\n";
					    $outputCode .= " call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= "end if\n";
					    $outputCode .= "allocate(transferred(transferredSize))\n";
					    $outputCode .= "transferred=transfer(self%".$variableName.",transferred)\n";
					    $outputCode .= "write (stateFile) transferredSize\n";
					    $outputCode .= "write (stateFile) transferred\n";
					    $outputCode .= "deallocate(transferred)\n";
					    $inputCode  .= "call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $inputCode  .= "read (stateFile) transferredSize\n";
					    $inputCode  .= "allocate(transferred(transferredSize))\n";
					    $inputCode  .= "read (stateFile) transferred\n";
					    $inputCode  .= "self%".$variableName."=transfer(transferred,self%".$variableName.")\n";
					    $inputCode  .= "deallocate(transferred)\n";
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
					    $labelUsed   = 1;
					    $outputCode .= "  if (allocated(self%".$variableName.")) then\n";
					    $outputCode .= "   if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    $outputCode .= "    write (label,'(i16)') sizeof(self%".$variableName.")\n";
					    $outputCode .= "    call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= "   end if\n";
					    $outputCode .= "   write (stateFile) .true.\n";
					    $outputCode .= "   write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n";
					    $outputCode .= "   write (stateFile) self%".$variableName."\n";
					    $outputCode .= "  else\n";
					    $outputCode .= "   write (stateFile) .false.\n";
					    $outputCode .= "  end if\n";
		    			    $inputCode  .= " read (stateFile) wasAllocated\n";
					    $inputCode  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
					    $inputCode  .= " if (wasAllocated) then\n";
					    $inputCode  .= "  call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $inputCode  .= "  allocate(storedShape(".$rank."))\n";
		    			    $inputCode  .= "  read (stateFile) storedShape\n";
					    $inputCode  .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
		    			    $inputCode  .= "  deallocate(storedShape)\n";
		    			    $inputCode  .= "  read (stateFile) self%".$variableName."\n";
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
		    $parentClassName = $classes{$parentClassName}->{'extends'};
		}
		foreach ( @staticVariables ) {
		    $labelUsed   = 1;
		    $outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
		    $outputCode .= "  write (label,'(i16)') sizeof(self%".$_.")\n";
		    $outputCode .= "  call Galacticus_Display_Message('storing \"".$_."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $outputCode .= " end if\n";
		}
		foreach ( @staticVariables ) {
		    $inputCode .= " call Galacticus_Display_Message('restoring \"".$_."\"',verbosity=verbosityWorking)\n";
		}
		$outputCode .= "  write (stateFile) ".join(", &\n  & ",map {"self%".$_} @staticVariables)."\n"
		    if ( scalar(@staticVariables) > 0 );
		$inputCode  .= "  read  (stateFile) ".join(", &\n  & ",map {"self%".$_} @staticVariables)."\n"
		    if ( scalar(@staticVariables) > 0 );
	    }
	}
	# Insert type-bindings.
	my $content = fill_in_string(<<'CODE', PACKAGE => 'code');
    !@  <objectMethods>
    !@   <object>{$className}</object>
    !@   <objectMethod>
    !@     <method>stateStore</method>
    !@     <type>void</type>
    !@     <arguments>\textcolor\{red\}\{\textless integer\textgreater\} stateFile\argin,\textcolor\{red\}\{\textless type((fgsl\_file))\textgreater\} fgslStateFile\argin,\textcolor\{red\}\{\textless integer\textgreater\} stateOperationID\argin</arguments>
    !@     <description>Store the state of this object to file.</description>
    !@   </objectMethod>
    !@   <objectMethod>
    !@     <method>stateRestore</method>
    !@     <type>void</type>
    !@     <arguments>\textcolor\{red\}\{\textless integer\textgreater\} stateFile\argin,\textcolor\{red\}\{\textless type((fgsl\_file))\textgreater\} fgslStateFile\argin,\textcolor\{red\}\{\textless integer\textgreater\} stateOperationID\argin</arguments>
    !@     <description>Restore the state of this object from file.</description>
    !@   </objectMethod>
    !@  </objectMethods>
    procedure :: stateStore   => {$className}StateStore
    procedure :: stateRestore => {$className}StateRestore
CODE
	my $binding =
	    [
	     {
		 type       => "code" ,
		 content    => $content,
		 firstChild => undef(),
		 sibling    => undef()
	     }
	    ];
	&Galacticus::Build::SourceTree::InsertPostContains($classes{$directive->{'class'}}->{'node'},$binding);
	# Record unused variables.
	push(@outputUnusedVariables,"fgslStateFile"                    )
	    unless ( $fgslStateFileUsed );
	push(@inputUnusedVariables ,"fgslStateFile"                    )
	    unless ( $fgslStateFileUsed );
	push(@outputUnusedVariables,"label"                            )
	    unless ( $labelUsed         );
	push(@outputUnusedVariables,"transferred"  , "transferredSize" )
	    unless ( $transferUsed      );
	push(@inputUnusedVariables ,"transferred"  , "transferredSize" )
	    unless ( $transferUsed      );
	# Join code fragments.
	my $functionCode = 
	    $outputCodeOpener.
	    (@outputUnusedVariables ? " !GCC\$ attributes unused :: ".join(", ",@outputUnusedVariables)."\n" : "").
	    $outputCode      .
	    $outputCodeCloser.
	    "\n"             .
	    $inputCodeOpener .
	    (@inputUnusedVariables ? " !GCC\$ attributes unused :: ".join(", ",@inputUnusedVariables )."\n" : "").
	    $inputCode       .
	    $inputCodeCloser .
	    "\n"             ;
	# Build the class restore functions.
	foreach $code::parentClassName ( sort(keys(%classes)) ) {
	    my $classRestoreCode = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine {$parentClassName}ClassRestore(self,stateFile)
 !% Restore the class of this object from file.
 use Galacticus_Error
 implicit none
 class  ({$parentClassName}), intent(inout), allocatable  :: self
 integer                    , intent(in   )               :: stateFile
 integer                                                  :: classIdentifier
 read (stateFile) classIdentifier
 select case (classIdentifier)
CODE
	    foreach my $childClassName( sort(keys(%classIdentifiers)) ) {
		my $parentClassName = $childClassName;
		while ( defined($parentClassName) ) {
		    if ( $parentClassName eq $code::parentClassName ) {
			# This class extends our parent class, so include it in this function.
			$classRestoreCode .= " case (".$classIdentifiers{$childClassName}.")\n";
			$classRestoreCode .= "  if (allocated(self)) deallocate(self)\n";
			$classRestoreCode .= "  allocate(".$childClassName." :: self)\n";
		    }
		    $parentClassName = $classes{$parentClassName}->{'extends'};
		}
	    }
	    $classRestoreCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
 case default
  call Galacticus_Error_Report('serialized object is of incorrect class')
 end select
 return
end subroutine {$parentClassName}ClassRestore
CODE
	    $functionCode .= $classRestoreCode."\n";
	    # Set visibility of class restore function.
	    &Galacticus::Build::SourceTree::SetVisibility($moduleNode,$code::parentClassName."ClassRestore","public");
	}
	# Insert code.
	my $treeTmp = &Galacticus::Build::SourceTree::ParseCode($functionCode,'null');
	&Galacticus::Build::SourceTree::ProcessTree($treeTmp);
	my $postContains =
	    [
	     {
		 type       => "code" ,
		 content    => &Galacticus::Build::SourceTree::Serialize($treeTmp),
		 firstChild => undef(),
		 sibling    => undef(),
		 parent     => undef(),
		 source     => "Galacticus::Build::SourceTree::Process::StateStorable::Process_StateStorable()",
		 line       => 1
	     }
	    ];
	&Galacticus::Build::SourceTree::InsertPostContains($directiveNode->{'parent'},$postContains);
    }
}

1;
