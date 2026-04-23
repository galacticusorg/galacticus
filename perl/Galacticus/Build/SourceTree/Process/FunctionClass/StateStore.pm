# Contains a Perl module which provides state store/restore code generation utilities for functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass::StateStore;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Process::FunctionClass::Utils qw(stripVariableName declarationRank);
use Exporter 'import';

our @EXPORT_OK = qw(stateStoreExplicitFunction stateStoreVariables generateAllocatableStateStoreCode);

# Alias shared state from the parent package.
our $stateStorables;
*stateStorables = \$Galacticus::Build::SourceTree::Process::FunctionClass::stateStorables;

sub stateStoreExplicitFunction {
    # Create state store/restore instructions for objects with explicit functions.
    my $nonAbstractClass  = shift();
    my $inputCode  = "";
    my $outputCode = "";
    my %modules;
    # Handle store.
    if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'variables'}) ) {
	foreach my $explicitVariable ( split(" ",$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'variables'}) ) {
	    if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'store'  }) ) {
		$outputCode .= "if (associated(self%".$explicitVariable.")) then\n";
		$outputCode .= " write (stateFile) .true.\n";
		$outputCode .= " call ".$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'store'  }."(self%".$explicitVariable.",stateFile,gslStateFile,stateOperationID)\n";
		$outputCode .= "else\n";
		$outputCode .= " write (stateFile) .false.\n";
		$outputCode .= "end if\n";
	    }
	    if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'restore'}) ) {
		$inputCode  .= "read (stateFile) wasAssociated\n";
		$inputCode  .= "if (wasAssociated) then\n";
		$inputCode  .= " call ".$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'restore'}."(self%".$explicitVariable.",stateFile,gslStateFile,stateOperationID)\n";
		$inputCode  .= "else\n";
		$inputCode  .= " nullify(self%".$explicitVariable.")\n";
		$inputCode  .= "end if\n";
	    }
	    $modules{$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'module'}} = 1
		if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'module' }) );
	}
    }
    return ($inputCode,$outputCode,%modules);
}

sub generateAllocatableStateStoreCode {
    # Generate state store/restore code for an allocatable variable.
    my $stateStore   = shift();
    my $stateStores  = shift();
    my $variableName = shift();
    my $rank         = shift();
    my $writeExpr    = shift();  # Fortran expression to write
    my $readExpr     = shift();  # Fortran expression to read into
    $stateStores->{'allocatablesFound'}  = 1;
    $stateStores->{'dimensionalsFound'}  = 1;
    $stateStores->{'stateFileUsed'}      = 1;
    $stateStores->{'labelUsed'}          = 1;
    $stateStore->{'outputCode'}         .= " if (allocated(self%".$variableName.")) then\n";
    $stateStore->{'outputCode'}         .= "  if (displayVerbosity() >= verbosityLevelWorking) then\n";
    # <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
    #  <description>
    #   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
    #  </description>
    # </workaround>
    $stateStore->{'outputCode'}         .= "   write (label,'(i16)') 0\n";
    $stateStore->{'outputCode'}         .= "   call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
    $stateStore->{'outputCode'}         .= "  end if\n";
    $stateStore->{'outputCode'}         .= "  write (stateFile) .true.\n"
	. "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
	. "  write (stateFile) ".$writeExpr."\n";
    $stateStore->{'outputCode'}         .= " else\n";
    $stateStore->{'outputCode'}         .= "  write (stateFile) .false.\n";
    $stateStore->{'outputCode'}         .= " end if\n";
    $stateStore->{'inputCode'}          .= " read (stateFile) wasAllocated\n";
    $stateStore->{'inputCode'}          .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
    $stateStore->{'inputCode'}          .= " if (wasAllocated) then\n";
    $stateStore->{'inputCode'}          .= "  call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
    $stateStore->{'inputCode'}          .= "  allocate(storedShape(".$rank."))\n";
    $stateStore->{'inputCode'}          .= "  read (stateFile) storedShape\n";
    $stateStore->{'inputCode'}          .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
    $stateStore->{'inputCode'}          .= "  deallocate(storedShape)\n";
    $stateStore->{'inputCode'}          .= "  read (stateFile) ".$readExpr."\n";
    $stateStore->{'inputCode'}          .= " end if\n";
}

sub stateStoreVariables {
    # Generate code to store/restore variables in functionClass objects.
    my  $stateStores        = shift();
    my  $stateStore         = shift();
    my  $class              = shift();
    my  $declarations       = shift();
    my  $explicitNamesFound = shift();
    foreach my $declaration ( &List::ExtraUtils::as_array($declarations) ) {
	# Identify variable type.
	if ( $declaration->{'intrinsic'} eq "procedure" || $declaration->{'intrinsic'} eq "final" ) {
	    # Type-bound procedure - nothing to do.
	} elsif ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) {
	    # Look for pointers to functionClasses.
	    (my $type = $declaration->{'type'}) =~ s/\s//g;
	    if (
		$declaration->{'intrinsic'} eq "class"
		&&
		(grep {$_ eq "pointer"}      @{$declaration   ->{'attributes'     }} )
		&&
		(grep {$_ eq $type    } keys(%{$stateStorables->{'functionClasses'}}))
		) {
		# Pointer to a functionClass object.
		foreach ( @{$declaration->{'variables'}} ) {
		    $stateStores->{'labelUsed'} = 1;
		    my $variableName = &stripVariableName($_);
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    $stateStore->{'outputCode'} .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    $stateStore->{'outputCode'} .= "  select type (c__ => self%".$variableName.")\n";
		    $stateStore->{'outputCode'} .= "  class is (".$declaration->{'type'}.")\n";
		    # <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
		    #  <description>
		    #   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
		    #  </description>
		    # </workaround>
		    $stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
		    #$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(c__)\n";
		    $stateStore->{'outputCode'} .= "  end select\n";
		    $stateStore->{'outputCode'} .= "  call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $stateStore->{'outputCode'} .= " end if\n";
		    $stateStore->{'inputCode'}  .= " call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
		    $stateStore->{'outputCode'} .= " call self%".$variableName."%stateStore  (stateFile,gslStateFile,stateOperationID)\n";
		    $stateStore->{'inputCode'}  .= " call self%".$variableName."%stateRestore(stateFile,gslStateFile,stateOperationID)\n";
		    $stateStores->{'stateFileUsed'}    = 1;
		    $stateStores->{'gslStateFileUsed'} = 1;
		}
	    } elsif (
		$declaration->{'intrinsic'} eq "type"
		&&
		$declaration->{'type'     } =~ m/^enumeration[a-z0-9_]+type/i
		) {
		# Enumeration.
		if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
		    # For allocatable variables we must first store the shape so that they can be reallocated on restore.
		    my $rank = &declarationRank($declaration);
		    foreach my $variableName ( @{$declaration->{'variables'}} ) {
			next
			    if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
			&generateAllocatableStateStoreCode($stateStore,$stateStores,$variableName,$rank,
			    "self%".$variableName."%ID","self%".$variableName."%ID");
		    }
		} else {
		    $stateStore->{'outputCode'} .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    foreach ( @{$declaration->{'variableNames'}} ) {
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(".$_.")\n";
			$stateStore->{'outputCode'} .= "  call displayMessage('storing \"".$_."\" with size '//trim(adjustl(label))//' bytes')\n";
		    }
		    $stateStore->{'outputCode'} .= " end if\n";
		    $stateStore->{'outputCode'} .= "  write (stateFile) ".join(",",map {"self%".$_."%ID"} @{$declaration->{'variableNames'}})."\n";
		    $stateStore->{'inputCode'}  .= "  read  (stateFile) ".join(",",map {"self%".$_."%ID"} @{$declaration->{'variableNames'}})."\n";
		}
	    } elsif (
		(  grep {$_->{'type'} eq $type    } &List::ExtraUtils::as_array($stateStorables->{'stateStorables'        }))
		||
		(  grep {$_           eq $type    } &List::ExtraUtils::as_array($stateStorables->{'functionClassInstances'}))
		){
		# This is a non-pointer object which is explicitly stateStorable.
		# Get presence of pointer attribute.
		my $isPointer = grep {$_ eq "pointer"} @{$declaration->{'attributes'}};
		# Get list of named pointer variables that are allowed.
		my @explicits = (defined($class) && exists($class->{'stateStorable'}->{'functionClass'}->{'variables'})) ? split(/\s*,\s*/,$class->{'stateStorable'}->{'functionClass'}->{'variables'}) : ();
		my $isFunctionClass = grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}};
		# Construct code to output.
		foreach ( @{$declaration->{'variables'}} ) {
		    my $variableName = &stripVariableName($_);
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    my $isExplicit = grep {lc($_) eq lc($variableName)} @explicits;
		    next
			unless ( (! $isPointer) || $isExplicit );
		    push(@{$explicitNamesFound},lc($variableName))
			if ( $isExplicit );
		    my $rank = &declarationRank($declaration);
		    $stateStores->{'rankMaximum'} = $rank
			if ( $rank > $stateStores->{'rankMaximum'} );
		    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
			# For allocatable variables we must first store the shape so that they can be reallocated on restore.
			$stateStores->{'allocatablesFound'}  = 1;
			$stateStores->{'dimensionalsFound'}  = 1
			    if ( $rank > 0 );
			$stateStore->{'outputCode'} .= " if (allocated(self%".$variableName.")) then\n";
			$stateStore->{'outputCode'} .= "  write (stateFile) .true.\n";
			$stateStore->{'outputCode'} .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
			    if ( $rank > 0 );
			$stateStore->{'inputCode'}  .= " read (stateFile) wasAllocated\n";
			$stateStore->{'inputCode'}  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
			$stateStore->{'inputCode'}  .= " if (wasAllocated) then\n";
			if ( $rank > 0 ) {
			    $stateStore->{'inputCode'}  .= "  allocate(storedShape(".$rank."))\n";
			    $stateStore->{'inputCode'}  .= "  read (stateFile) storedShape\n";
			}
			if ( $declaration->{'intrinsic'} eq "class" ) {
			    (my $storable) = grep {$_->{'type'} eq $type} @{$stateStorables->{'stateStorables'}};
			    my $functionName = $type."ClassRestore".($rank > 0 ? $rank."D" : "");
			    $stateStores->{'stateRestoreModules'}->{$storable->{'module'}.",only:".$functionName} = 1;
			    $stateStore->{'inputCode'}  .= "  call ".$functionName."(self%".$variableName.",stateFile".($rank > 0 ? ",storedShape" : "").")\n";
			} else {
			    $stateStore->{'inputCode'}  .= "  allocate(self%".$variableName.($rank > 0 ? "(".join(",",map {"storedShape(".$_.")"} 1..$rank).")" : "").")\n";
			}
			if ( $rank > 0 ) {
			    $stateStore->{'inputCode'}  .= "  deallocate(storedShape)\n";
			}
		    }
		    for(my $i=1;$i<=$rank;++$i) {
			$stateStore->{'outputCode'} .= (" " x $i)."do i".$i."=lbound(self%".$variableName.",dim=".$i."),ubound(self%".$variableName.",dim=".$i.")\n";
			$stateStore->{'inputCode'}  .= (" " x $i)."do i".$i."=lbound(self%".$variableName.",dim=".$i."),ubound(self%".$variableName.",dim=".$i.")\n";
		    }
		    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
		    $stateStores->{'labelUsed'}   = 1;
		    $stateStore->{'outputCode'} .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    if ( $declaration->{'intrinsic'} eq "class" ) {
			$stateStore->{'outputCode'} .= "  select type (c__ => self%".$variableName.$arrayElement.")\n";
			$stateStore->{'outputCode'} .= "  class is (".$declaration->{'type'}.")\n";
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(c__)\n";
			$stateStore->{'outputCode'} .= "  end select\n";
		    } else {
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
		    }
		    $stateStore->{'outputCode'} .= "  call displayMessage('storing \"".$variableName.$arrayElement."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $stateStore->{'outputCode'} .= " end if\n";
		    $stateStore->{'inputCode'}  .= " call displayMessage('restoring \"".$variableName.$arrayElement."\"',verbosity=verbosityLevelWorking)\n";
		    $stateStore->{'inputCode'}  .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateRestore(stateFile,gslStateFile".($isFunctionClass ? ",stateOperationID" : "").")\n";
		    $stateStore->{'outputCode'} .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateStore  (stateFile,gslStateFile".($isFunctionClass ? ",stateOperationID" : ",storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.")).")\n";
		    for(my $i=1;$i<=$rank;++$i) {
			$stateStore->{'outputCode'} .= (" " x ($rank+1-$i))."end do\n";
			$stateStore->{'inputCode'}  .= (" " x ($rank+1-$i))."end do\n";
		    }
		    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
			$stateStore->{'inputCode'}  .= " end if\n";
			$stateStore->{'outputCode'} .= " else\n";
			$stateStore->{'outputCode'} .= "  write (stateFile) .false.\n";
			$stateStore->{'outputCode'} .= " end if\n";
		    }
		    $stateStores->{'stateFileUsed'}    = 1;
		    $stateStores->{'gslStateFileUsed'} = 1;
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
		my $rank = &declarationRank($declaration);
		foreach my $variableName ( @{$declaration->{'variables'}} ) {
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    &generateAllocatableStateStoreCode($stateStore,$stateStores,$variableName,$rank,
			"self%".$variableName,"self%".$variableName);
		}
	    } else {
		# Statically-sized variable.
		foreach ( @{$declaration->{'variables'}} ) {
		    my $variableName = &stripVariableName($_);
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    my $store = 1;
		    if ( defined($class) && exists($class->{'stateStorable'}) && exists($class->{'stateStorable'}->{'restoreTo'}) ) {
			foreach ( &List::ExtraUtils::as_array($class->{'stateStorable'}->{'restoreTo'}) ) {
			    my @variables = split(/\s*,\s*/,$_->{'variables'});
			    if ( grep {lc($_) eq lc($variableName)} @variables ) {
				$store = 0;
				$stateStore->{'inputCode'} .= " self%".$variableName."=".$_->{'state'}."\n";
			    }
			}
		    }
		    push(@{$stateStore->{'staticVariables'}},$variableName)
			if ( $store );
		}
	    }
	}
	# Check for a custom state store/restore.
	$stateStore->{'hasCustomStateStore'  } = 1
	    if
	    (
	     $declaration->{'intrinsic'} eq "procedure"
	     &&
	     $declaration->{'variables'}->[0] =~ m/^stateStore=>/
	    );
	$stateStore->{'hasCustomStateRestore'} = 1
	    if
	    (
	     $declaration->{'intrinsic'} eq "procedure"
	     &&
	     $declaration->{'variables'}->[0] =~ m/^stateRestore=>/
	    );
    }
}

1;
