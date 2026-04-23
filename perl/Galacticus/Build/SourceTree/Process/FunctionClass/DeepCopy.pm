# Contains a Perl module which provides deep copy code generation utilities for functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass::DeepCopy;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Process::FunctionClass::Utils qw(stripVariableName declarationRank);
use Exporter 'import';

our @EXPORT_OK = qw(deepCopyDeclarations deepCopyCopiedSelfBlock generateAssignmentAllocatableCode);

# Alias shared state from the parent package.
our $stateStorables;
our $deepCopyActions;
*stateStorables  = \$Galacticus::Build::SourceTree::Process::FunctionClass::stateStorables;
*deepCopyActions = \$Galacticus::Build::SourceTree::Process::FunctionClass::deepCopyActions;

sub deepCopyCopiedSelfBlock {
    # Append the copiedSelf select-type deep copy block to assignments.
    my $deepCopy         = shift();
    my $name             = shift();
    my $declaration      = shift();
    my $nonAbstractClass = shift();
    my $indent           = shift() // "";
    $deepCopy->{'assignments'} .= $indent."if (associated(self%".$name."\%copiedSelf)) then\n";
    $deepCopy->{'assignments'} .= $indent." select type(s => self%".$name."\%copiedSelf)\n";
    $deepCopy->{'assignments'} .= $indent." ".$declaration->{'intrinsic'}." is (".$declaration->{'type'}.")\n";
    $deepCopy->{'assignments'} .= $indent."  destination%".$name." => s\n";
    $deepCopy->{'assignments'} .= $indent." class default\n";
    $deepCopy->{'assignments'} .= $indent."  call Error_Report('copiedSelf has incorrect type'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
    $deepCopy->{'assignments'} .= $indent." end select\n";
    $deepCopy->{'assignments'} .= $indent." call self%".$name."\%copiedSelf\%referenceCountIncrement()\n";
    $deepCopy->{'assignments'} .= $indent."else\n";
    $deepCopy->{'assignments'} .= $indent." allocate(destination%".$name.",mold=self%".$name.")\n";
}

sub generateAssignmentAllocatableCode {
    # Generate code to assign an allocatable (or association-forced) member variable.
    my $assignment     = shift();
    my $declaration    = shift();
    my $name           = shift();
    my $allocated      = shift();  # "allocated" or "associated"
    my $rankMaximumRef = shift();
    $assignment->{'code'} .= "    if (".$allocated."(self%".$name.")) deallocate(self%".$name.")\n";
    $assignment->{'code'} .= "    if (".$allocated."(from%".$name.")) then\n";
    my $rank   = &declarationRank($declaration);
    my $bounds = $rank > 0
	? "(".join(",",map {"lbound(from%".$name.",dim=".$_."):ubound(from%".$name.",dim=".$_.")"} 1..$rank).")"
	: "";
    $assignment->{'code'} .= "      allocate(self%".$name.$bounds.",mold=from%".$name.")\n";
    if ( $rank > 0 && $declaration->{'intrinsic'} eq "type" ) {
	# <workaround type="gfortran" PR="46897" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=46897">
	#   <seeAlso type="gfortran" PR="57696" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=57696"/>
	#   <description>
	#     Type-bound defined assignment not done because multiple part array references would occur in intermediate expressions.
	#   </description>
	# </workaround>
	for(my $i=1;$i<=$rank;++$i) {
	    $assignment->{'code'} .= "      do i".$i."__=lbound(from%".$name.",dim=".$i."),ubound(from%".$name.",dim=".$i.")\n";
	}
	my $indices = join(",",map {"i".$_."__"} 1..$rank);
	$assignment->{'code'} .= "      self%".$name."(".$indices.")=from%".$name."(".$indices.")\n";
	for(my $i=1;$i<=$rank;++$i) {
	    $assignment->{'code'} .= "      end do\n";
	}
	$$rankMaximumRef = $rank
	    if ( $rank > $$rankMaximumRef );
    } else {
	# Simple assignment.
	$assignment->{'code'} .= "      self%".$name."=from%".$name."\n";
    }
    $assignment->{'code'} .= "    end if\n";
}

sub deepCopyDeclarations {
    # Process variable declarations from a node for deep copy.
    my $class              =   shift() ;
    my $nonAbstractClass   =   shift() ;
    my $node               =   shift() ;
    my $declarations       =   shift() ;
    my @ignore             = @{shift()};
    my $lineNumber         =   shift() ;
    my $deepCopy           =   shift() ;
    my $foundDeepCopyNames =   shift() ;
    foreach my $declaration ( &List::ExtraUtils::as_array($declarations) ) {
	# Deep copy of functionClass objects.
	(my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
	    if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
	if
	    (
	     $declaration->{'intrinsic'} eq "class"
	     &&
	     (grep {$_ eq $type    } (keys(%{$stateStorables->{'functionClasses'}}),@{$stateStorables->{'functionClassInstances'}}))
	     &&
	     grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }}
	    )
	{
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		my $name = &stripVariableName($object);
		next
		    if ( grep {lc($_) eq lc($name)} @ignore );
		$deepCopy->{'resetCode'   } .= "if (associated(self%".$name.")) call self%".$name."%deepCopyReset   ()\n";
		$deepCopy->{'finalizeCode'} .= "if (associated(self%".$name.")) call self%".$name."%deepCopyFinalize()\n";
		$deepCopy->{'assignments' } .= "nullify(destination%".$name.")\n";
		$deepCopy->{'assignments' } .= "if (associated(self%".$name.")) then\n";
		# Undo the reference count increment that occurred as a result of the `destination=self` assignment.
		$deepCopy->{'needReferenceCount'} = 1;
		$deepCopy->{'assignments' } .= " referenceCount__=self%".$name."\%referenceCountDecrement()\n";
		&deepCopyCopiedSelfBlock($deepCopy,$name,$declaration,$nonAbstractClass," ");
		$deepCopy->{'assignments' } .= "  call self%".$name."%deepCopy(destination%".$name.")\n";
		$deepCopy->{'assignments' } .= "  self%".$name."%copiedSelf => destination%".$name."\n";
		$deepCopy->{'assignments' } .= "  call destination%".$name."%autoHook()\n";
		$deepCopy->{'assignments' } .= " end if\n";
		$deepCopy->{'assignments' } .= "end if\n";
	    }
	};
	# Deep copy of objects with explicit deep copy actions.
	if
	    (
	     (
	      $declaration->{'intrinsic'} eq "class"
	      ||
	      $declaration->{'intrinsic'} eq "type"
	     )
	     &&
	     (grep {$_->{'type'} eq $type} &List::ExtraUtils::as_array($deepCopyActions->{'deepCopyActions'}))
	    ) {
		my $isAllocatable = grep {$_ eq "allocatable"} @{$declaration->{'attributes'}};
		my $isPointer     = grep {$_ eq "pointer"    } @{$declaration->{'attributes'}};
		my $rank          = &declarationRank($declaration);
		$deepCopy->{'rankMaximum'} = $rank
		    if ( $rank > $deepCopy->{'rankMaximum'} );
		foreach my $variableName ( @{$declaration->{'variableNames'}} ) {
		    $deepCopy->{'assignments'} .= "if (allocated(self%".$variableName.")) then\n"
			if ( $isAllocatable );
		    $deepCopy->{'assignments'} .= "if (associated(self%".$variableName.")) then\n"
			if ( $isPointer     );
		    for(my $i=1;$i<=$rank;++$i) {
			$deepCopy->{'assignments'} .= (" " x $i)."do i".$i."=lbound(self%".$variableName.",dim=".$i."),ubound(self%".$variableName.",dim=".$i.")\n";
		    }
		    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
		    $deepCopy->{'assignments'} .= (" " x $rank)."call destination%".$variableName.$arrayElement."%deepCopyActions()\n";
		    for(my $i=1;$i<=$rank;++$i) {
			$deepCopy->{'assignments'} .= (" " x ($rank+1-$i))."end do\n";
		    }
		    $deepCopy->{'assignments'} .= "end if\n"
			if ( $isAllocatable || $isPointer );
		}
	}
	# Deep copy of HDF5 objects.
	if
	    (
	     $declaration->{'intrinsic'} eq "type"
	     &&
	     $declaration->{'type'     } =~ m/^\s*hdf5object\s*$/i
	    ) {
		$deepCopy->{'modules'}->{'HDF5_Access'} = 1;
		$deepCopy->{'assignments'} .= "!\$ call hdf5Access%set  ()\n";
		$deepCopy->{'assignments'} .= "call self%".$_."%deepCopy(destination%".$_.")\n"
		    foreach ( @{$declaration->{'variables'}} );
		$deepCopy->{'assignments'} .= "!\$ call hdf5Access%unset()\n";
	}
	# Deep copy of non-(class,pointer) functionClass objects.
	if ( exists($class->{'deepCopy'}->{'functionClass'}) ) {
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		my $name = &stripVariableName($object);
		if ( grep {lc($_) eq lc($name)} split(/\s*,\s*/,$class->{'deepCopy'}->{'functionClass'}->{'variables'}) ) {
		    push(@{$foundDeepCopyNames},$name);
		    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
			$deepCopy->{'assignments' } .= "nullify(destination%".$name.")\n";
			$deepCopy->{'assignments' } .= "if (associated(self%".$name.")) then\n";
			$deepCopy->{'resetCode'   } .= "if (associated(self%".$name.")) then\n";
			$deepCopy->{'finalizeCode'} .= "if (associated(self%".$name.")) then\n";
			$deepCopy->{'needReferenceCount'} = 1;
			$deepCopy->{'assignments' } .= "referenceCount__=self%".$name."\%referenceCountDecrement()\n";
			&deepCopyCopiedSelfBlock($deepCopy,$name,$declaration,$nonAbstractClass);
		    }
		    $deepCopy->{'resetCode'   } .= "call self%".$name."%deepCopyReset   ()\n";
		    $deepCopy->{'finalizeCode'} .= "call self%".$name."%deepCopyFinalize()\n";
		    $deepCopy->{'assignments' } .= "call self%".$name."%deepCopy(destination%".$name.")\n";
		    $deepCopy->{'assignments' } .= "self%".$name."%copiedSelf => destination%".$name."\n";
		    $deepCopy->{'assignments' } .= "call destination%".$name."%autoHook()\n";
		    $deepCopy->{'assignments' } .= "end if\n"
			if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} );
		    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
			$deepCopy->{'assignments' } .= "end if\n";
			$deepCopy->{'resetCode'   } .= "end if\n";
			$deepCopy->{'finalizeCode'} .= "end if\n";
		    }
		}
	    }
	}
	# Perform any increments.
	if ( exists($class->{'deepCopy'}->{'increment'}) ) {
	    my @increments = map {{variable => $_}} split(/\s*,\s*/,$class->{'deepCopy'}->{'increment'}->{'variables'});
	    foreach ( @increments ) {
		($_->{'host'} = $_->{'variable'}) =~ s/^([^%]+)%.+/$1/;
	    }
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		my $name = &stripVariableName($object);
		foreach my $increment ( @increments ) {
		    if ( lc($increment->{'host'}) eq lc($name) ) {
			$deepCopy->{'assignments'} .= "!\$omp atomic\n"
			    if ( exists($class->{'deepCopy'}->{'increment'}->{'atomic'}) && $class->{'deepCopy'}->{'increment'}->{'atomic'} eq "yes" );
			$deepCopy->{'assignments'} .= "destination\%".$increment->{'variable'}."=destination\%".$increment->{'variable'}."+1\n";
		    }
		}
	    }
	}
	# Perform any sets.
	if ( exists($class->{'deepCopy'}->{'setTo'}) ) {
	    my @setTos = map {{variable => $_}} split(/\s*,\s*/,$class->{'deepCopy'}->{'setTo'}->{'variables'});
	    foreach ( @setTos ) {
		($_->{'host'} = $_->{'variable'}) =~ s/^([^%]+)%.+/$1/;
	    }
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		my $name = &stripVariableName($object);
		foreach my $setTo ( @setTos ) {
		    if ( lc($setTo->{'host'}) eq lc($name) ) {
			$deepCopy->{'assignments'} .= "destination\%".$setTo->{'variable'}."=".$class->{'deepCopy'}->{'setTo'}->{'value'}."\n";
		    }
		}
	    }
	}
	# Perform any explicit deep copies.
	if ( exists($class->{'deepCopy'}->{'deepCopy'}) ) {
	    my @deepCopies = split(/\s*,\s*/,$class->{'deepCopy'}->{'deepCopy'}->{'variables'});
	    foreach my $object ( @{$declaration->{'variableNames'}} ) {
		foreach my $deepCopy ( @deepCopies ) {
		    if ( lc($object) eq lc($deepCopy) ) {
			$deepCopy->{'assignments'} .= "nullify(destination\%".$object.")\n";
			$deepCopy->{'assignments'} .= "allocate(destination\%".$object.",mold=self\%".$object.")\n";
			if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'copy'}) ) {
			    $deepCopy->{'modules'}->{$class->{'deepCopy'}->{'deepCopy'}->{'module'}} = 1
				if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'module'}) );
			    $deepCopy->{'assignments'} .= "if (associated(self\%".$object.")) call ".$class->{'deepCopy'}->{'deepCopy'}->{'copy'}."(self\%".$object.",destination\%".$object.")\n";
			} else {
			    $deepCopy->{'assignments'} .= "if (associated(self\%".$object.")) call self\%".$object."\%deepCopy(destination\%".$object.")\n";
			}
			if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'reset'}) ) {
			    $deepCopy->{'resetModules'}->{$class->{'deepCopy'}->{'deepCopy'}->{'module'}} = 1
				if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'module'}) );
			    $deepCopy->{'resetCode'} .= "if (associated(self\%".$object.")) call ".$class->{'deepCopy'}->{'deepCopy'}->{'reset'}."(self\%".$object.")\n";
			}
			if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'finalize'}) ) {
			    $deepCopy->{'finalizeModules'}->{$class->{'deepCopy'}->{'deepCopy'}->{'module'}} = 1
				if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'module'}) );
			    $deepCopy->{'finalizeCode'} .= "if (associated(self\%".$object.")) call ".$class->{'deepCopy'}->{'deepCopy'}->{'finalize'}."(self\%".$object.")\n";
			}
		    }
		}
	    }
	}
    }
}

1;
