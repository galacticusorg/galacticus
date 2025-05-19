# Contains a Perl module which provides various serialization/deserialization functions for tree nodes.

package Galacticus::Build::Components::TreeNodes::Serialization;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodeSerialization =>
     {
	 functions =>
	     [
	      \&Tree_Node_Serialize_ASCII,
	      \&Tree_Node_Serialize_XML  ,
	      \&Tree_Node_Serialize_Raw  ,
	      \&Tree_Node_Deserialize_Raw
	     ]
     }
    );

sub Tree_Node_Serialize_ASCII {
    # Generate a function to produce an ASCII serialization of a tree node object suitable for writing to the terminal.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeASCII",
	description => "Serialize node content to ASCII.",
	modules     =>
	    [
	     "ISO_Varying_String",
	     "Display",
	     "String_Handling"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationVerbosityLevelType",
		 variables  => [ "verbosityLevel" ],
		 attributes => [ "intent(in   )", "optional" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationVerbosityLevelType",
		 variables  => [ "verbosityLevel_" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "message" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=22",
		 variables  => [ "label" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
verbosityLevel_=verbosityLevelStandard
if (present(verbosityLevel)) verbosityLevel_=verbosityLevel
message='Dumping node '
message=message//self%index()
call displayIndent(message,verbosityLevel_)
message='host tree: '
if (associated(self%hostTree)) then
 message=message//self%hostTree%index
else
 message=message//'unhosted'
end if
call displayMessage(message,verbosityLevel_)
call displayIndent('pointers',verbosityLevel_)
CODE
    foreach $code::pointer ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ) {
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
message='{" " x (14-length($pointer))}{$pointer}: '
message=message//self%{$pointer}%index()
call displayMessage(message,verbosityLevel_)
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call displayUnindent('done',verbosityLevel_)
call displayIndent('state',verbosityLevel_)
CODE
    foreach $code::state ( "isPhysicallyPlausible", "isSolvable" ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
message='{" " x (22-length($tate))}{$state}: '
write (label,'(l1)') self%{$state}
message=message//trim(adjustl(label))
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call displayMessage(message,verbosityLevel_)
call displayUnindent('done',verbosityLevel_)
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%serializeASCII(verbosityLevel_)
  end do
end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call displayUnindent('done',verbosityLevel_)
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeASCII"
	}
	);
}

sub Tree_Node_Serialize_XML {
    # Generate a function to produce an XML serialization of a tree node object.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeXML",
	description => "Serialize tree node content as XML.",
	modules     =>
	    [
	     "ISO_Varying_String",
	     "Display",
	     "String_Handling"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=20",
		 variables  => [ "idLabel", "treeLabel" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$omp critical(Node_XML_Dump)
write (  idLabel,'(i20)') self         %index()
write (treeLabel,'(i20)') self%hostTree%index
write (fileHandle,'(a,a,a,a,a)') ' <node tree="',trim(adjustl(treeLabel)),'" id="',trim(adjustl(idLabel)),'" >'
write (fileHandle,'(a)') '  <pointer>'
CODE
    foreach $code::pointer ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (idLabel,'(i20)') self%{$pointer}%index()
write (fileHandle,'(a,a,a)') '   <{$pointer}>',trim(adjustl(idLabel)),'</{$pointer}>'
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,'(a)') '  </pointer>'
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%serializeXML(fileHandle)
  end do
end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,*) ' </node>'
!$omp end critical(Node_XML_Dump)
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeXML"
	}
	);
}

sub Tree_Node_Serialize_Raw {
    # Generate a function to produce a raw (binary) serialization of a tree node object.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeRaw",
	description => "Serialize all content of a tree node to a raw (binary) file.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle) self%isPhysicallyPlausible
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle) allocated(self%component{ucfirst($class->{'name'})})
if (allocated(self%component{ucfirst($class->{'name'})})) then
  select type (component => self%component{ucfirst($class->{'name'})}(1))
  type is (nodeComponent{ucfirst($class->{'name'})})
    write (fileHandle) .false.
  class is (nodeComponent{ucfirst($class->{'name'})})
    write (fileHandle) .true.
  end select
  write (fileHandle) size(self%component{ucfirst($class->{'name'})})
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%serializeRaw(fileHandle)
  end do
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeRaw"
	}
	);
}

sub Tree_Node_Deserialize_Raw {
    # Generate a function to produce an raw (binary) desserialization of a tree node object.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeDeserializeRaw",
	description => "Deserialize a tree node object from a raw (binary) file.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i", "componentCount" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "isAllocated" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileHandle) self%isPhysicallyPlausible
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileHandle) isAllocated
if (isAllocated) then
  read (fileHandle) isAllocated
  read (fileHandle) componentCount
  if (allocated(self%component{ucfirst($class->{'name'})})) deallocate(self%component{ucfirst($class->{'name'})})
  if (isAllocated) then
    allocate(self%component{ucfirst($class->{'name'})}(componentCount),source=default{ucfirst($class->{'name'})}Component)
  else
    allocate(self%component{ucfirst($class->{'name'})}(componentCount),source={ucfirst($class->{'name'})}Class)
  end if
  select type (self)
  type is (treeNode)
    do i=1,componentCount
      self%component{ucfirst($class->{'name'})}(i)%hostNode => self
    end do
  end select
  do i=1,componentCount
    call self%component{ucfirst($class->{'name'})}(i)%deserializeRaw(fileHandle)
  end do
else
   if (allocated(self%component{ucfirst($class->{'name'})})) deallocate(self%component{ucfirst($class->{'name'})})
   allocate(self%component{ucfirst($class->{'name'})}(1))
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "deserializeRaw"
	}
	);
}

1;
