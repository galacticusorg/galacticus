# Contains a Perl module which handles serialization of component implementations.

package Galacticus::Build::Components::Implementations::Serialization;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw($fullyQualifiedNameLengthMax $implementationPropertyNameLengthMax &isIntrinsic);
use Galacticus::Build::Components::DataTypes;
use Galacticus::Build::Components::Classes::MetaProperties;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationsSerialization =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_Serialize_ASCII,
	      \&Implementation_Serialize_XML  ,
	      \&Implementation_Serialize_Raw  ,
	      \&Implementation_Deserialize_Raw
	     ]
     }
    );

sub Implementation_Serialize_ASCII {
    # Generate a function to serialize component implementations to ASCII.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."SerializeASCII",
	description => "Serialize the contents of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component to ASCII.",
	modules     =>
	    [
	     "Display",
	     "ISO_Varying_String",
	     "String_Handling"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationVerbosityLevelType",
		 variables  => [ "verbosityLevel" ],
		 attributes => [ "intent(in   )" ]
	     }
	    ],
	content     => ""
    };
    # Add character variables.
    push
	(
	 @{$function->{'variables'}},
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
	);
    # Add counter variables for any rank-1 and meta-properties.
    push
	(
	 @{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "i", "j" ]
	 }
	);
    # Define format labels for different data types.
    %code::formatLabel = 
	(
	 "double"      => "'(e22.16)'",
	 "integer"     => "'(i8)'"    ,
	 "longInteger" => "'(i16)'"   ,
	 "logical"     => "'(l1)'"
	);
    # Skip null components.
    unless ( $code::member->{'name'} eq "null" ) {
	# Serialize the parent type if necessary.
	if ( exists($code::member->{'extends'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%serializeASCII(verbosityLevel)
CODE
	}
	$code::padding = " " x ($fullyQualifiedNameLengthMax-length($code::class->{'name'}));
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call displayIndent('{$class->{'name'}}: {$padding.$member->{'name'}}',verbosityLevel)
CODE
	foreach $code::property ( map {! $_->{'attributes'}->{'isVirtual'} ? $_ : ()} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	    $code::nameLength = length($code::property->{'name'});
	    if ( $code::property->{'data'}->{'rank'} == 0 ) {
		if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (label,{$formatLabel{$property->{'data'}->{'type'}}}) self%{$property->{'name'}}Data
message='{$property->{'name'}}: '//repeat(' ',propertyNameLengthMax-{$nameLength})//label
call displayMessage(message,verbosityLevel)
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
message='{$property->{'name'}}:'
call displayIndent(message,verbosityLevel)
call self%{$property->{'name'}}Data%dump(verbosityLevel)
call displayUnindent('end',verbosityLevel)
CODE
		}
	    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
		if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%{$property->{'name'}}Data)
   write (label,'(i3)') i
   message='{$property->{'name'}}: '//repeat(' ',propertyNameLengthMax-{$nameLength})//trim(label)
   write (label,{$formatLabel{$property->{'data'}->{'type'}}}) self%{$property->{'name'}}Data(i)
   message=message//': '//label
   call displayMessage(message,verbosityLevel)
end do
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%{$property->{'name'}}Data)
   write (label,'(i3)') i
   message='{$property->{'name'}}: '//repeat(' ',propertyNameLengthMax-{$nameLength})//trim(label)
   call displayIndent(message,verbosityLevel)
   call self%{$property->{'name'}}Data(i)%dump(verbosityLevel)
   call displayUnindent('end',verbosityLevel)
end do
CODE
		}
	    }
	}
    }
    # Serialize meta-properties.
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    $code::label    = $metaPropertyType->{'label'};
	    $code::rank     = $metaPropertyType->{'rank' };
	    $code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	    if      ( $metaPropertyType->{'label'} eq "float"   ) {
		$code::format = $code::formatLabel{'double'};
	    } elsif ( $metaPropertyType->{'label'} eq "integer" ) {
		$code::format = $code::formatLabel{'integer'};
	    } elsif ( $metaPropertyType->{'label'} eq "longInteger" ) {
		$code::format = $code::formatLabel{'longInteger'};
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_ASCII: unknown meta-property type");
	    }
	    if      ( $metaPropertyType->{'rank'} == 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 do i=1,size({$class->{'name'}.$prefix}MetaPropertyNames)
  write (label,{$format}) self%{$prefix}MetaProperties(i)
  message=trim({$class->{'name'}.$prefix}MetaPropertyNames(i))//': '//repeat(' ',propertyNameLengthMax-len_trim({$class->{'name'}.$prefix}MetaPropertyNames(i)))//label
  call displayMessage(message,verbosityLevel)
 end do
end if
CODE
	    } elsif ( $metaPropertyType->{'rank'} == 1 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 do i=1,size({$class->{'name'}.$prefix}MetaPropertyNames)
  do j=1,size( self%{$prefix}MetaProperties(i)%values)
   write (label,'(i3)') j
   message=trim({$class->{'name'}.$prefix}MetaPropertyNames(i))//': '//repeat(' ',propertyNameLengthMax-len_trim({$class->{'name'}.$prefix}MetaPropertyNames(i)))//trim(label)
   write (label,{$format}) self%{$prefix}MetaProperties(i)%values(j)
   message=message//': '//label
   call displayMessage(message,verbosityLevel)
  end do
 end do
end if
CODE
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_ASCII: unsupported meta-property rank")
	    }
	}
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call displayUnindent('done',verbosityLevel)
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "serializeASCII", 
	}
	);
}

sub Implementation_Serialize_XML {
    # Generate a function to serialize component implementations to XML.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."SerializeXML",
	description => "Serialize the contents of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component to XML.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     }
	    ],
	content     => ""
    };
    # Add counter variables for any rank-1 and meta-properties.
    push
	(
	 @{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "i", "j" ]
	 }
	);
    # Define format labels for different data types.
    %code::formatLabel = 
	(
	 "double"      => "(e12.6)",
	 "integer"     => "(i8)"   ,
	 "longInteger" => "(i16)"  ,
	 "logical"     => "(l1)"
	);
    # Generate the code.
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(", ",@unused)}
CODE
    }
    # Serialize the parent type if necessary.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%serializeXML(fileHandle)
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,'(a)') '  <{$class->{'name'}} type="{$member->{'name'}}">'
CODE
    foreach $code::property ( map {! $_->{'attributes'}->{'isVirtual'} || $_->{'data'}->{'rank'} == 0 ? $_ : ()} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	if ( ! $code::property->{'attributes'}->{'isVirtual'} ) {
	    if ( $code::property->{'data'}->{'rank'} == 0 ) {
		if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,'(a,{$formatLabel{$property->{'data'}->{'type'}}},a)') '   <{$property->{'name'}}>',self%{$property->{'name'}}Data,'</{$property->{'name'}}>'
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,'(a)') '   <{$property->{'name'}}>'
write (fileHandle,'(a)') '   </{$property->{'name'}}>'
CODE
		}
	    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
		if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%{$property->{'name'}}Data)
   write (fileHandle,'(a,{$formatLabel{$property->{'data'}->{'type'}}},a)') '   <{$property->{'name'}}>',self%{$property->{'name'}}Data(i),'</{$property->{'name'}}>'
end do
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%{$property->{'name'}}Data)
   write (fileHandle,'(a)') '   <{$property->{'name'}}>'
   write (fileHandle,'(a)') '   </{$property->{'name'}}>'
end do
CODE
		}
	    }
	} elsif ( $code::property->{'attributes'}->{'isGettable'} && $code::property->{'rank'} == 0 && &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,'(a,{$formatLabel{$property->{'data'}->{'type'}}},a)') '   <{$property->{'name'}}>',self%{$property->{'name'}}(),'</{$property->{'name'}}>'
CODE
	}
    }
    # Serialize meta-properties.
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    $code::label    = $metaPropertyType->{'label'};
	    $code::rank     = $metaPropertyType->{'rank' };
	    $code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	    if      ( $metaPropertyType->{'label'} eq "float"   ) {
		$code::format = $code::formatLabel{'double'};
	    } elsif ( $metaPropertyType->{'label'} eq "integer" ) {
		$code::format = $code::formatLabel{'integer'};
	    } elsif ( $metaPropertyType->{'label'} eq "longInteger" ) {
		$code::format = $code::formatLabel{'longInteger'};
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_XML: unknown meta-property type");
	    }
	    if      ( $metaPropertyType->{'rank'} == 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 do i=1,size(({$class->{'name'}.$prefix}MetaPropertyNames))
  write (fileHandle,'(a,a,a,{$format},a,a,a)') '   <'//char({$class->{'name'}.$prefix}MetaPropertyNames(i))//'>',self%{$prefix}MetaProperties(i),'</'//char({$class->{'name'}.$prefix}MetaPropertyNames(i))//'>'
 end do
end if
CODE
	    } elsif ( $metaPropertyType->{'rank'} == 1 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 do i=1,size(({$class->{'name'}.$prefix}MetaPropertyNames))
  do j=1,size(self%{$prefix}MetaProperties(i)%values)
   write (fileHandle,'(a,a,a,{$format},a,a,a)') '   <'//char({$class->{'name'}.$prefix}MetaPropertyNames(i))//'>',self%{$prefix}MetaProperties(i)%values(j),'</'//char({$class->{'name'}.$prefix}MetaPropertyNames(i))//'>'
  end do
 end do
end if
CODE
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_XML: unsupported meta-property rank")
	    }
	}
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle,'(a)') '  </{$class->{'name'}}>'
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "serializeXML", 
	}
	);
}

sub Implementation_Serialize_Raw {
    # Generate a function to serialize component implementations to raw (binary) file.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."SerializeRaw",
	description => "Serialize the contents of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component to raw (binary) file.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     }
	    ],
	content     => ""
    };
    # Add a counter variable for any rank-1 properties and meta-properties.
    push
	(
	 @{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    # Generate the code.
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(", ",@unused)}
CODE
    }
    # Serialize the parent type if necessary.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%serializeRaw(fileHandle)
CODE
    }
    foreach $code::property ( grep {! $_->{'attributes'}->{'isVirtual'}} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle) self%{$property->{'name'}}Data
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%{$property->{'name'}}Data%dumpRaw(fileHandle)
CODE
	    }
	} elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle) allocated(self%{$property->{'name'}}Data)
if (allocated(self%{$property->{'name'}}Data)) then
   write (fileHandle) size(self%{$property->{'name'}}Data)
CODE
	    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileHandle) self%{$property->{'name'}}Data
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   do i=1,size(self%{$property->{'name'}}Data)
      call self%{$property->{'name'}}Data(i)%dumpRaw(fileHandle)
   end do
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
	}
    }
    # Serialize meta-properties.
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    $code::label    = $metaPropertyType->{'label'};
	    $code::rank     = $metaPropertyType->{'rank' };
	    $code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	    if      ( $metaPropertyType->{'label'} eq "float"   ) {
		$code::format = $code::formatLabel{'double'};
	    } elsif ( $metaPropertyType->{'label'} eq "integer" ) {
		$code::format = $code::formatLabel{'integer'};
	    } elsif ( $metaPropertyType->{'label'} eq "longInteger" ) {
		$code::format = $code::formatLabel{'longInteger'};
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_Raw: unknown meta-property type");
	    }
	    if      ( $metaPropertyType->{'rank'} == 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) write (fileHandle) self%{$prefix}MetaProperties
CODE
	    } elsif ( $metaPropertyType->{'rank'} == 1 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 do i=1,size({$class->{'name'}.$prefix}MetaPropertyNames)
  write (fileHandle) size(self%{$prefix}MetaProperties(i)%values)
  write (fileHandle) self%{$prefix}MetaProperties(i)%values
 end do
end if
CODE
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_Raw: unsupported meta-property rank")
	    }
	}
    }
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "serializeRaw", 
	}
	);
}

sub Implementation_Deserialize_Raw {
    # Generate a function to deserialize component implementations from raw (binary) file.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."DeserializeRaw",
	description => "Deserialize the contents of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component from raw (binary) file.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     }
	    ],
	content     => ""
    };
    # Add a counter variable for any rank-1 properties and meta-properties.
    push
	(
	 @{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "i", "metaPropertySize" ]
	 }
	);
    # Add variables required for array reads.
    push
	(
	 @{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "arraySize" ]
	 },
	 {
	     intrinsic  => "logical",
	     variables  => [ "isAllocated" ]
	 }
	)
	if ( grep {! $_->{'attributes'}->{'isVirtual'} && $_->{'data'}->{'rank'} == 1} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) );
    # Generate the code.
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(", ",@unused)}
CODE
    }
    # Deserialize the parent type if necessary.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%deserializeRaw(fileHandle)
CODE
    }
    foreach $code::property ( grep {! $_->{'attributes'}->{'isVirtual'}} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileHandle) self%{$property->{'name'}}Data
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%{$property->{'name'}}Data%readRaw(fileHandle)
CODE
	    }
	} elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileHandle) isAllocated
if (isAllocated) then
   read (fileHandle) arraySize
CODE
	    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   allocate(self%{$property->{'name'}}Data(arraySize))
   read (fileHandle) self%{$property->{'name'}}Data
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   allocate(self%{$property->{'name'}}Data(arraySize))
   do i=1,arraySize)
      call self%{$property->{'name'}}Data(i)%readRaw(fileHandle)
   end do
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
	}
    }
    # Deserialize meta-properties.
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    $code::label    = $metaPropertyType->{'label'};
	    $code::rank     = $metaPropertyType->{'rank' };
	    $code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	    if      ( $metaPropertyType->{'label'} eq "float"   ) {
		$code::format = $code::formatLabel{'double'};
	    } elsif ( $metaPropertyType->{'label'} eq "integer" ) {
		$code::format = $code::formatLabel{'integer'};
	    } elsif ( $metaPropertyType->{'label'} eq "longInteger" ) {
		$code::format = $code::formatLabel{'longInteger'};
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_Raw: unknown meta-property type");
	    }
	    if      ( $metaPropertyType->{'rank'} == 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 allocate(self%{$prefix}MetaProperties(size({$class->{'name'}.$prefix}MetaPropertyNames)))
 read (fileHandle) self%{$prefix}MetaProperties
end if
CODE
	    } elsif ( $metaPropertyType->{'rank'} == 1 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames  )) then
 allocate(self%{$prefix}MetaProperties(size({$class->{'name'}.$prefix}MetaPropertyNames)))
 do i=1,size({$class->{'name'}.$prefix}MetaPropertyNames)
  read (fileHandle) metaPropertySize
  allocate(self%{$prefix}MetaProperties(i)%values(metaPropertySize))
  read (fileHandle) self%{$prefix}MetaProperties(i)%values
 end do
end if
CODE
	    } else {
		die("Galacticus::Build::Components::Implementations::Serialization::Implementation_Serialize_Raw: unsupported meta-property rank")
	    }
	}
    }
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "deserializeRaw", 
	}
	);
}

1;
