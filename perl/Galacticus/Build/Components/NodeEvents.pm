# Contains a Perl module which builds the nodeEvent class.

package Galacticus::Build::Components::NodeEvents;
use strict;
use warnings;
no warnings 'once';
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::Components::Utils;
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     nodeEvents => 
     {
	 types     =>
	     [
	      \&Build_Node_Event_Class
	     ],
         interfaces =>
	     [
	      \&Node_Event_Task_Interface	      ,
	      \&Node_Event_Merge_Time_Set_Interface
	     ],
	 functions =>
	     [
	      \&Node_Event_Serialize_Raw              ,
	      \&Node_Event_Deserialize_Raw            ,
	      \&Node_Event_Deserialize_Raw_Polymorphic
	     ]
     }
    );

sub Build_Node_Event_Class {
    # Build the "nodeEvent" class.
    my $build = shift();
    # Define the node event classes.
    @{$build->{'nodeEventClasses'}} =
	(
	 {
	     name        => "nodeEvent",
	     description => "Base class for events attached to nodes.",
	     data        =>
		 [
		  {
		      intrinsic  => "integer",
		      type       => "kind=kind_int8",
		      attributes => [ "public" ],
		      variables  => [ "ID" ]
		  },
		  {
		      intrinsic  => "type",
		      type       => "treeNode",
		      attributes => [ "pointer", "public" ],
		      variables  => [ "node => null()" ]
		  },
		  {
		      intrinsic  => "double precision",
		      attributes => [ "public" ],
		      variables  => [ "time" ]
		  },
		  {
		      intrinsic  => "class",
		      type       => "nodeEvent",
		      attributes => [ "public", "pointer" ],
		      variables  => [ "next => null()" ]
		  },
		  {
		      intrinsic  => "procedure",
		      type       => "nodeEventTask",
		      attributes => [ "public", "pointer" ],
		      variables  => [ "task" ]
		  }
		 ]
	 },
	 {
	     name        => "nodeEventBranchJump",
	     description => "Class for branch jump events attached to nodes.",
	     extends     => "nodeEvent",
	     data        => []
	 },
	 {
	     name        => "nodeEventSubhaloPromotion",
	     description => "Class for subhalo promotion events attached to nodes.",
	     extends     => "nodeEvent",
	     data        => []
	 },
	 {
	     name        => "nodeEventBranchJumpInterTree",
	     description => "Class for inter-tree branch jump events attached to nodes.",
	     extends     => "nodeEvent",
	     data        =>
		 [
		  {
		      intrinsic  => "integer",
		      type       => "kind=c_size_t",
		      attributes => [ "public" ],
		      variables  => [ "splitForestUniqueID" ]
		  },
		  {
		      intrinsic  => "integer",
		      type       => "kind=kind_int8",
		      attributes => [ "public" ],
		      variables  => [ "pairedNodeID" ]
		  },
		  {
		      intrinsic  => "logical",
		      attributes => [ "public" ],
		      variables  => [ "isPrimary", "hasSecondary" ]
		  },
		  {
		      intrinsic  => "class",
		      type       => "*",
		      attributes => [ "pointer", "public" ],
		      variables  => [ "creator" ]
		  },
		  {
		      intrinsic  => "procedure",
		      type       => "nodeEventInterTreeMergeTimeSet",
		      attributes => [ "pointer", "nopass", "public" ],
		      variables  => [ "mergeTimeSet" ]
		  }
		 ]
	 },
	 {
	     name        => "nodeEventSubhaloPromotionInterTree",
	     description => "Class for inter-tree subhalo promotion events attached to nodes.",
	     extends     => "nodeEvent",
	     data        =>
		 [
		  {
		      intrinsic  => "integer",
		      type       => "kind=c_size_t",
		      attributes => [ "public" ],
		      variables  => [ "splitForestUniqueID" ]
		  },
		  {
		      intrinsic  => "integer",
		      type       => "kind=kind_int8",
		      attributes => [ "public" ],
		      variables  => [ "pairedNodeID" ]
		  },
		  {
		      intrinsic  => "logical",
		      attributes => [ "public" ],
		      variables  => [ "isPrimary" ]
		  },
		  {
		      intrinsic  => "class",
		      type       => "*",
		      attributes => [ "pointer", "public" ],
		      variables  => [ "creator" ]
		  },
		  {
		      intrinsic  => "procedure",
		      type       => "nodeEventInterTreeMergeTimeSet",
		      attributes => [ "pointer", "nopass", "public" ],
		      variables  => [ "mergeTimeSet" ]
		  }
		 ]
	 }
	);
    # Build the classes.
    foreach my $class ( @{$build->{'nodeEventClasses'}} ) {
	# Create the type.
	$build->{'types'}->{$class->{'name'}} = {
	    name           => $class->{'name'       },
	    comment        => $class->{'description'},
	    isPublic       => "true",
	    dataContent    => $class->{'data'       },
	};
	$build->{'types'}->{$class->{'name'}}->{'extends'} = $class->{'extends'}
	    if ( exists($class->{'extends'}) );
    }
}

sub Node_Event_Task_Interface {
    # Define an interface for "nodeEvent" task functions.
    my $build = shift();
    # Define the interface.
    $build->{'interfaces'}->{'nodeEventTask'} =
    {
	name      => "nodeEventTask"                  ,
	comment   => "Interface for node event tasks.",
	intrinsic => "logical"                        ,
	data      =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeEvent",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "thisEvent" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "pointer", "intent(inout)" ],
		 variables  => [ "thisNode" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "deadlockStatus" ]
	     }
	    ]		
    };
}

sub Node_Event_Merge_Time_Set_Interface {
    # Define an interface for "nodeEvent" merge time set functions.
    my $build = shift();
    # Define the interface.
    $build->{'interfaces'}->{'nodeEventInterTreeMergeTimeSet'} =
    {
	name      => "nodeEventInterTreeMergeTimeSet"                               ,
	comment   => "Interface for node event inter tree merge time set functions.",
	intrinsic => "void"                                                         ,
	data      =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "*",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]	 
	     },
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "nodeSatellite", "nodeHost" ]
	     }
	    ]		
    };
}

sub Node_Event_Deserialize_Raw {
    # Deserialize a nodeEvent object from a raw (binary) file.
    my $build = shift();
    # Iterate over node event classes.
    foreach $code::class ( @{$build->{'nodeEventClasses'}} ) {
	# Build the function.
	my $function =
	{
	    type        => "void",
	    name        => $code::class->{'name'}."DeserializeRaw",
	    description => "Deserialize a {\\normalfont \\ttfamily ".$code::class->{'name'}."} object from raw file.",
	    content     => "",
	    modules     =>
		[
		 "ISO_C_Binding"
		],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => $code::class->{'name'},
		     variables  => [ "self" ],
		     attributes => [ "intent(inout)" ]
		 },
		 {
		     intrinsic  => "integer",
		     variables  => [ "fileUnit" ],
		     attributes => [ "intent(in   )" ]
		 }
		]
	};
	push(
	    @{$function->{'variables'}},
	    {
		intrinsic  => "integer",
		variables  => [ "pointerAssociated" ]
	    },
	    {
		intrinsic  => "type",
		type       => "c_funptr",
		variables  => [ "functionLocation" ]
	    }
	    )
	    if ( $code::class->{'name'} eq "nodeEvent" );
	if ( exists($code::class->{'extends'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
! Read the parent class.
call self%{$class->{'extends'}}%deserializeRaw(fileUnit)
CODE
	}
	foreach $code::data ( @{$code::class->{'data'}} ) {
	    unless ( grep {$_ eq "pointer" } @{$code::data->{'attributes'}} ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileUnit) {join(",",map {"self%".$_} @{$data->{'variables'}})}
CODE
	    }
	}
	# The task function pointer is handled by transfering the memory address to an integer array.
	if ( $code::class->{'name'} eq "nodeEvent" ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileUnit) pointerAssociated
if (pointerAssociated == 1) then
   read (fileUnit) functionLocation
   call c_f_ProcPointer(functionLocation,self%task)
else
   self%task => null()
end if
CODE
	}
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{$code::class->{'name'}}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => "deserializeRaw"
	    }
	    );
    }
}

sub Node_Event_Serialize_Raw {
    # Serialize a nodeEvent object to a raw (binary) file.
    my $build = shift();
    # Iterate over node event classes.
    $code::classCount = -1;
    foreach $code::class ( @{$build->{'nodeEventClasses'}} ) {
	# Build the function.
	my $function =
	{
	    type        => "void",
	    name        => $code::class->{'name'}."SerializeRaw",
	    description => "Serialize a {\\normalfont \\ttfamily ".$code::class->{'name'}."} object to raw file.",
	    modules     =>
		[
		 "ISO_C_Binding"
		],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => $code::class->{'name'},
		     variables  => [ "self" ],
		     attributes => [ "intent(in   )" ]
		 },
		 {
		     intrinsic  => "integer",
		     variables  => [ "fileUnit" ],
		     attributes => [ "intent(in   )" ]
		 },
		 {
		     intrinsic  => "logical",
		     variables  => [ "includeType" ],
		     attributes => [ "intent(in   ), optional" ]
		 }
		]
	};
	push(
	    @{$function->{'variables'}},
	    {
		intrinsic  => "type",
		type       => "c_funptr",
		variables  => [ "functionLocation" ]
	    }
	    )
	    if ( $code::class->{'name'} eq "nodeEvent" );
	++$code::classCount;
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
! Write an integer indicating the type of this event if requested.
if (.not.present(includeType).or.includeType) write (fileUnit) {$classCount}
CODE
	if ( exists($code::class->{'extends'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
! Serialize the parent class.
call self%{$class->{'extends'}}%serializeRaw(fileUnit,.false.)
CODE
	}
	foreach $code::data ( @{$code::class->{'data'}} ) {
	    unless ( grep {$_ eq "pointer" } @{$code::data->{'attributes'}} ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
write (fileUnit) {join(",",map {"self%".$_} @{$data->{'variables'}})}
CODE
	    }
	}
	# The task function pointer is handled by transfering the memory address to an integer array.
	if ( $code::class->{'name'} eq "nodeEvent" ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (associated(self%task)) then
   functionLocation=c_FunLoc(self%task)
   write (fileUnit) 1
   write (fileUnit) functionLocation
else
   write (fileUnit) 0
end if
CODE
	}
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{$code::class->{'name'}}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => "serializeRaw", 
	    }
	    );
    }
}

sub Node_Event_Deserialize_Raw_Polymorphic {
    # Deserialize a polymorphic nodeEvent object from a raw (binary) file.
    my $build = shift();
    # Build the function.
    my $function =
    {
	type        => "class(nodeEvent), pointer => event",
	name        => "nodeEventBuildFromRaw",
	description => "Build a {\\normalfont \\ttfamily nodeEvent} class object from a raw dump file.",
	variables   =>
	    [
	     {
		 intrinsic  => "integer",
		 variables  => [ "fileUnit" ],
		 attributes => [ "intent(in   )" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "classType" ]
	     }
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
read (fileUnit) classType
select case (classType)
CODE
    $code::classCount = -1;
    foreach $code::class ( @{$build->{'nodeEventClasses'}} ) {
	++$code::classCount;
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
case ({$classCount})
  allocate({$class->{'name'}} :: event)
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
case default
   call Galacticus_Error_Report('unknown class type'//\{introspection:location\})
end select
call event%deserializeRaw(fileUnit)
CODE
    # Insert into the function list.
    push(
	@{$build->{'functions'}},
	$function
	);
}

1;
