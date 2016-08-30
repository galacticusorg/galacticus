# Contains a Perl module which provides various output-related functions for tree nodes.

package Galacticus::Build::Components::TreeNodes::Output;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodeOutput =>
     {
	 functions =>
	     [
	      \&Tree_Node_Output_Count,
	      \&Tree_Node_Output_Names,
	      \&Tree_Node_Output
	     ]
     }
    );

sub Tree_Node_Output_Count {
    # Generate a function to return a count of the number of properties to be output from a tree node.
    my $build = shift();   
    my $function =
    {
	type        => "void",
	name        => "treeNodeOutputCount",
	description => "Increment the count of properties to output for a {\\normalfont \\ttfamily treeNode}.",
	content     => "",
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
		 attributes => [ "intent(inout)" ],
		 variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };    
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%outputCount(integerPropertyCount,doublePropertyCount,time,instance=i)
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
	    name        => "outputCount"
	}
	);
}

sub Tree_Node_Output_Names {
    # Generate a function to return names of properties to be output from a tree node.
    my $build = shift();   
    my $function =
    {
	type        => "void",
	name        => "treeNodeOutputNames",
	description => "Establish the names of properties to output for a {\\normalfont \\ttfamily treeNode}.",
	content     => "",
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
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerPropertyNames", "integerPropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "integerPropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "doubleProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "doublePropertyNames", "doublePropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "doublePropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };    
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance=i)
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
	    name        => "outputNames"
	}
	);
}

sub Tree_Node_Output {
    # Generate a function to populate output buffers with data from a tree node.
    my $build = shift();   
    my $function =
    {
	type        => "void",
	name        => "treeNodeOutput",
	description => "Populate output buffers with properties to output for a {\\normalfont \\ttfamily treeNode}.",
	content     => "",
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
	     attributes => [ "intent(inout)" ], 
	     variables  => [ "integerProperty", "integerBufferCount" ]
	 },
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     attributes => [ "intent(inout)", "dimension(:,:)" ],
	     variables  => [ "integerBuffer" ]
	 },
	 {
	     intrinsic  => "integer", 
	     attributes => [ "intent(inout)" ], 
	     variables  => [ "doubleProperty", "doubleBufferCount" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(inout)", "dimension(:,:)" ],
	     variables  => [ "doubleBuffer" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ], 
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	    ]
    };    
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,instance=i)
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
	    name        => "output"
	}
	);
}

1;
