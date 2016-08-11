# Contains a Perl module which provides various ODE solver-related functions for component implementations.

package ODESolver;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use Text::Template 'fill_in_string';
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationODESolver =>
     {
	 functions =>
	     [
	      \&Implementation_ODE_Name_From_Index
	     ]
     }
    );

sub Implementation_ODE_Name_From_Index {
    # Generate a function to return the name of a property given the index of that property in the serialization of a component
    # implementation.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
	    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
	    my $function =
	    {
		type        => "type(varying_string) => name",
		name        => $implementationTypeName."NameFromIndex",
		description => "Return the name of the property of given index for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class.",
		modules     =>
		    [
		     "ISO_Varying_String"
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
			 intrinsic  => "integer",
			 attributes => [ "intent(inout)" ],
			 variables  => [ "count" ]
		     }
		    ]
	    };
	    # Determine if "self" will be used. It is used iff the implementation extends another implementation, or if any
	    # property is evolveable and is not a rank-0 double.
	    undef(@code::unused);
	    push(@code::unused,"self")
		unless (
		    exists($code::member->{'extends'})
		    ||
		    grep
		    {	
			! $_->{'attributes'}->{'isVirtual'  }
			&&	    
			  $_->{'data'      }->{'isEvolvable'}
			&&
			    (
			     $_->{'data'}->{'rank'       } >  0 
			     ||
			     $_->{'data'}->{'type'       } ne "double"
			    )			
		    }
		    &ExtraUtils::hashList($code::member->{'properties'}->{'property'})
		);
	    # Determine if "count" will be used. It is used iff the implementation extends another implementation, or if any
	    # property is evolveable.
	    push(@code::unused,"count")
		unless (
		    exists($code::member->{'extends'})
		    ||
		    grep
		    {		    
			! $_->{'attributes'}->{'isVirtual'  }
			&&	    
			  $_->{'data'      }->{'isEvolvable'}
		    }
		    &ExtraUtils::hashList($code::member->{'properties'}->{'property'})
		);
	    # Build the function.
	    $function->{'content'}  = "";
	    if ( scalar(@code::unused) > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@unused)}
CODE
	    }
	    # If this component is an extension, first call on the extended type.
	    if ( exists($code::member->{'extends'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
name=self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%nameFromIndex(count)
if (count <= 0) return
CODE
	    }
	    # Iterate over properties.
	    foreach $code::property ( &ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
		# Only evolvable, non-virtual properties are included in the ODE solver.
		next
		    unless
		    (
		     ! $code::property->{'attributes'}->{'isVirtual'  }
		     &&
		       $code::property->{'data'      }->{'isEvolvable'}
		    );
		# Find condition for count update. For allocatable properties, condition is that the object be allocated. For
		# non-allocatable properties, always update the count.
		$code::condition = 
		    $code::property->{'data'}->{'rank'} > 0 
		    ? 
		    "if (allocated(self%".$code::property->{'name'}."Data)) "
		    :
		    "";
		# Find the size of the object. Rank-0 double properties always have a count of 1. For other rank-0 types, call
		# their serialization count method. Rank>0 must be double, so simply use the array length.
		$code::count     = 
		    $code::property->{'data'}->{'rank'} > 0 
		    ?
		    "size(self%".$code::property->{'name'}."Data)"
		    :
		    (
		     $code::property->{'data'}->{'type'} eq "double"
		     ?
		     "1" 
		     : 
		     "self%".$code::property->{'name'}."Data%serializeCount()"
		    );
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$condition}count=count-{$count}
CODE
	        $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (count <= 0) then
  name='{$class->{'name'}}:{$member->{'name'}}:{$property->{'name'}}'
  return
end if
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
name='?'
CODE
	    # Add the function to the functions list.
	    push(
		@{$build->{'functions'}},
		$function
		);
	    # Insert a type-binding for this function into the treeNode type.
	    push(
		@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
		{
		    type        => "procedure", 
		    name        => "nameFromIndex", 
		    function    => $implementationTypeName."NameFromIndex",
		    description => "Return the name of the property of given index for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class.",
		    returnType  => "\\textcolor{red}{\\textless varying\\_string\\textgreater}", 
		    arguments   => "\\intzero\\ index\\argin"
		}
		);	    
	}
    }
}

1;
