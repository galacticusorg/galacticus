# Contains a Perl module which handles the object hierarchy.

package Galacticus::Build::Components::Hierarchy;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     hierarchy => 
     {
	 functions =>
	     [
	      \&Hierarchy_Initialization,
	      \&Hierarchy_Finalization
	     ]
     }
    );

sub Hierarchy_Initialization {
    # Generate a function to initialize the node class hierarchy.
    my $build = shift();
    # Generate the function.
    my %outputConditions;
    my $function =
    {
	type        => "void",
	name        => "nodeClassHierarchyInitialize",
	description => "Initialize the \\glc\\\ node/component class hierarchy.",
	modules     => 
	    [ 
	      "Input_Parameters"  ,
	      "ISO_Varying_String",
	      "Error"
	    ],
	variables   =>
	[
	     {
		 intrinsic  => "type",
		 type       => "inputParameters",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "parameters_" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "componentSelection", "message" ]
	     },
	     map
	     {
		 {
		     intrinsic  => "type",
		     type       => "nodeComponent".$_,
		     variables  => [ "default".$_."Component" ]
		 }
	     }
	     @{$build->{'componentIdList'}}
	    ],
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
if (hierarchyInitialized == 0) then
  !$omp critical (Node_Class_Hierarchy_Initialize)
  if (hierarchyInitialized == 0) then
    ! Parameters controlling output.
CODE
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$code::defaultImplementation = $code::class->{'defaultImplementation'};
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    ! Insert a function call to get the parameter controlling the choice of implementation for this class.
    if (.not.parameters_%isPresent('component{ucfirst($class->{'name'})}')) call parameters_%addParameter('component{ucfirst($class->{'name'})}','{$defaultImplementation}')
    !![
    <inputParameter>
      <name>component{ucfirst($class->{'name'})}</name>
      <variable>componentSelection</variable>
      <source>parameters_</source>
      <defaultValue>var_str('{$defaultImplementation}')</defaultValue>
      <description>Specifies the implementation to be used for the {$class->{'name'}} component of nodes.</description>
      <type>string</type>
      <cardinality>1</cardinality>
    </inputParameter>
    !!]
CODE
    	foreach $code::component ( @{$code::class->{'members'}} ) {
    	    $code::fullName  = ucfirst($code::class->{'name'}).ucfirst($code::component->{'name'});
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    if (componentSelection == '{$component->{'name'}}') then
       allocate(default{ucfirst($class->{'name'})}Component,source=default{$fullName}Component)
        nodeComponent{$fullName}IsActiveValue=.true.
CODE
	    until ( $code::fullName eq "" ) {
		if ( exists($build->{'components'}->{$code::fullName}->{'extends'}) ) {
		    $code::fullName = ucfirst($build->{'components'}->{$code::fullName}->{'extends'}->{'class'}).ucfirst($build->{'components'}->{$code::fullName}->{'extends'}->{'name'});
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
	nodeComponent{$fullName}IsActiveValue=.true.
CODE
		} else {
		    $code::fullName = "";
		}
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    end if
CODE
	    # Insert code to read and parameters controlling outputs.
	    foreach $code::property ( &List::ExtraUtils::hashList($code::component->{'properties'}->{'property'}) ) {
		# Check for output and output condition.
		if (
		    exists($code::property->{'output'}               )                       &&
		    exists($code::property->{'output'}->{'condition'})                       &&
		           $code::property->{'output'}->{'condition'} =~ m/\[\[([^\]]+)\]\]/
		    )
		{		    
		    $code::parameterName = $1;
		    unless ( exists($outputConditions{$code::parameterName}) ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
     !![
     <inputParameter>
      <name>{$parameterName}</name>
      <variable>{$parameterName}</variable>
      <source>parameters_</source>
      <defaultValue>.false.</defaultValue>
      <attachedTo>module</attachedTo>
      <description>Specifies whether the \{\\normalfont \\ttfamily {$property->{'name'}}\} method of the \{\\normalfont \\ttfamily {$component->{'name'}}\} implemention of the \{\\normalfont \\ttfamily {$componentClass}\} component class should be output.</description>
      <type>string</type>
      <cardinality>1</cardinality>
     </inputParameter>
     !!]
CODE
			# Add a module-scope variable to store the output status of this property.
			push(
			    @{$build->{'variables'}},
			    {
				intrinsic => "logical",
				variables => [ $code::parameterName ]
			    }
			    );
			$outputConditions{$code::parameterName} = 1;
		    }
		}
	    }
    	}       
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    if (.not.allocated(default{ucfirst($class->{'name'})}Component)) then
      message='unrecognized class "'//componentSelection//'" for "{$class->{'name'}}" component'
      message=message//char(10)//'  available methods are:'
      {
       join(
            "\n      ",
            map 
             {"message=message//char(10)//'    ".$_->{'name'}."'"} 
             @{$class->{'members'}}
           )
      }
      call Error_Report(message//\{introspection:location\})
    end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    ! Increment the count of hierarchy initializations.
    hierarchyInitialized=hierarchyInitialized+1
  end if
  !$omp end critical (Node_Class_Hierarchy_Initialize)
end if
CODE
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
}

sub Hierarchy_Finalization {
    # Generate a finalization function for the node hierarchy.
    $code::build = shift();
    # Generate the function.
    my $function =
    {
	type        => "void",
	name        => "nodeClassHierarchyFinalize",
	description => "Finalize the \\glc\\\ node/component class hierarchy."
    };
    # Generate the function code.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$omp critical (Node_Class_Hierarchy_Initialize)
hierarchyInitialized=hierarchyInitialized-1
if (hierarchyInitialized == 0) then
 {join(" ",map {"deallocate(default".$_->{'name'}."Component)\n"} &List::ExtraUtils::hashList($build->{'componentClasses'}))}
 call massDistributionsDestroy()
end if
!$omp end critical (Node_Class_Hierarchy_Initialize)
CODE
    # Add the function to the functions list.
    push(
	@{$code::build->{'functions'}},
	$function
	);
}

1;
