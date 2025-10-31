# Contains a Perl module which handles deferred attributes of properties.

package Galacticus::Build::Components::Properties::Deferred;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Text::Template 'fill_in_string';
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::NullFunctions qw(createNullFunction);
use Galacticus::Build::Components::Properties::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertiesDeferred => 
     {
	 propertyIteratedFunctions =>
	     [
	      \&Properties_Deferred_Pointers      ,
	      \&Properties_Deferred_Get_Functions ,
	      \&Properties_Deferred_Set_Functions ,
	      \&Properties_Deferred_Rate_Functions
	     ]
     }
    );

# Record of previously created pointers.
my %createdPointers;

sub Properties_Deferred_Pointers {
    # Generate procedure pointers for deferred attributes of properties.
    my $build    = shift();
    my $class    = shift();
    my $member   = shift();
    my $property = shift();
    # Skip properties which are not deferred.
    return
	if ( $property->{'attributes' }->{'isDeferred'} eq "" );
    # Determine where to attach.
    my $attachTo = $class->{'name'}.ucfirst($member->{'name'});
    # Iterate over attributes.
    foreach my $attribute ( split(/:/,$property->{'attributes' }->{'isDeferred'}) ) {	
	# Determine function name.
	my $functionLabel = lcfirst($attachTo).ucfirst($property->{'name'}).ucfirst($attribute);
	# Determine if this attribute is deferred and has not yet had a procedure pointer created.
	next
	    unless ( exists($Galacticus::Build::Components::Properties::Utils::attributeAdjective{$attribute}) && $property->{'attributes' }->{$Galacticus::Build::Components::Properties::Utils::attributeAdjective{$attribute}} && ! exists($createdPointers{$functionLabel}) );
	# Construct the template function.
	my $template =
	    $attribute eq "get" 
	    ?
	    $class->{'name'}.ucfirst($member->{'name'}).ucfirst($property->{'name'}).ucfirst($attribute)
	    :
	    &createNullFunction($build,{selfType => $class->{'name'}, attribute => $attribute, property => $property, intent => "inout"});
	# Generate the procedure pointer and a boolean to indicate if is has been attached.
	push(
	    @{$build->{'variables'}},
	    {
		intrinsic  => "procedure",
		type       => $template,
		attributes => [ "pointer" ],
		variables  => [ $functionLabel."Deferred" ]
	    },
	    {
		intrinsic  => "logical",
		variables  => [ $functionLabel."IsAttchdVl=.false." ]
	    },
	    );
	# Record that this procedure pointer has been created.
	$createdPointers{$functionLabel} = 1;
    }
}

sub Properties_Deferred_Get_Functions {
    # Generate wrapper functions to call deferred property get functions.
    my $build       = shift();
    $code::class    = shift();
    $code::member   = shift();
    $code::property = shift();
    # Skip unless the property is gettable, the get function is deferred, and must be built.
    return
	unless (
	    (grep {$_ eq "get"} split(":",$code::property->{'attributes' }->{'isDeferred'}))
	    &&
	                                  $code::property->{'attributes' }->{'isGettable'}           
	    &&	    
	                                  $code::property->{'getFunction'}->{'build'     }
	);
    # Get properties of the data type needed.
    (my $functionTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property);
    my $functionType = 
	                                                                    $functionTypeDescriptor->{'intrinsic' }            .
	(exists($functionTypeDescriptor->{'type'      }) ? "(" .            $functionTypeDescriptor->{'type'      }  .")" : "").
	(exists($functionTypeDescriptor->{'attributes'}) ? ", ".join(", ",@{$functionTypeDescriptor->{'attributes'}})     : "");
    # Build the function.
    my $function =
    {
	type        => $functionType." => propertyValue",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'}).ucfirst($code::property->{'name'})."Get",
	description => "Get the value of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of the {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component using a deferred function.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }	     
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
propertyValue={$class->{'name'}.ucfirst($member->{'name'}).ucfirst($property->{'name'})}GetDeferred(self)
CODE
    # Bind this function to the relevant type.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	{
	    type       => "procedure",
	    name       => $code::property->{'name'}, 
	    descriptor => $function
	}
	);
    # Generate an attacher function.
    &Generate_Deferred_Function_Attacher($code::member,$code::property,$build,"get");
}

sub Properties_Deferred_Set_Functions {
    # Generate wrapper functions to call deferred property set functions.
    my $build       = shift();
    $code::class    = shift();
    $code::member   = shift();
    $code::property = shift();
    # Skip unless the property is settable, the set function is deferred, and must be built.
    return
	unless (
	    (grep {$_ eq "set"} split(":",$code::property->{'attributes' }->{'isDeferred'}))
	    &&
	                                  $code::property->{'attributes' }->{'isSettable'}           
	    &&	    
	                                  $code::property->{'setFunction'}->{'build'     }
	);
    # Get properties of the data type needed.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property,matchOnly => 1);
    push(@{$propertyTypeDescriptor->{'variables' }},"setValue"      );
    push(@{$propertyTypeDescriptor->{'attributes'}},"intent(in   )" );
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'}).ucfirst($code::property->{'name'})."Set",
	description => "Set the value of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of the {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component using a deferred function.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     $propertyTypeDescriptor
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
call {$class->{'name'}.ucfirst($member->{'name'}).ucfirst($property->{'name'})}SetDeferred(self,setValue)
CODE
    # Bind this function to the relevant type.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	{
	    type       => "procedure",
	    name       => $code::property->{'name'}."Set",
	    descriptor => $function
	}
	);
    # Generate an attacher function.
    &Generate_Deferred_Function_Attacher($code::member,$code::property,$build,"set");
}

sub Properties_Deferred_Rate_Functions {
    # Generate wrapper functions to call deferred property rate functions.
    my $build       = shift();
    $code::class    = shift();
    $code::member   = shift();
    $code::property = shift();
    # Skip unless the property is ratetable, and the rate function is deferred.
    return
	unless (
	    (grep {$_ eq "rate"} split(":",$code::property->{'attributes' }->{'isDeferred'}))
	    &&
	                                   $code::property->{'attributes' }->{'isEvolvable'}
	);
    # Get properties of the data type needed.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property,matchOnly => 1);
    push(@{$propertyTypeDescriptor->{'variables' }},"setValue"      );
    push(@{$propertyTypeDescriptor->{'attributes'}},"intent(in   )" );
    # Determine the class to which this method binds.
    my $class       = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    $code::attachTo =                         $code::class->{'name'} .ucfirst($code::member->{'name'});
    # Skip if this function has been created already.
    return
	if ( grep {$_->{'name'} eq $code::property->{'name'}."Rate"} @{$build->{'types'}->{$class}->{'boundFunctions'}} );	    
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'}).ucfirst($code::property->{'name'})."Rate",
	description => "Accumulate the rate of change of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of the {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component using a deferred function.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $class,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     $propertyTypeDescriptor,
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(inout)", "optional" ],
		 variables  => [ "interrupt" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "interruptTask",
		 attributes => [ "pointer", "optional", "intent(inout)" ],
		 variables  => [ "interruptProcedure" ]
	     }
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
call {$attachTo.ucfirst($property->{'name'})}RateDeferred(self,setValue,interrupt,interruptProcedure)
CODE
    # Bind this function to the relevant type.
    push(
	@{$build->{'types'}->{$class}->{'boundFunctions'}},
	{
	    type       => "procedure",
	    name       => $code::property->{'name'}."Rate",
	    descriptor => $function
	}
	);    
    # Generate an attacher function.
    &Generate_Deferred_Function_Attacher($code::member,$code::property,$build,"rate");
}

sub Generate_Deferred_Function_Attacher {
    # Generate functions to attach a function to a deferred property and to query the attachment state.
    my $component    = shift();
    my $property     = shift();
    my $build        = shift();
    my $method       = shift();
    # Determine a method-specific suffix for the function.
    my $methodSuffix = $method eq "get" ? "" : ucfirst($method);
    # Get the component fully-qualified, class and method names.
    my $componentClassName = $component->{'class'             };
    my $componentName      = $component->{'fullyQualifiedName'};
    my $propertyName       = $property ->{'name'              };
    # Construct the function name.
    $code::functionLabel = lcfirst($componentName).ucfirst($propertyName).ucfirst($method);
    my $functionName  = $code::functionLabel."Function";    
    # Skip if this function was already created.
    return
	if ( grep {$_->{'name'} eq $propertyName.$methodSuffix."Function"} @{$build->{'types'}->{"nodeComponent".ucfirst($componentName)}->{'boundFunctions'}} );
    # Construct the attacher function.
    my $attachFunction =
    {
	type        => "void",
	name        => $functionName,
	description => "Set the function to be used for the {\\normalfont \\ttfamily ".$method."} method of the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$componentName."} component.",
	variables   =>
	    [
	     {
		 intrinsic  => "procedure",
		 type       =>
		     $method eq "get"
		     ? $componentName.ucfirst($propertyName).ucfirst($method) 
		     : 
		     &createNullFunction($build,{selfType => $component->{'class'}, attribute => $method, property => $property, intent => "inout"}),
		 variables  => [ "deferredFunction" ],
		 isArgument => 1
	     } 
	    ]
    };
    $attachFunction->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
{$functionLabel}Deferred        => deferredFunction
{$functionLabel}IsAttchdVl =  .true.
CODE
    # Construct the attachment status function.
    my $attachStatusFunction =
    {
	type        => "logical",
	name        => $code::functionLabel."IsAttached",
	description => "Return true if the deferred function used to ".$method." the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$componentName."} component class has been attached.",
    };
    $attachStatusFunction->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
{$functionLabel}IsAttached={$functionLabel}IsAttchdVl
CODE
    # Bind this function to the relevant type.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($componentName)}->{'boundFunctions'}},
	{
	    type       => "procedure"                             , 
	    pass       => "nopass"                                , 
	    name       => $propertyName.$methodSuffix."Function"  , 
	    descriptor => $attachFunction
	},
	{
	    type       => "procedure"                             , 
	    pass       => "nopass"                                , 
	    name       => $propertyName.$methodSuffix."IsAttached", 
	    descriptor => $attachStatusFunction
	}
    );
}

1;
