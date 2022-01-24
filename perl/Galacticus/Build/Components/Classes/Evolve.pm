# Contains a Perl module which handles evolution of component class meta-properties during build.

package Galacticus::Build::Components::Classes::Evolve;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Text::Template 'fill_in_string';
use Galacticus::Build::Components::Utils qw(&offsetName);

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classesEvolve => 
     {
	 classIteratedFunctions =>
	     [
	      \&Build_Rate_Functions    ,
	      \&Build_Scale_Functions   ,
	      \&Build_Inactive_Functions
	     ]
     }
    );

sub Build_Rate_Functions {
    # Build rate setting functions for evolvable meta-properties.
    my $build = shift();
    my $class = shift();
    # Build the function.
    my $classTypeName = "nodeComponent".ucfirst($class->{'name'});
    my $function =
    {
	type        => "void",
	name        => $class->{'name'}."FloatRank0MetaPropertyRate",
	description => "Accumulate to the rate of change of the indexed rank-0 float meta-property of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $classTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "metaPropertyID" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "setValue" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "optional", "intent(inout)" ],
		 variables  => [ "interrupt" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "interruptTask",
		 attributes => [ "optional", "intent(inout)", "pointer" ],
		 variables  => [ "interruptProcedure" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "offset" ]
	     }
	    ]
    };
    # Build the function.
    if ( grep {$class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: interrupt, interruptProcedure, self
CODE
	$code::className          = $class->{'name'};
	$code::offsetNameAll      = &offsetName('all'     ,$class->{'name'},'floatRank0MetaProperties');
	$code::offsetNameActive   = &offsetName('active'  ,$class->{'name'},'floatRank0MetaProperties');
	$code::offsetNameInactive = &offsetName('inactive',$class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.{$className}FloatRank0MetaPropertyCreator(metaPropertyID)) call metaPropertyNoCreator('{$className}',char({$className}FloatRank0MetaPropertyLabels(metaPropertyID)),'float',0)
if (rateComputeState == propertyTypeAll          ) then
 offset={$offsetNameAll}(metaPropertyID)
else if (rateComputeState == propertyTypeActive  ) then
 if (     nodeInactives({$offsetNameAll}(metaPropertyID))) return
 offset={$offsetNameActive}(metaPropertyID)
else if (rateComputeState == propertyTypeInactive) then
 if (.not.nodeInactives({$offsetNameAll}(metaPropertyID))) return
 offset={$offsetNameInactive}(metaPropertyID)
else
 return
end if
nodeRates(offset)=nodeRates(offset)+setValue
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "floatRank0metaPropertyRate"
	}
	);
}

sub Build_Scale_Functions {
    # Build scale setting functions for evolvable meta-properties.
    my $build = shift();
    my $class = shift();
    # Build the function.
    my $classTypeName = "nodeComponent".ucfirst($class->{'name'});
    my $function =
    {
	type        => "void",
	name        => $class->{'name'}."FloatRank0MetaPropertyScale",
	description => "Set the absolute scale of the rank-0 float indexed meta-property of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $classTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "metaPropertyID" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "setValue" ]
	     }
	    ]
    };
    # Build the function.
    if ( grep {$class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$code::className  = $class->{'name'};
	$code::offsetName = &offsetName('all',$class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
if (.not.{$className}FloatRank0MetaPropertyCreator(metaPropertyID)) call metaPropertyNoCreator('{$className}',char({$className}FloatRank0MetaPropertyLabels(metaPropertyID)),'float',0)
nodeScales({$offsetName}(metaPropertyID))=setValue
CODE
	}
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "floatRank0MetaPropertyScale"
	}
	);  
}

sub Build_Inactive_Functions {
    # Build functions to indicate variables which are inactive (i.e. do not appear on the right-hand side of any differential equation being solved) for meta-properties.
    my $build = shift();
    my $class = shift();
    # Build the function.
    my $classTypeName = "nodeComponent".ucfirst($class->{'name'});
    my $function =
    {
	type        => "void",
	name        => $class->{'name'}."FloatRank0MetaPropertyJcbnZr",
	description => "Indicate that the indexed rank-0 float meta-property of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class is inactive for differential equation solving.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $classTypeName,
		 attributes => [ "intent(inout)" ],
		 ,variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "metaPropertyID" ]
	     }
	    ]
    };
    # Build the function.
    if ( grep {$class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$code::offsetName = &offsetName('all',$class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
nodeInactives({$offsetName}(metaPropertyID))=.true.
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "floatRank0MetaPropertyInactive"
	}
	);  
}

1;
