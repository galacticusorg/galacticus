# Contains a Perl module which handles creation of null functions.

package Galacticus::Build::Components::NullFunctions;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use NestedMap;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Data::Dumper;
use Galacticus::Build::Components::Utils qw(&isIntrinsic %intrinsicTypes %intrinsicNulls);
use Galacticus::Build::Components::DataTypes;
use Exporter qw(import);
our $VERSION = 1.00;
our @EXPORT_OK = qw(createNullFunction);

# A record of fingerprints of all null functions created.
my %nullFunctionFingerprints;

sub createNullFunction {
    # Return the name of a null function matching the given descriptor. Creates the function if necessary.
    my $build      = shift();
    my $descriptor = shift();
    # Create a unique fingerprint for this null function.
    my $fingerprint = join(":",map {$descriptor->{$_}} ('selfType','attribute','intent')).":".join(":",map {$descriptor->{'property'}->{$_}} ('type','rank'));
    # Create a name for the function.
    my $functionName = 
	"null"                                                                                   .
	join("",map {ucfirst($descriptor              ->{$_})} ('selfType','attribute','intent')).
	join("",map {ucfirst($descriptor->{'property'}->{$_})} ('type'    ,'rank'              ));
    # If a null function with this fingerprint has already been created, we are done. Otherwise, record the fingerprint.
    return $functionName
	if ( exists($nullFunctionFingerprints{$fingerprint}) );
    $nullFunctionFingerprints{$fingerprint} = 1;
    # Construct the data type for the self argument.
    my $selfType = "nodeComponent".($descriptor->{'selfType'} eq "generic" ? "" : $descriptor->{'selfType'});
    # Construct data type for property.
    my $propertyDescriptor;
    if ( &isIntrinsic($descriptor->{'property'}->{'type'}) ) {
	$propertyDescriptor->{'intrinsic'} = $intrinsicTypes{$descriptor->{'property'}->{'type'}};
    } else {
	$propertyDescriptor->{'intrinsic'} = "type";
	$propertyDescriptor->{'type'     } = $descriptor->{'property'}->{'type'};
    }
    # Construct variables needed.
    my @variables;
    my @modules;
    my $returnType;
    if ( $descriptor->{'attribute'} eq "rate" ) {
	$propertyDescriptor->{'variables' } = [ "setValue" ];
	$propertyDescriptor->{'attributes'} = [ "intent(in   )", ($descriptor->{'property'}->{'rank'} > 0 ? "dimension(".join(",",(":") x $descriptor->{'property'}->{'rank'}).")" : ()) ];
	$returnType = "void";
	@variables =
	    (
	     {
		 intrinsic  => "class",
		 type       => $selfType,
		 attributes => [ "intent(".$descriptor->{'intent'}.")" ],
		 variables  => [ "self" ]
	     },
	     $propertyDescriptor,
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(inout)", "optional" ],
		 variables  => [ "interrupt" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "interruptTask", 
		 attributes => [ "intent(inout)", "optional", "pointer" ],
		 variables  => [ "interruptProcedure" ]
	     }
	    );
    } elsif ( $descriptor->{'attribute'} eq "set" || $descriptor->{'attribute'} eq "scale" ) {
	$propertyDescriptor->{'variables' } = [ "setValue" ];
	$propertyDescriptor->{'attributes'} = [ "intent(in   )", ($descriptor->{'property'}->{'rank'} > 0 ? "dimension(".join(",",(":") x $descriptor->{'property'}->{'rank'}).")" : ()) ];
	$returnType = "void";
	@variables =
	    (
	     {
		 intrinsic  => "class",
		 type       => $selfType,
		 attributes => [ "intent(".$descriptor->{'intent'}.")" ],
		 variables  => [ "self" ]
	     },
	     $propertyDescriptor
	    );
	@modules = ( "Error" );
    } elsif ( $descriptor->{'attribute'} eq "get" ) {
	$propertyDescriptor->{'variables' } = [ "setValue" ];
	$propertyDescriptor->{'attributes'} = [ "intent(in   )", ($descriptor->{'property'}->{'rank'} > 0 ? "dimension(".join(",",(":") x $descriptor->{'property'}->{'rank'}).")" : ()) ];
	$returnType = 
	    (exists($propertyDescriptor->{'type'}) ? "type(".$propertyDescriptor->{'type'}.")" : $propertyDescriptor->{'intrinsic'}).
	    ($descriptor->{'property'}->{'rank'} > 0 ? ", dimension(".join(",",(":") x $descriptor->{'property'}->{'rank'})."), allocatable" : "").
	    " => getValue";
	@variables =
	    (
	     {
		 intrinsic  => "class",
		 type       => $selfType,
		 attributes => [ "intent(".$descriptor->{'intent'}.")" ],
		 variables  => [ "self" ]
	     }
	    );
    } elsif ( $descriptor->{'attribute'} eq "analytic" ) {
	$returnType = "void";
	@variables =
	    (
	     {
		 intrinsic  => "class",
		 type       => $selfType,
		 attributes => [ "intent(".$descriptor->{'intent'}.")" ],
		 variables  => [ "self" ]
	     }
	    );
    } else {
	die("createNullFunction: attribute '".$descriptor->{'attribute'}."' not supported");
    }
    # Construct the function.
    my $function =
    {
	type        => $returnType,
	name        => $functionName,
	description => "A null ".$descriptor->{'attribute'}." rate function for a rank ".$descriptor->{'property'}->{'rank'}." {\\normalfont \\ttfamily ".lc($selfType)."} class.\n",
	variables   => \@variables,
	content     => "!\$GLC attributes unused :: ".join(", ",map {@{$_->{'variables'}}} @variables)."\n"
    };
    # Add modules if any required.
    $function->{'modules'} = \@modules
	if ( @modules );
    # Add null return value for get functions.
    if ( $descriptor->{'attribute'} eq "get" ) {
	$function->{'content'} .=
	    "getValue=".
	    (
	     &isIntrinsic($descriptor->{'property'}->{'type'})	    
	     ?
	     $intrinsicNulls{$descriptor->{'property'}->{'type'}}
	     :
	     "null".ucfirst($descriptor->{'property'}->{'type'}).$descriptor->{'property'}->{'rank'}."d"
	    ).
	    "\n";
    }
    # Add error for set functions.
    if ( $descriptor->{'attribute'} eq "set" ) {
	$function->{'content'} .= "call Error_Report('attempt to set value in null component'//{introspection:location})\n";
    }  
    # Insert into the function list (if it is used).
    push(
	@{$build->{'functions'}},
	$function
	);
    return $functionName;
}

1;
