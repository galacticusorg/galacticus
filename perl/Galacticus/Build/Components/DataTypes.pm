# Contains a Perl module which handles data types for the component build system.

package Galacticus::Build::Components::DataTypes;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use LaTeX::Encode;
use Galacticus::Build::Components::Utils qw(%intrinsicTypes);

sub dataObjectPrimitiveName {
    # Construct and return the name and attributes of the primitive data class to use for data of given type and rank.
    my $dataObject = shift();
    my %options;
    (%options) = @_
	if ( $#_ >= 1 );
    # Validate input.
    die "DataTypes::dataObjectPrimitveName: no 'type' specifier present"
	unless ( exists($dataObject->{'type'}) );
    die "DataTypes::dataObjectPrimitveName: no 'rank' or 'shape' specifier present"
	unless ( exists($dataObject->{'rank'}) || exists($dataObject->{'shape'}) );
    die "DataTypes::dataObjectPrimitveName: can not have both 'rank' and 'shape' specifiers present"
	if     ( exists($dataObject->{'rank'}) && exists($dataObject->{'shape'}) );
    # Construct name, type, and attributes.
    my $name = 
	exists($intrinsicTypes{$dataObject->{'type'}}) 
	?
	       $intrinsicTypes{$dataObject->{'type'}} 
        :
	"type(".                      $dataObject->{'type'}.")";    
    my $type = join("",map {ucfirst($_)} split(" ",$dataObject->{'type'}));
    my @attributes;    
    if ( exists($dataObject->{'rank'}) && $dataObject->{'rank'} > 0 ) {
	push(@attributes,"dimension(".join(",",(":") x $dataObject->{'rank'}).")");
	push(@attributes,"allocatable" )
	    unless ( exists($options{'matchOnly'}) && $options{'matchOnly'} );
    }
    if ( exists($dataObject->{'shape'}) ) {
	push(@attributes,"dimension(".$dataObject->{'shape'}.")");
    }
    my $attributeList = scalar(@attributes) > 0 ? ", ".join(", ",@attributes) : "";
    return ($name,$type,$attributeList);
}

sub dataObjectDocName {
    # Construct and return the name of the object to use in documentation for data of given type and rank.
    my $dataObject = shift();
    # Validate input.
    foreach ( "type" ) {
	die "DataTypes::dataObjectPrimitveName: no '".$_."' specifier present"
	    unless ( exists($dataObject->{$_}) );
    }
    # Determine the data object's rank.
    my $rank = exists($dataObject->{'rank'}) ? $dataObject->{'rank'} : 0;
    # Construct the documentation.
    return
	"\\textcolor{red}{\\textless ".
	(
	 exists              ($intrinsicTypes{$dataObject->{'type'}})
	 ?
	         latex_encode($intrinsicTypes{$dataObject->{'type'}})
	 :
	 "type(".latex_encode(                       $dataObject->{'type'} ).")"
	).
	(
	 $rank > 0
	 ?
	 "[".join(",",(":") x $rank)."]"
	 :
	 ""
	).
	"\\textgreater}";
}

sub dataObjectName {
    # Construct and return the name of the object to use for data of given type and rank.
    my $dataObject = shift;
    # Create the object name.
    my $name = "nodeData";
    if ( exists($intrinsicTypes{$dataObject->{'type'}}) ) {
	$name .= join("",map {ucfirst($_)} split(" ",$intrinsicTypes{$dataObject->{'type'}}));
    } else {
	$name .= ucfirst($dataObject->{'type'});
    }
    if ( exists($dataObject->{'rank'}) ) {
	$name .= "Scalar";
    } elsif ( $dataObject->{'type'} ne "void" ) {	
	$name .= $dataObject->{'rank'}."D";
    }
    $name .= "Evolvable"
	if ( $dataObject->{'isEvolvable'} );
    $name =~ s/\s//g;
    return $name;
}

sub dataObjectDefinition {
    # Construct and return the name and attributes of the primitive data class to use for data of given type and rank.
    my $dataObject = shift();
    my %options;
    (%options) = @_
	if ( $#_ >= 1 );
    # Variables to store the object name and attributes.
    my $intrinsicName;
    my $type         ;
    my $label        ;
    my @attributes   ;
    # Validate.
    die "dataObjectDefinition: no 'type' specifier present"
	unless ( exists($dataObject->{'type'}) );
    # Construct properties.
    if ( exists($intrinsicTypes{$dataObject->{'type'}}) ) {
	$intrinsicName = $intrinsicTypes{$dataObject->{'type'}};
    } else {
	$intrinsicName =                               "type"  ;
	$type          =                 $dataObject->{'type'} ;
    }
    $label =ucfirst($dataObject->{'type'});
    if ( exists($dataObject->{'rank'}) ) {
	if ( $dataObject->{'rank'} > 0 ) {
	    push(@attributes,"dimension(".join(",",(":") x $dataObject->{'rank'}).")");
	    push(@attributes,"allocatable" )
		unless ( exists($options{'matchOnly'}) && $options{'matchOnly'} == 1 );
	}
    } else {
	die "dataObjectDefinition: no 'rank' specifier present";
    }
    # Construct the definitions.
    my $dataDefinition;
    $dataDefinition  ->{'intrinsic' }  = $intrinsicName;
    $dataDefinition  ->{'type'      }  = $type
	if ( defined($type) );
    @{$dataDefinition->{'attributes'}} = @attributes
	if ( @attributes    );
    # Return the data definition and label.
    return ($dataDefinition,$label);
}

1;
