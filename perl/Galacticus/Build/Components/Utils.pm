# Contains a Perl module which implements hooks for the component build system.

package Utils;
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
use Data::Dumper;
use List::Util qw(max);

# Define a hash into which modules can insert their hooks.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     utils => 
     {
	 gather =>
	     [
	      \&Label_Lengths
	     ]
     }
    );

# Global verbosity level.
our $verbosityLevel = 1;

# Boolean labels.
our @booleanLabel = ( "false", "true" );

# Intrinsic types.
our %intrinsicTypes =
    (
     "integer"     => "integer"                ,
     "longInteger" => "integer(kind=kind_int8)",
     "logical"     => "logical"                ,
     "double"      => "double precision"       ,
     "void"        => "void"
    );

# Maximum lengths of labels (used for formatting).
our $classNameLengthMax;
our $implementationNameLengthMax;
our $fullyQualifiedNameLengthMax;

sub Label_Lengths {
    # Determine the lengths of various types of label for use in formatting.
    my $build = shift();
    # Find maximum lengths.
    $classNameLengthMax          = max map {                      length($_->{'class'})} &ExtraUtils::hashList($build->{'components'});
    $implementationNameLengthMax = max map {length($_->{'name' })                      } &ExtraUtils::hashList($build->{'components'});
    $fullyQualifiedNameLengthMax = max map {length($_->{'name' })+length($_->{'class'})} &ExtraUtils::hashList($build->{'components'});
    # Report.
    if ( $verbosityLevel >= 1 ) {
	print "         --> Maximum label lengths:\n";
	print "            -->           Class: ".$classNameLengthMax         ."\n";
	print "            -->  Implementation: ".$implementationNameLengthMax."\n";
	print "            --> Fully-qualified: ".$fullyQualifiedNameLengthMax."\n";
    }
}

sub isIntrinsic {
    # Return true if the given type matches an intrinsic type.
    my $type = shift();
    return (grep {$_ eq $type} keys(%intrinsicTypes)) == 1 ? 1 : 0;
}

sub isOutputIntrinsic {
    # Return true if the given type matches an outputtable intrinsic type.
    my $type = shift();
    return (grep {$_ eq $type} ("double","integer","longInteger")) == 1 ? 1 : 0;
}

sub padClass {
    # Pad a class name to give nicely aligned formatting in the output code.
    return &padGeneric($classNameLengthMax,@_);
}

sub padImplementation {
    # Pad an implementation name to give nicely aligned formatting in the output code.
    return &padGeneric($implementationNameLengthMax,@_);
}

sub padFullyQualified {
    # Pad a fully-qualified name to give nicely aligned formatting in the output code.
    return &padGeneric($fullyQualifiedNameLengthMax,@_);
}

sub padGeneric {
    # Pad a string to give nicely aligned formatting in the output code.
    my $length     = shift;
    my $text       = shift;
    my @extraPad   = @{$_[0]};
    my $padLength  = $length+$extraPad[0];
    $padLength     = $extraPad[1]
	if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

1;
