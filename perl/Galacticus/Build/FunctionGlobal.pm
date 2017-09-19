# Contains a Perl module which implements processing of "functionGlobal" directives in the Galacticus build system.

package Galacticus::Build::FunctionGlobal;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use DateTime;
use Data::Dumper;
use Galacticus::Build::Hooks;
use List::ExtraUtils;
use Fortran::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Hooks::moduleHooks = 
    (
     %Galacticus::Build::Hooks::moduleHooks,
     functionGlobalEstablish => {parse => \&FunctionGlobal_Parse_Directive, generate => \&FunctionGlobal_Establish_Generate_Output},
     functionGlobalPointers  => {parse => \&FunctionGlobal_Parse_Directive , generate => \&FunctionGlobal_Pointers_Generate_Output }
    );

sub FunctionGlobal_Parse_Directive {
    # Parse content for a "functionGlobal" directive.
    my $buildData = shift;
    $buildData->{'functionGlobals'}->{$buildData->{'currentDocument'}->{'unitName'}}->{'name'     } = 
	$buildData->{'currentDocument'}->{'unitName' };
    $buildData->{'functionGlobals'}->{$buildData->{'currentDocument'}->{'unitName'}}->{'type'     } = 
	$buildData->{'currentDocument'}->{'type'     };
    $buildData->{'functionGlobals'}->{$buildData->{'currentDocument'}->{'unitName'}}->{'arguments'} =
	$buildData->{'currentDocument'}->{'arguments'}
    if ( exists($buildData->{'currentDocument'}->{'arguments'}) );    
}

sub FunctionGlobal_Establish_Generate_Output {
    # Generate output for a "functionGlobalEstablish" directive.
    my $buildData = shift;

    # Assert that we have a file name.
    die("Galacticus::Build::FunctionGlobal::FunctionGlobal_Establish_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::FunctionGlobal\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n";

    # Output pointer assignments.
    foreach ( keys(%{$buildData->{'functionGlobals'}}) ) {
	$buildData->{'content'}.= $buildData->{'functionGlobals'}->{$_}->{'name'}."_ => ".$buildData->{'functionGlobals'}->{$_}->{'name'}."\n";
    }

}

sub FunctionGlobal_Pointers_Generate_Output {
    # Generate output for a "functionGlobalPointers" directive.
    my $buildData = shift;
    # Assert that we have a file name.
    die("Galacticus::Build::FunctionGlobal::FunctionGlobal_Pointers_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );
    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::FunctionGlobal\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n";
    # Output pointer definitions.
    foreach ( keys(%{$buildData->{'functionGlobals'}}) ) {
	$buildData->{'content'} .= "procedure(".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null), pointer :: ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_ => ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null\n";
    }
    $buildData->{'content'} .= "\ncontains\n\n";
    # Output pointer definitions.
    foreach ( keys(%{$buildData->{'functionGlobals'}}) ) {
	my $opener = "subroutine";
	my $closer = "subroutine";
	if ( $buildData->{'functionGlobals'}->{$_}->{'type'} ne "void" ) {
	    $opener = $buildData->{'functionGlobals'}->{$_}->{'type'}." function";
	    $closer =                                                  "function";
	}
	my @names;
	foreach ( &List::ExtraUtils::as_array($buildData->{'functionGlobals'}->{$_}->{'arguments'}) ) {
	    if ( $_ =~ m/::(.*)$/ ) {
		push(@names,&Fortran::Utils::Extract_Variables($1));
	    }
	}
	$buildData->{'content'} .= " ".$opener." ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null(".join(",",@names).")\n";
	$buildData->{'content'} .= "  use Galacticus_Error\n";
	foreach ( &List::ExtraUtils::as_array($buildData->{'functionGlobals'}->{$_}->{'arguments'}) ) {
	    $buildData->{'content'} .= "  ".$_."\n";
	}
	$buildData->{'content'} .= "  !GCC\$ attributes unused :: ".join(",",@names)."\n";
	if ( $buildData->{'functionGlobals'}->{$_}->{'type'} eq "double precision" ) {
	    $buildData->{'content'} .= "  ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null=0.0d0\n";
	}
	$buildData->{'content'} .= "  call Galacticus_Error_Report('global functions have not been initialized'//{introspection:location})\n";
	$buildData->{'content'} .= " end ".$closer." ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null\n";
    }
}

1;
