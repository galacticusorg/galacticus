# Contains a Perl module which implements processing of "functionGlobal" directives in the Galacticus build system.

package Galacticus::Build::FunctionGlobal;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Galacticus::Build::Hooks;
use List::ExtraUtils;
use Fortran::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Hooks::moduleHooks = 
    (
     %Galacticus::Build::Hooks::moduleHooks,
     functionGlobalEstablish => {parse => \&FunctionGlobal_Parse_Directive, generate => \&FunctionGlobal_Establish_Generate_Output},
     functionGlobalPointers  => {parse => \&FunctionGlobal_Parse_Directive, generate => \&FunctionGlobal_Pointers_Generate_Output  }
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
    $buildData->{'functionGlobals'}->{$buildData->{'currentDocument'}->{'unitName'}}->{'module'   } =
	$buildData->{'currentDocument'}->{'module'   }
    if ( exists($buildData->{'currentDocument'}->{'module'   }) );    
}

sub FunctionGlobal_Establish_Generate_Output {
    # Generate output for a "functionGlobalEstablish" directive.
    my $buildData = shift;

    # Assert that we have a file name.
    die("Galacticus::Build::FunctionGlobal::FunctionGlobal_Establish_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::FunctionGlobal\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";

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
    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::FunctionGlobal\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    # Output pointer definitions.
    foreach ( keys(%{$buildData->{'functionGlobals'}}) ) {
	$buildData->{'content'} .= "procedure(".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null), pointer :: ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_ => ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null\n";
    }
    $buildData->{'content'} .= "\ncontains\n\n";
    # Output pointer definitions.
    foreach ( keys(%{$buildData->{'functionGlobals'}}) ) {
	my $opener;
	my $closer;
	if ( $buildData->{'functionGlobals'}->{$_}->{'type'} eq "void" ) {
	    $opener = "subroutine";
	    $closer = "subroutine";
	} else {
	    $opener = "function";
	    $closer = "function";
	}
	my @names;
	foreach ( &List::ExtraUtils::as_array($buildData->{'functionGlobals'}->{$_}->{'arguments'}) ) {
	    if ( $_ =~ m/::(.*)$/ ) {
		push(@names,&Fortran::Utils::Extract_Variables($1));
	    }
	}
	$buildData->{'content'} .= " ".$opener." ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null(".join(",",@names).")\n";
	$buildData->{'content'} .= "  use Error\n";
	if ( exists($buildData->{'functionGlobals'}->{$_}->{'module'}) ) {
	    foreach my $module ( &List::ExtraUtils::as_array($buildData->{'functionGlobals'}->{$_}->{'module'}) ) {
		(my $moduleName = $module) =~ s/\s*,.*//;
		$buildData->{'content'} .= "use".($moduleName eq "ISO_C_Binding" ? ", intrinsic" : "")." :: ".$module."\n";
	    }
	}
	if ( $buildData->{'functionGlobals'}->{$_}->{'type'} ne "void" ) {
	    $buildData->{'content'} .= "  ".$buildData->{'functionGlobals'}->{$_}->{'type'}." :: ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null\n";
	}
	foreach ( &List::ExtraUtils::as_array($buildData->{'functionGlobals'}->{$_}->{'arguments'}) ) {
	    $buildData->{'content'} .= "  ".$_."\n";
	}
	$buildData->{'content'} .= "  !\$GLC attributes unused :: ".join(",",@names)."\n"
	    unless ( scalar(@names) == 0 );
	if ( $buildData->{'functionGlobals'}->{$_}->{'type'} eq "double precision" ) {
	    $buildData->{'content'} .= "  ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null=0.0d0\n";
	} elsif ( $buildData->{'functionGlobals'}->{$_}->{'type'} =~ m/,\s*pointer/ ) {
	    $buildData->{'content'} .= "  ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null => null()\n";
	}
	$buildData->{'content'} .= "  call Error_Report('global functions have not been initialized'//{introspection:location})\n";
	$buildData->{'content'} .= " end ".$closer." ".$buildData->{'functionGlobals'}->{$_}->{'name'}."_Null\n";
    }
}

1;
