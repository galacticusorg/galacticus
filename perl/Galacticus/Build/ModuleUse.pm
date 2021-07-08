# Contains a Perl module which implements processing of "moduleUse" directives in the Galacticus build system.

package Galacticus::Build::ModuleUse;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Galacticus::Build::Hooks;

# Insert hooks for our functions.
%Galacticus::Build::Hooks::moduleHooks = 
    (
     %Galacticus::Build::Hooks::moduleHooks,
     moduleUse => {parse => \&ModuleUse_Parse_Directive, generate => \&ModuleUse_Generate_Output}
    );

sub ModuleUse_Parse_Directive {
    # Parse content for a "moduleUse" directive.
    my $buildData = shift;

    # Assert that we have a module name, file name and code type.
    die("Galacticus::Build::ModuleUse::ModuleUse_Parse_Directive: no moduleName present"      )
	unless ( exists($buildData->{'moduleName'     }) );
    die("Galacticus::Build::ModuleUse::ModuleUse_Parse_Directive: no currentFileName present" )
	unless ( exists($buildData->{'currentFileName'}) );
    die("Galacticus::Build::ModuleUse::ModuleUse_Parse_Directive: no codeType present"        )
	unless ( exists($buildData->{'codeType'       }) );

    # Determine the type of file being processed.
    if ( $buildData->{'codeType'} eq "fortran" ) {
	my $include = 1;
	if ( exists($buildData->{'exclude'}) ) {
	    $include = 0
		if ( $buildData->{'exclude'} eq $buildData->{'moduleName'} );
	}
	my $unitName;
	my $moduleName = $buildData->{'moduleName'};
	if ( exists($buildData->{'currentDocument'}->{'unitName'}) ) {
	    $unitName = $buildData->{'currentDocument'}->{'unitName'};
	    # If a globalized-version is to be called, append the suffix.
	    if ( exists($buildData->{'currentDocument'}->{'useGlobal'}) && $buildData->{'currentDocument'}->{'useGlobal'} eq "yes" ) {
		$unitName   .= "_";
		$moduleName  = "Functions_Global";
	    }
	} else {
	    $unitName = "__all__";
	}
	push(@{$buildData->{'moduleUses'}->{$moduleName}},$unitName)
	    if ( $include == 1 );	
    } elsif ( $buildData->{'codeType'} eq "c" ) {
	(my $leafName = $buildData->{'currentFileName'}) =~ s/.*?([^\/]+)\.c(pp)??$/$1.o/;
	push(@{$buildData->{'cUses'}},$leafName);
    }
}

sub ModuleUse_Generate_Output {
    # Generate output for a "moduleUse" directive.
    my $buildData = shift;

    # Assert that we have a file name.
    die("Galacticus::Build::ModuleUse::ModuleUse_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::ModuleUse\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";

    # Iterate over all modules, and add them to the content.
    foreach my $module ( sort(keys(%{$buildData->{'moduleUses'}})) ) {
	$buildData->{'content'} .= "use ".$module.((grep {$_ eq "__all__"} @{$buildData->{'moduleUses'}->{$module}}) ? "" : ", only : ".join(", ",@{$buildData->{'moduleUses'}->{$module}}))."\n";
    }
    foreach my $leafName ( @{$buildData->{'cUses'}} ) {
	$buildData->{'content'} .= "!: ".$ENV{'BUILDPATH'}."/".$leafName."\n";
    }

}

1;
