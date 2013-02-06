# Contains a Perl module which implements processing of "moduleUse" directives in the Galacticus build system.

package ModuleUse;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use Switch;
use DateTime;
use Data::Dumper;
require Galacticus::Build::Hooks;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
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
    switch ( $buildData->{'codeType'} ) {
	case ( "fortran" ) {
	    my $include = 1;
	    if ( exists($buildData->{'exclude'}) ) {
		$include = 0
		    if ( $buildData->{'exclude'} eq $buildData->{'moduleName'} );
	    }
	    $buildData->{'moduleUses'}->{$buildData->{'moduleName'}} = "use ".$buildData->{'moduleName'}."\n"
		if ( $include == 1 );
	}
	case ( "c" ) {
	    (my $leafName = $buildData->{'currentFileName'}) =~ s/.*?([^\/]+)\.c(pp)??$/$1.o/;
	    $buildData->{'moduleUses'}->{$buildData->{'moduleName'}} = "!: ./work/build/".$leafName."\n";
	}
    }
}

sub ModuleUse_Generate_Output {
    # Generate output for a "moduleUse" directive.
    my $buildData = shift;

    # Assert that we have a file name.
    die("Galacticus::Build::ModuleUse::ModuleUse_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::ModuleUse\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n";

    # Iterate over all modules, and add them to the content.
    foreach my $module ( keys(%{$buildData->{'moduleUses'}}) ) {
	$buildData->{'content'} .= $buildData->{'moduleUses'}->{$module};
    }
}

1;
