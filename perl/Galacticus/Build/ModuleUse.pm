# Contains a Perl module which implements processing of "moduleUse" directives in the Galacticus build system.

package Galacticus::Build::ModuleUse;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use DateTime;
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
	$buildData->{'moduleUses'}->{$buildData->{'moduleName'}} = "use ".$buildData->{'moduleName'}."\n"
	    if ( $include == 1 );
    } elsif ( $buildData->{'codeType'} eq "c" ) {
	(my $leafName = $buildData->{'currentFileName'}) =~ s/.*?([^\/]+)\.c(pp)??$/$1.o/;
	$buildData->{'moduleUses'}->{$buildData->{'moduleName'}} = "!: ".$ENV{'BUILDPATH'}."/".$leafName."\n";
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
