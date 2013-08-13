# Contains a Perl module which implements processing of "label" directives in the Galacticus build system.

package Labels;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use DateTime;
use Text::Table;
require Galacticus::Build::Hooks;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     label => {parse => \&Labels_Parse_Directive, generate => \&Labels_Generate_Output}
    );

sub Labels_Parse_Directive {
    # Parse content for a "label" directive.
    my $buildData = shift;

    # Assert that we have a prefix, currentDocument and directive.
    die("Galacticus::Build::Labels::Labels_Parse_Directive: no prefix present"         )
	unless ( exists($buildData->{'prefix'}                    ) );
    die("Galacticus::Build::Labels::Labels_Parse_Directive: no directive present"      )
	unless ( exists($buildData->{'directive'}                 ) );
    die("Galacticus::Build::Labels::Labels_Parse_Directive: no currentDocument present")
	unless ( exists($buildData->{'currentDocument'}           ) );
    die("Galacticus::Build::Labels::Labels_Parse_Directive: no label present"          )
	unless ( exists($buildData->{'currentDocument'}->{'label'}) );

    # Generate a predix for this method.
    ++$buildData->{'labelCount'}->{$buildData->{'prefix'}};

    # Create a table which will be used to build the label code.
    $buildData->{'labels'}->{$buildData->{'directive'}} = 
	Text::Table->new(
	    {
		align => "left"
	    },
	    {
		is_sep => 1,
		body   => " "
	    },
	    {
		align => "left"
	    },
	    {
		is_sep => 1,
		body   => ""
	    },
	    {
		align => "left"
	    },
	    {
		is_sep => 1,
		body   => ""
	    },
	    {
		align => "right"
	    }
	)
	unless ( exists($buildData->{'labels'}->{$buildData->{'directive'}}) );

    # Generate content for this label and insert into the table.
    my @labelContent = 
	(
	 "integer, public, parameter ::",
	 $buildData->{'prefix'    }.  $buildData->{'currentDocument'}->{'label'},
	 "=",
	 $buildData->{'labelCount'}->{$buildData->{'prefix'         }}
	);
    $buildData->{'labels'}->{$buildData->{'directive'}}->add(@labelContent);

}

sub Labels_Generate_Output {
    # Generate output for a "moduleUse" directive.
    my $buildData = shift;

    # Assert that we have a file name and directive present.
    die("Galacticus::Build::ModuleUse::ModuleUse_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );
    die("Galacticus::Build::Labels::Labels_Parse_Directive: no directive present")
	unless ( exists($buildData->{'directive'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::Labels\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n";

    # Iterate over all labels, and add them to the content.
    foreach my $directive ( keys(%{$buildData->{'labels'}}) ) {
	$buildData->{'content'} .= $buildData->{'labels'}->{$directive}->table();
    }

}

1;
