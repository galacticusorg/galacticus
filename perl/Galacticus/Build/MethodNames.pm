# Contains a Perl module which implements processing of "methodNames" directives in the Galacticus build system.

package MethodNames;
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
use Data::Dumper;
require Galacticus::Build::Hooks;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     methodNames => {parse => \&MethodNames_Parse_Directive, generate => \&MethodNames_Generate_Output}
    );

sub MethodNames_Parse_Directive {
    # Parse content for a "methodNames" directive.
    my $buildData = shift;

    # Assert that we have a file name.
    die("Galacticus::Build::MethodNames::MethodNames_Parse_Directive: no currentFileName present" )
	unless ( exists($buildData->{'currentFileName'}) );
    # Store the directive name.
    $buildData->{'methodNamesDirective'} = $buildData->{'directive'};
    # Store the function name.
    $buildData->{'methodNamesFunction'} = $buildData->{'function'};
    # Store whether the input method is scalar or array.
    $buildData->{'methodNamesRank'} = 0;
    $buildData->{'methodNamesRank'} = $buildData->{'methodRank'}
        if ( exists($buildData->{'methodRank'}) );
    # Store the method name.
    $buildData->{'methodNames'}->{$buildData->{'currentDocument'}->{'methodName'}} = 1
        if ( exists($buildData->{'currentDocument'}->{'methodName'}) );
}

sub MethodNames_Generate_Output {
    # Generate output for a "methodNames" directive.
    my $buildData = shift;

    # Assert that we have a file name.
    die("Galacticus::Build::MethodNames::MethodNames_Generate_Output: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );

    # Generate a timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::MethodNames\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";
    $buildData->{'content'} .= "!  Time: ".$now."\n";

    # Begin the function.
    $buildData->{'content'} .= "subroutine ".$buildData->{'methodNamesFunction'}."(inputMethod)\n";
    $buildData->{'content'} .= "   use Galacticus_Display\n";
    $buildData->{'content'} .= "   use String_Handling\n";
    $buildData->{'content'} .= "   use ISO_Varying_String\n";
    $buildData->{'content'} .= "   implicit none\n";
    $buildData->{'content'} .= "   type(varying_string), intent(in   )";
    $buildData->{'content'} .= ", dimension(:)"
	if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= " :: inputMethod\n";
    $buildData->{'content'} .= "   type(varying_string)";
    $buildData->{'content'} .= ", dimension(size(inputMethod))"
	if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= " :: guessedMethod\n";
    $buildData->{'content'} .= "   logical";
    $buildData->{'content'} .= ", dimension(size(inputMethod))"
	if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= " :: methodMatched\n";
    $buildData->{'content'} .= "   type(varying_string), dimension(".keys(%{$buildData->{'methodNames'}}).") :: availableMethods\n";
    $buildData->{'content'} .= "   type(varying_string) :: message\n";
    $buildData->{'content'} .= "   integer :: i,j,distanceMinimum\n";

    # Construct array of method names.
    my $i = 0;
    foreach my $methodName ( keys(%{$buildData->{'methodNames'}}) ) {
	++$i;
	$buildData->{'content'} .= "   availableMethods(".$i.")='".$methodName."'\n";
    }

    # Check if input methods are matched.
    my $index = "";
    $index = "(i)"
	if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= "do i=1,size(inputMethod)\n"
        if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= "   methodMatched".$index."=any(availableMethods == inputMethod".$index.")\n";
    $buildData->{'content'} .= "   if (.not.methodMatched".$index.") then\n";
    $buildData->{'content'} .= "      distanceMinimum=1000000\n";
    $buildData->{'content'} .= "      do j=1,size(availableMethods)\n";
    $buildData->{'content'} .= "         if (String_Levenshtein_Distance(char(inputMethod".$index."),char(availableMethods(j))) < distanceMinimum) then\n";
    $buildData->{'content'} .= "            distanceMinimum=String_Levenshtein_Distance(char(inputMethod".$index."),char(availableMethods(j)))\n";
    $buildData->{'content'} .= "            guessedMethod".$index."=availableMethods(j)\n";
    $buildData->{'content'} .= "         end if\n";
    $buildData->{'content'} .= "      end do\n";
    $buildData->{'content'} .= "   end if\n";
    $buildData->{'content'} .= "end do\n"
        if ( $buildData->{'methodNamesRank'} == 1 );

    # Output help.
    if ( $buildData->{'methodNamesRank'} == 1 ) {
	$buildData->{'content'} .= "   if (.not.all(methodMatched)) then\n";
    } else {
	$buildData->{'content'} .= "   if (.notmethodMatched) then\n";
    }
    $buildData->{'content'} .= "      message='ERROR: parameter \"".$buildData->{'methodNamesDirective'}."\" value does not correspond to an allowed value'\n";
    $buildData->{'content'} .= "      call Galacticus_Display_Message(message)\n";
    $buildData->{'content'} .= "      message='HELP: allowed values for \"".$buildData->{'methodNamesDirective'}."\" are:'\n";
    $buildData->{'content'} .= "      call Galacticus_Display_Message(message)\n";
    $buildData->{'content'} .= "      do i=1,size(availableMethods)\n";
    $buildData->{'content'} .= "         message='   '//availableMethods(i)\n";
    $buildData->{'content'} .= "         call Galacticus_Display_Message(message)\n";
    $buildData->{'content'} .= "      end do\n";
    $buildData->{'content'} .= "      message='HELP: possibly what you meant was:'\n";
    $buildData->{'content'} .= "      call Galacticus_Display_Message(message)\n";
    $buildData->{'content'} .= "do i=1,size(inputMethod)\n"
        if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= "      message='    '//inputMethod".$index."//' => '//guessedMethod".$index."\n";
    $buildData->{'content'} .= "      call Galacticus_Display_Message(message)\n";
    $buildData->{'content'} .= "end do\n"
        if ( $buildData->{'methodNamesRank'} == 1 );
    $buildData->{'content'} .= "   end if\n";
    
    # Close the function.
    $buildData->{'content'} .= "  return\n";
    $buildData->{'content'} .= "end subroutine ".$buildData->{'methodNamesFunction'}."\n";

}

1;
