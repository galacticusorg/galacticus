# Contains a Perl module which implements processing of enumeration directives.

package Enumerations;
use strict;
use warnings;
use utf8;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Data::Dumper;
use XML::Simple;
use LaTeX::Encode;
require List::ExtraUtils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::processHooks{'enumerations'} = \&Process_Enumerations;

sub Process_Enumerations {
    # Get the tree.
    my $tree      = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "enumeration" && ! $node->{'directive'}->{'processed'} ) {
	    # Assert that our parent is a module or file (for now).
	    die("Process_Enumerations: parent node must be a module or file")
		unless ( $node->{'parent'}->{'type'} eq "module" || $node->{'parent'}->{'type'} eq "file" );
	    # Generate source code for the enumeration.
	    $node->{'directive'}->{'processed'} =  1;
	    my $enumerationSource     ;
	    my $i                 = -1;
	    $enumerationSource .= "  ! Auto-generated enumeration\n";
	    $enumerationSource .= "  integer, parameter, public :: ".$node->{'directive'}->{'name'}.ucfirst($_->{'label'})."=".++$i."\n"
		foreach ( &ExtraUtils::as_array($node->{'directive'}->{'entry'}) );
	    $enumerationSource .= "  ! End auto-generated enumeration\n\n";
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code"            ,
		content    => $enumerationSource,
		firstChild => undef()
	    };
	    # Insert the node.
	    &SourceTree::InsertAfterNode($node,[$newNode]);
	    # Construct encode functions as necessary.
	    if ( exists($node->{'directive'}->{'encodeFunction'}) && $node->{'directive'}->{'encodeFunction'} eq "yes" ) {
		# Generate function code.
		my $functionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Encode";
		my $encodeFunction;
		$encodeFunction .= "\n";
		$encodeFunction .= "  ! Auto-generated enumeration function\n";
		$encodeFunction .= "  integer function ".$functionName."(name,includesPrefix)\n";
		$encodeFunction .= "    !% Decode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration from a string, returning the appropriate identifier.\n";
		$encodeFunction .= "    use Galacticus_Error\n";
		$encodeFunction .= "    implicit none\n\n";
		$encodeFunction .= "    character(len=*), intent(in   )           :: name\n";
		$encodeFunction .= "    logical         , intent(in   ), optional :: includesPrefix\n";
		$encodeFunction .= "    logical                                   :: includesPrefix_\n\n";
		$encodeFunction .= "    includesPrefix_=.true.\n";
		$encodeFunction .= "    if (present(includesPrefix)) includesPrefix_=includesPrefix\n";
		for(my $j=0;$j<2;++$j) {
		    if ( $j == 0 ) {
			$encodeFunction .= "    if (includesPrefix_) then\n";
		    } else {
			$encodeFunction .= "    else\n";
		    }
		    $encodeFunction .= "      select case (trim(name))\n";
		    my $i = -1;
		    foreach ( &ExtraUtils::as_array($node->{'directive'}->{'entry'}) ) {
			$encodeFunction .= "      case ('";
			if ( $j == 0 ) {
			    $encodeFunction .= $node->{'directive'}->{'name'}.ucfirst($_->{'label'});
			} else {
			    $encodeFunction .=                                        $_->{'label'} ;
			}
			$encodeFunction .= "')\n";
			$encodeFunction .= "  	".$functionName."=".++$i."\n";
		    }
		    $encodeFunction .= "      case default\n";
		    $encodeFunction .= "      call Galacticus_Error_Report('".$functionName."','unrecognized mass specifier ['//trim(name)//']')\n";
		    $encodeFunction .= "      end select\n";
		}
		$encodeFunction .= "    end if\n";
		$encodeFunction .= "    return\n";
		$encodeFunction .= "  end function ".$functionName."\n";
		$encodeFunction .= "  ! End auto-generated enumeration function\n";
		# Create a new node.
		my $newNode =
		{
		    type       => "block"         ,
		    content    => $encodeFunction,
		    firstChild => undef(),
		};
		&SourceTree::BuildTree($newNode);
		# Insert into the module.
		&SourceTree::InsertPostContains($node->{'parent'},[$newNode]);
		# Set the visibility.
		&SourceTree::SetVisibility($node->{'parent'},$functionName,"public");
	    }
	    # Create documentation.
	    system("mkdir -p doc/enumerations/definitions");
	    open(my $defHndl,">doc/enumerations/definitions/".$node->{'directive'}->{'name'}.".tex");
	    print $defHndl "\\subsubsection{\\large {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."}}\\hypertarget{ht:AutoEnumerations".ucfirst($node->{'directive'}->{'name'})."}{}\\label{sec:AutoEnumerations".ucfirst($node->{'directive'}->{'name'})."}\\index{enumerations!".$node->{'directive'}->{'name'}."\@{\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."}}\n\n";
	    print $defHndl "\\begin{tabular}{rp{130mm}}\n";
	    print $defHndl "Description: & ".$node->{'directive'}->{'description'}." \\\\\n";
	    my $moduleNode = $node;
	    do {
		$moduleNode = $moduleNode->{'parent'};
	    } until ( ! $moduleNode ||  $moduleNode->{'type'} eq "module" );
	    if ( $moduleNode ) {
		print $defHndl "Provided by: & {\\normalfont \\ttfamily module} \\hyperlink{".$moduleNode->{'parent'}->{'name'}.":".lc($moduleNode->{'name'})."}{\\normalfont \\ttfamily ".latex_encode($moduleNode->{'name'})."} \\\\\n";
	    }
	    my $first = 1;
	    foreach ( &ExtraUtils::as_array($node->{'directive'}->{'entry'}) ) {
		print $defHndl "Members:"
		    if ( $first == 1 );
		print $defHndl " & {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}.ucfirst($_->{'label'})."}\\\\\n";
		$first = 0;
	    }
	    print $defHndl "\\end{tabular}\n";
	    close($defHndl);
	    system("mkdir -p doc/enumerations/specifiers");
	    open(my $specHndl,">doc/enumerations/specifiers/".$node->{'directive'}->{'name'}.".tex");
	    print $specHndl "\\def\\enum".ucfirst($node->{'directive'}->{'name'})."{\\textcolor{red}{\\hyperlink{ht:AutoEnumerations".ucfirst($node->{'directive'}->{'name'})."}{\\textless ".$node->{'directive'}->{'name'}."\\textgreater}}}\n";
	    close($specHndl);
	}
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
