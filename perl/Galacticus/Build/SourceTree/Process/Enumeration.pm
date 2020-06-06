# Contains a Perl module which implements processing of enumeration directives.

package Galacticus::Build::SourceTree::Process::Enumeration;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use LaTeX::Encode;
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Process::SourceIntrospection;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'enumerations'} = \&Process_Enumerations;

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
	    my $visibility = exists($node->{'directive'}->{'visibility'}) ? $node->{'directive'}->{'visibility'} : "public";
	    my $validator  = exists($node->{'directive'}->{'validator' }) ? $node->{'directive'}->{'validator' } : "no"    ;
	    my $enumerationSource;
	    my $indexing          = exists($node->{'directive'}->{'indexing'}) ? $node->{'directive'}->{'indexing'} : 0;
	    my $i                 = $indexing-1;
	    $enumerationSource .= "  ! Auto-generated enumeration\n";
	    $enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}.ucfirst($_->{'label'})."=".++$i."\n"
		foreach ( &List::ExtraUtils::as_array($node->{'directive'}->{'entry'}) );
	    my $enumerationCount   = $i+1-$indexing;
	    if ( $validator eq "yes" ) {
		$enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}."Min  =0\n";
		$enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}."Max  =".$i."\n";
		$enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}."Count=".$enumerationCount."\n";
	    }
	    $enumerationSource .= "  ! End auto-generated enumeration\n\n";
	    # Create and insert new nodes.
	    my $enumerationTree = &Galacticus::Build::SourceTree::ParseCode($enumerationSource,"Galacticus::Build::SourceTree::Process::Enumeration()");
	    my @enumerationNodes = &Galacticus::Build::SourceTree::Children($enumerationTree);
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,\@enumerationNodes);
	    # Construct validator function as necessary.
	    if ( $validator eq "yes" ) {
		# Generate function code.
		my $functionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."IsValid";
		my $validatorFunction;
		$validatorFunction .= "\n";
		$validatorFunction .= "  ! Auto-generated enumeration function\n";
		$validatorFunction .= "  logical function ".$functionName."(enumerationValue)\n";
		$validatorFunction .= "    !% Validate a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration value.\n";
		$validatorFunction .= "    implicit none\n\n";
		$validatorFunction .= "    integer, intent(in   ) :: enumerationValue\n";
		$validatorFunction .= "    ".$functionName."=(enumerationValue >= 0 .and. enumerationValue < ".$node->{'directive'}->{'name'}."Count)\n";
		$validatorFunction .= "    return\n";
		$validatorFunction .= "  end function ".$functionName."\n";
		$validatorFunction .= "  ! End auto-generated enumeration function\n";
		# Insert into the module.
		my $validatorTree = &Galacticus::Build::SourceTree::ParseCode($validatorFunction,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @validatorNodes = &Galacticus::Build::SourceTree::Children($validatorTree);
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@validatorNodes);
		# Set the visibility.
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$functionName,$visibility);
	    }
	    # Construct encode functions as necessary.
	    if ( exists($node->{'directive'}->{'encodeFunction'}) && $node->{'directive'}->{'encodeFunction'} eq "yes" ) {
		# Generate function code.
		my $encodeFunctionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Encode";
		my $decodeFunctionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Decode";
		my $function;
		$function .= "\n";
		$function .= "  ! Auto-generated enumeration function\n";
		$function .= "  integer function ".$encodeFunctionName."(name,includesPrefix)\n";
		$function .= "    !% Encode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration from a string, returning the appropriate identifier.\n";
		$function .= "    use Galacticus_Error\n";
		$function .= "    implicit none\n\n";
		$function .= "    character(len=*), intent(in   )           :: name\n";
		$function .= "    logical         , intent(in   ), optional :: includesPrefix\n";
		$function .= "    logical                                   :: includesPrefix_\n\n";
		$function .= "    includesPrefix_=.true.\n";
		$function .= "    if (present(includesPrefix)) includesPrefix_=includesPrefix\n";
		for(my $j=0;$j<2;++$j) {
		    if ( $j == 0 ) {
			$function .= "    if (includesPrefix_) then\n";
		    } else {
			$function .= "    else\n";
		    }
		    $function .= "      select case (trim(name))\n";
		    my $i = $indexing-1;
		    foreach ( &List::ExtraUtils::as_array($node->{'directive'}->{'entry'}) ) {
			$function .= "      case ('";
			if ( $j == 0 ) {
			    $function .= $node->{'directive'}->{'name'}.ucfirst($_->{'label'});
			} else {
			    $function .=                                        $_->{'label'} ;
			}
			$function .= "')\n";
			$function .= "        ".$encodeFunctionName."=".++$i."\n";
		    }
		    $function .= "      case default\n";
		    $function .= "      ".$encodeFunctionName."=-1\n";
		    $function .= "      call Galacticus_Error_Report('unrecognized enumeration member ['//trim(name)//']'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		    $function .= "      end select\n";
		}
		$function .= "    end if\n";
		$function .= "    return\n";
		$function .= "  end function ".$encodeFunctionName."\n\n";
		$function .= "  function ".$decodeFunctionName."(enumerationValue,includePrefix)\n";
		$function .= "    !% Decode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration to a string.\n";
		$function .= "    use ISO_Varying_String\n";
		$function .= "    use Galacticus_Error\n";
		$function .= "    implicit none\n\n";
		$function .= "    type   (varying_string)                          :: ".$decodeFunctionName."\n";
		$function .= "    integer                , intent(in   )           :: enumerationValue\n";
		$function .= "    logical                , intent(in   ), optional :: includePrefix\n";
		for(my $j=0;$j<2;++$j) {
		    if ( $j == 0 ) {
			$function .= "    if (present(includePrefix).and.includePrefix) then\n";
			$function .= "      ".$decodeFunctionName."='".$node->{'directive'}->{'name'}."'\n";
		    } else {
			$function .= "    else\n";
			$function .= "      ".$decodeFunctionName."=''\n";
		    }
		    my $i = $indexing-1;
		    $function .= "    select case(enumerationValue)\n";
		    foreach ( &List::ExtraUtils::as_array($node->{'directive'}->{'entry'}) ) {
			$function .= "    case (".++$i.")\n";
			$function .= "       ".$decodeFunctionName."=".$decodeFunctionName."//'";
			if ( $j == 0 ) {
			    $function .= ucfirst($_->{'label'});
			} else {
			    $function .=         $_->{'label'} ;
			}
			$function .= "'\n";
		    }
		    $function .= "    case default\n";
		    $function .= "      call Galacticus_Error_Report('invalid enumeration value'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		    $function .= "    end select\n";
		}
		$function .= "    end if\n";
		$function .= "    return\n";
		$function .= "  end function ".$decodeFunctionName."\n";
		$function .= "  ! End auto-generated enumeration function\n";

		# Insert into the module.
		my $encodeTree = &Galacticus::Build::SourceTree::ParseCode($function,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @encodeNodes = &Galacticus::Build::SourceTree::Children($encodeTree);
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@encodeNodes);
		# Set the visibility.
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$encodeFunctionName,$visibility);
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$decodeFunctionName,$visibility);
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
	    foreach ( &List::ExtraUtils::as_array($node->{'directive'}->{'entry'}) ) {
		print $defHndl "Members:"
		    if ( $first == 1 );
		print $defHndl " & {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}.latex_encode(ucfirst($_->{'label'}))."}\\\\\n";
		$first = 0;
	    }
	    print $defHndl "\\end{tabular}\n";
	    close($defHndl);
	    system("mkdir -p doc/enumerations/specifiers");
	    open(my $specHndl,">doc/enumerations/specifiers/".$node->{'directive'}->{'name'}.".tex");
	    print $specHndl "\\def\\enum".ucfirst($node->{'directive'}->{'name'})."{\\textcolor{red}{\\hyperlink{ht:AutoEnumerations".ucfirst($node->{'directive'}->{'name'})."}{\\textless ".$node->{'directive'}->{'name'}."\\textgreater}}}\n";
	    close($specHndl);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
