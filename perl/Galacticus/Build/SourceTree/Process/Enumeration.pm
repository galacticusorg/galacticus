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
use List::Util qw(max);
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
	    $enumerationSource .= "  type, extends(enumerationType) :: enumeration".$node->{'directive'}->{'name'}."Type\n";
	    $enumerationSource .= "  contains\n";
	    $enumerationSource .= "    !![\n";
	    $enumerationSource .= "    <methods>\n";
	    $enumerationSource .= "      <method method=\"operator(==)\" description=\"Test the equality of two members of the enumeration.\"/>\n";
	    $enumerationSource .= "    </methods>\n";
	    $enumerationSource .= "    !!]\n";
	    $enumerationSource .= "    procedure ::                  enumeration".$node->{'directive'}->{'name'}."IsEqual\n";
	    $enumerationSource .= "    generic   :: operator(==) =>  enumeration".$node->{'directive'}->{'name'}."IsEqual\n";
	    $enumerationSource .= "  end type enumeration".$node->{'directive'}->{'name'}."Type\n";
	    $enumerationSource .= "  type(enumeration".$node->{'directive'}->{'name'}."Type), parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}.ucfirst($_->{'label'})."=enumeration".$node->{'directive'}->{'name'}."Type(".++$i.")\n"
		foreach ( &List::ExtraUtils::as_array($node->{'directive'}->{'entry'}) );
	    my $enumerationCount   = $i+1-$indexing;
	    if ( $validator eq "yes" ) {
		$enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}."Min  =".$indexing        ."\n";
		$enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}."Max  =".$i               ."\n";
		$enumerationSource .= "  integer, parameter, ".$visibility." :: ".$node->{'directive'}->{'name'}."Count=".$enumerationCount."\n";
	    }
	    $enumerationSource .= "  ! End auto-generated enumeration\n\n";
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    "Enumerations" =>
		    {
			intrinsic => 0,
			only      => {enumerationType => 1}
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'}                                                    ,$usesNode  );
	    &Galacticus::Build::SourceTree::SetVisibility             ($node->{'parent'},"enumeration".$node->{'directive'}->{'name'}."Type",$visibility);
	    # Create and insert new nodes.
	    my $enumerationTree = &Galacticus::Build::SourceTree::ParseCode($enumerationSource,"Galacticus::Build::SourceTree::Process::Enumeration()");
	    my @enumerationNodes = &Galacticus::Build::SourceTree::Children($enumerationTree);
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,\@enumerationNodes);
	    # Construct equality operator.
	    my $functionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."IsEqual";
	    my $equalityFunction;
	    $equalityFunction .= "\n";
	    $equalityFunction .= "  ! Auto-generated enumeration function\n";
	    $equalityFunction .= "  pure elemental logical function ".$functionName."(enumerationA,enumerationB) result(isEqual)\n";
	    $equalityFunction .= "    !!{\n";
	    $equalityFunction .= "    Validate a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration value.\n";
	    $equalityFunction .= "    !!}\n";
	    $equalityFunction .= "    implicit none\n\n";
	    $equalityFunction .= "    class(enumeration".$node->{'directive'}->{'name'}."Type), intent(in   ) :: enumerationA, enumerationB\n";
	    $equalityFunction .= "    isEqual=enumerationA%ID == enumerationB%ID\n";
	    $equalityFunction .= "    return\n";
	    $equalityFunction .= "  end function ".$functionName."\n";
	    $equalityFunction .= "  ! End auto-generated enumeration function\n";
	    # Insert into the module.
	    my $equalityTree = &Galacticus::Build::SourceTree::ParseCode($equalityFunction,"Galacticus::Build::SourceTree::Process::Enumeration()");
	    my @equalityNodes = &Galacticus::Build::SourceTree::Children($equalityTree);
	    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@equalityNodes);
	    # Construct validator function as necessary.
	    if ( $validator eq "yes" ) {
		# Generate function code.
		my $functionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."IsValid";
		my $validatorFunction;
		$validatorFunction .= "\n";
		$validatorFunction .= "  ! Auto-generated enumeration function\n";
		$validatorFunction .= "  logical function ".$functionName."(enumerationValue)\n";
		$validatorFunction .= "    !!{\n";
		$validatorFunction .= "    Validate a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration value.\n";
		$validatorFunction .= "    !!}\n";
		$validatorFunction .= "    implicit none\n\n";
		$validatorFunction .= "    type(enumeration".$node->{'directive'}->{'name'}."Type), intent(in   ) :: enumerationValue\n";
		$validatorFunction .= "    ".$functionName."=(enumerationValue%ID >= ".$node->{'directive'}->{'name'}."Min .and. enumerationValue%ID <= ".$node->{'directive'}->{'name'}."Max)\n";
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
	    # Construct encode function as necessary.
	    if ( exists($node->{'directive'}->{'encodeFunction'}) && $node->{'directive'}->{'encodeFunction'} eq "yes" ) {
		# Generate function code.
		my $encodeFunctionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Encode";
		my $onError;
		$onError = $node->{'directive'}->{'errorValue'}
		    if ( exists($node->{'directive'}->{'errorValue'}) );
		my $interface;
		$interface .= " interface ".$encodeFunctionName."\n";
		$interface .= "  module procedure ".$encodeFunctionName."Char\n";
		$interface .= "  module procedure ".$encodeFunctionName."VarStr\n";
		$interface .= " end interface ".$encodeFunctionName."\n\n";
		$interface .= " interface ".$encodeFunctionName."ID\n";
		$interface .= "  module procedure ".$encodeFunctionName."IDChar\n";
		$interface .= "  module procedure ".$encodeFunctionName."IDVarStr\n";
		$interface .= " end interface ".$encodeFunctionName."ID\n\n";
		my $function;
		$function .= "\n";
		$function .= "  ! Auto-generated enumeration functions\n";
		$function .= "  integer function ".$encodeFunctionName."IDVarStr(name,includesPrefix)\n";
		$function .= "    !!{\n";
		$function .= "    Encode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration from a string, returning the appropriate identifier ID.\n";
		$function .= "    !!}\n";
		$function .= "    use :: ISO_Varying_String\n";
		$function .= "    implicit none\n\n";
		$function .= "    type   (varying_string), intent(in   )           :: name\n";
		$function .= "    logical                , intent(in   ), optional :: includesPrefix\n";
		$function .= "    ".$encodeFunctionName."IDVarStr=".$encodeFunctionName."ID(char(name),includesPrefix)\n";
		$function .= "    return\n";
		$function .= "  end function ".$encodeFunctionName."IDVarStr\n\n";
		$function .= "  integer function ".$encodeFunctionName."IDChar(name,includesPrefix)\n";
		$function .= "    !!{\n";
		$function .= "    Encode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration from a string, returning the appropriate identifier ID.\n";
		$function .= "    !!}\n";
		$function .= "    use :: ISO_Varying_String\n";
		$function .= "    implicit none\n\n";
		$function .= "    character(len=*), intent(in   )           :: name\n";
		$function .= "    logical         , intent(in   ), optional :: includesPrefix\n";
		$function .= "    type(enumeration".$node->{'directive'}->{'name'}."Type) :: member\n";
		$function .= "    member=".$encodeFunctionName."(name,includesPrefix)\n";
		$function .= "    ".$encodeFunctionName."IDChar=member%ID\n";
		$function .= "    return\n";
		$function .= "  end function ".$encodeFunctionName."IDChar\n\n";
		$function .= "  function ".$encodeFunctionName."VarStr(name,includesPrefix".($onError ? "" : ",status").")\n";
		$function .= "    !!{\n";
		$function .= "    Encode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration from a string, returning the appropriate identifier.\n";
		$function .= "    !!}\n";
		$function .= "    use :: ISO_Varying_String\n";
		$function .= "    implicit none\n\n";
		$function .= "    type   (enumeration".$node->{'directive'}->{'name'}."Type) :: ".$encodeFunctionName."VarStr\n";
		$function .= "    type   (varying_string), intent(in   )           :: name\n";
		$function .= "    logical                , intent(in   ), optional :: includesPrefix\n";
		$function .= "    integer                , intent(  out), optional :: status\n"
		    unless ( $onError );
		$function .= "    ".$encodeFunctionName."VarStr=".$encodeFunctionName."(char(name),includesPrefix".($onError ? "" : ",status").")\n";
		$function .= "    return\n";
		$function .= "  end function ".$encodeFunctionName."VarStr\n\n";
		$function .= "  function ".$encodeFunctionName."Char(name,includesPrefix".($onError ? "" : ",status").")\n";
		$function .= "    !!{\n";
		$function .= "    Encode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration from a string, returning the appropriate identifier.\n";
		$function .= "    !!}\n";
		$function .= "    use :: Error             , only : Error_Report, errorStatusSuccess, errorStatusFail\n"
		    unless ( $onError );
		$function .= "    use :: ISO_Varying_String, only : var_str     , operator(//)\n";
		$function .= "    implicit none\n\n";
		$function .= "    type   (enumeration".$node->{'directive'}->{'name'}."Type) :: ".$encodeFunctionName."Char\n";
		$function .= "    character(len=*), intent(in   )           :: name\n";
		$function .= "    logical         , intent(in   ), optional :: includesPrefix\n";
		$function .= "    integer         , intent(  out), optional :: status\n"
		    unless ( $onError );
		$function .= "    logical                                   :: includesPrefix_\n\n";
		$function .= "    includesPrefix_=.true.\n";
		$function .= "    if (present(includesPrefix)) includesPrefix_=includesPrefix\n";
		$function .= "    if (present(status)) status=errorStatusSuccess\n"
		    unless ( $onError );
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
			$function .= "        ".$encodeFunctionName."Char=enumeration".$node->{'directive'}->{'name'}."Type(".++$i.")\n";
		    }
		    $function .= "      case default\n";
		    if ( $onError ) {
			$function .= "      ".$encodeFunctionName."Char=enumeration".$node->{'directive'}->{'name'}."Type(".$onError.")\n";
		    } else {
			$function .= "      ".$encodeFunctionName."Char=enumeration".$node->{'directive'}->{'name'}."Type(-1)\n";
			$function .= "      if (present(status)) then\n";
			$function .= "         status=errorStatusFail\n";
			$function .= "      else\n";
			$function .= "         call Error_Report(var_str('unrecognized enumeration member [')//trim(name)//']'//enumeration".ucfirst($node->{'directive'}->{'name'})."Describe()//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";	
			$function .= "      end if\n";
		    }
		    $function .= "      end select\n";
		}
		$function .= "    end if\n";
		$function .= "    return\n";
		$function .= "  end function ".$encodeFunctionName."Char\n\n";
		$function .= "  ! End auto-generated enumeration functions\n";
		# Insert into the module.
		my $encodeTree = &Galacticus::Build::SourceTree::ParseCode($function,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @encodeNodes = &Galacticus::Build::SourceTree::Children($encodeTree);
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@encodeNodes);
		my $interfaceTree = &Galacticus::Build::SourceTree::ParseCode($interface,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @interfaceNodes = &Galacticus::Build::SourceTree::Children($interfaceTree);
		&Galacticus::Build::SourceTree::InsertPreContains($node->{'parent'},\@interfaceNodes);
		# Set the visibility.
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$encodeFunctionName,$visibility);
	    }
	    # Construct decode function as necessary.
	    if ( exists($node->{'directive'}->{'decodeFunction'}) && $node->{'directive'}->{'decodeFunction'} eq "yes" ) {
		# Generate function code.
		my $decodeFunctionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Decode";
		my $interface;
		$interface .= " interface ".$decodeFunctionName."\n";
		$interface .= "  module procedure ".$decodeFunctionName."Enumerator\n";
		$interface .= "  module procedure ".$decodeFunctionName."ID\n";
		$interface .= " end interface ".$decodeFunctionName."\n\n";
		my $function;
		$function .= "\n";
		$function .= "  ! Auto-generated enumeration function\n";
		$function .= "  function ".$decodeFunctionName."Enumerator(enumerationValue,includePrefix)\n";
		$function .= "    !!{\n";
		$function .= "    Decode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration to a string.\n";
		$function .= "    !!}\n";
		$function .= "    use ISO_Varying_String\n";
		$function .= "    implicit none\n\n";
		$function .= "    type   (varying_string)                                                             :: ".$decodeFunctionName."Enumerator\n";
		$function .= "    type   (enumeration".$node->{'directive'}->{'name'}."Type), intent(in   )           :: enumerationValue\n";
		$function .= "    logical                                                   , intent(in   ), optional :: includePrefix\n\n";
		$function .= "    ".$decodeFunctionName."Enumerator=".$decodeFunctionName."(enumerationValue%ID,includePrefix)\n";
		$function .= "    return\n";
		$function .= "  end function ".$decodeFunctionName."Enumerator\n";
		$function .= "  function ".$decodeFunctionName."ID(enumerationValue,includePrefix)\n";
		$function .= "    !!{\n";
		$function .= "    Decode a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration to a string.\n";
		$function .= "    !!}\n";
		$function .= "    use ISO_Varying_String\n";
		$function .= "    use Error\n"
		    unless ( exists($node->{'directive'}->{'errorValue'}) );
		$function .= "    implicit none\n\n";
		$function .= "    type   (varying_string)                          :: ".$decodeFunctionName."ID\n";
		$function .= "    integer                , intent(in   )           :: enumerationValue\n";
		$function .= "    logical                , intent(in   ), optional :: includePrefix\n";
		for(my $j=0;$j<2;++$j) {
		    if ( $j == 0 ) {
			$function .= "    if (present(includePrefix).and.includePrefix) then\n";
			$function .= "      ".$decodeFunctionName."ID='".$node->{'directive'}->{'name'}."'\n";
		    } else {
			$function .= "    else\n";
			$function .= "      ".$decodeFunctionName."ID=''\n";
		    }
		    my $i = $indexing-1;
		    $function .= "    select case(enumerationValue)\n";
		    foreach ( &List::ExtraUtils::as_array($node->{'directive'}->{'entry'}) ) {
			$function .= "    case (".++$i.")\n";
			$function .= "       ".$decodeFunctionName."ID=".$decodeFunctionName."ID//'";
			if ( $j == 0 ) {
			    $function .= ucfirst($_->{'label'});
			} else {
			    $function .=         $_->{'label'} ;
			}
			$function .= "'\n";
		    }
		    $function .= "    case default\n";
		    if ( exists($node->{'directive'}->{'errorValue'}) ) {
			$function .= "      ".$decodeFunctionName."ID=".$decodeFunctionName."ID//'Error'\n";
		    } else {
			$function .= "      call Error_Report('invalid enumeration value'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		    }
		    $function .= "    end select\n";
		}
		$function .= "    end if\n";
		$function .= "    return\n";
		$function .= "  end function ".$decodeFunctionName."ID\n";
		$function .= "  ! End auto-generated enumeration function\n";
		# Insert into the module.
		my $decodeTree = &Galacticus::Build::SourceTree::ParseCode($function,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @decodeNodes = &Galacticus::Build::SourceTree::Children($decodeTree);
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@decodeNodes);
		my $interfaceTree = &Galacticus::Build::SourceTree::ParseCode($interface,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @interfaceNodes = &Galacticus::Build::SourceTree::Children($interfaceTree);
		&Galacticus::Build::SourceTree::InsertPreContains($node->{'parent'},\@interfaceNodes);
		# Set the visibility.
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$decodeFunctionName,$visibility);
	    }
	    # Create description function.
	    if ( exists($node->{'directive'}->{'decodeFunction'}) && $node->{'directive'}->{'decodeFunction'} eq "yes" ) {
	    	my $functionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Description";
		my $interface;
		$interface .= " interface ".$functionName."\n";
		$interface .= "  module procedure ".$functionName."Enumerator\n";
		$interface .= "  module procedure ".$functionName."ID\n";
		$interface .= " end interface ".$functionName."\n\n";
		my $descriptorFunctionWrapper;
		$descriptorFunctionWrapper .= "\n";
		$descriptorFunctionWrapper .= "  ! Auto-generated enumeration function\n";
		$descriptorFunctionWrapper .= "  function ".$functionName."Enumerator(enumerationValue)\n";
		$descriptorFunctionWrapper .= "    !!{\n";
		$descriptorFunctionWrapper .= "    Return a description of a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration member.\n";
		$descriptorFunctionWrapper .= "    !!}\n";
		$descriptorFunctionWrapper .= "    use ISO_Varying_String\n";
		$descriptorFunctionWrapper .= "    implicit none\n\n";
		$descriptorFunctionWrapper .= "    type   (varying_string)                                                             :: ".$functionName."Enumerator\n";
		$descriptorFunctionWrapper .= "    type   (enumeration".$node->{'directive'}->{'name'}."Type), intent(in   )           :: enumerationValue\n\n";
		$descriptorFunctionWrapper .= "    ".$functionName."Enumerator=".$functionName."(enumerationValue%ID)\n";
		$descriptorFunctionWrapper .= "    return\n";
		$descriptorFunctionWrapper .= "  end function ".$functionName."Enumerator\n";
		my $descriptorFunction;
		$descriptorFunction .= "  function ".$functionName."ID(enumerationValue) result(description)\n";
		$descriptorFunction .= "    !!{\n";
		$descriptorFunction .= "    Return a description of a {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration value.\n";
		$descriptorFunction .= "    !!}\n";
		$descriptorFunction .= "    use :: ISO_Varying_String, only : varying_string, assignment(=)\n";
		$descriptorFunction .= "    implicit none\n";
		$descriptorFunction .= "    type(varying_string) :: description\n";
		$descriptorFunction .= "    integer, intent(in) :: enumerationValue\n\n";
		my $description      = "    select case (enumerationValue)\n";
		my @entries          = &List::ExtraUtils::as_array($node->{'directive'}->{'entry'});
		my $i = $indexing-1;
		foreach my $entry ( @entries ) {
		    $description .= "   case (".++$i.")\n";
		    $description .= "    description='".(exists($entry->{'description'}) ? $entry->{'description'} : "")."'\n";
		}
		$description        .= "    end select\n";
		$description        .= "    return\n";
		$descriptorFunction .= "  end function ".$functionName."ID\n";
		$descriptorFunction .= "  ! End auto-generated enumeration function\n";
		# Insert into the module.
		my $descriptorWrapperTree  = &Galacticus::Build::SourceTree::ParseCode($descriptorFunctionWrapper,"Galacticus::Build::SourceTree::Process::Enumeration()", instrument => 0);
		my $descriptorTree         = &Galacticus::Build::SourceTree::ParseCode($descriptorFunction       ,"Galacticus::Build::SourceTree::Process::Enumeration()", instrument => 0);
		my @descriptorNodes        = &Galacticus::Build::SourceTree::Children ($descriptorTree                                                                                    );
		my @descriptorWrapperNodes = &Galacticus::Build::SourceTree::Children ($descriptorWrapperTree                                                                             );
		my $newNode = $descriptorNodes[0];
		while ( $newNode->{'type'} ne "function" ) {
		    $newNode = $newNode->{'sibling'};
		}
		$newNode = $newNode->{'firstChild'};
		while ( defined($newNode->{'sibling'}) ) {
		    $newNode = $newNode->{'sibling'};
		}
		my $describeNode =
		{
		    type => "code",
		    content => $description
		};
		&Galacticus::Build::SourceTree::InsertAfterNode($newNode,[$describeNode]);
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@descriptorNodes       );
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@descriptorWrapperNodes);
		my $interfaceTree = &Galacticus::Build::SourceTree::ParseCode($interface,"Galacticus::Build::SourceTree::Process::Enumeration()");
		my @interfaceNodes = &Galacticus::Build::SourceTree::Children($interfaceTree);
		&Galacticus::Build::SourceTree::InsertPreContains($node->{'parent'},\@interfaceNodes);
		# Set the visibility.
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$functionName,$visibility);
	    }
	    # Create describe function.
	    {
	    	my $functionName = "enumeration".ucfirst($node->{'directive'}->{'name'})."Describe";
		my $descriptorFunction;
		$descriptorFunction .= "\n";
		$descriptorFunction .= "  ! Auto-generated enumeration function\n";
		$descriptorFunction .= "  function ".$functionName."() result(description)\n";
		$descriptorFunction .= "    !!{\n";
		$descriptorFunction .= "    Return a description of the {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."} enumeration.\n";
		$descriptorFunction .= "    !!}\n";
		$descriptorFunction .= "    use :: ISO_Varying_String, only : varying_string, var_str, operator(//)\n";
		$descriptorFunction .= "    implicit none\n";
		$descriptorFunction .= "    type(varying_string) :: description\n\n";
		my $description      = "    description=var_str(char(10))//\"Enumeration '".$node->{'directive'}->{'name'}."' has the following members:\"\n";
		my @entries       = &List::ExtraUtils::as_array($node->{'directive'}->{'entry'});
		my $lengthMaximum = max map {length($_->{'label'})} @entries;
		for(my $i=0;$i<scalar(@entries);++$i) {
		    my $entry     = $entries[$i];
		    my $separator = $i == scalar(@entries)-1 ? "." : ";";
		    $description .= "    description=description//char(10)//\"   ".(" " x ($lengthMaximum-length($entry->{'label'}))).$entry->{'label'}.(exists($entry->{'description'}) ? ": ".$entry->{'description'}.$separator : "")."\"\n";
		}
		$description        .= "    \n";
		$description        .= "    return\n";
		$descriptorFunction .= "  end function ".$functionName."\n";
		$descriptorFunction .= "  ! End auto-generated enumeration function\n";
		# Insert into the module.
		my $validatorTree = &Galacticus::Build::SourceTree::ParseCode($descriptorFunction,"Galacticus::Build::SourceTree::Process::Enumeration()", instrument => 0);
		my @validatorNodes = &Galacticus::Build::SourceTree::Children($validatorTree);
		my $newNode = $validatorNodes[0];
		while ( $newNode->{'type'} ne "function" ) {
		    $newNode = $newNode->{'sibling'};
		}
		$newNode = $newNode->{'firstChild'};
		while ( defined($newNode->{'sibling'}) ) {
		    $newNode = $newNode->{'sibling'};
		}
		my $describeNode =
		{
		    type => "code",
		    content => $description
		};
		&Galacticus::Build::SourceTree::InsertAfterNode($newNode,[$describeNode]);
		&Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@validatorNodes);
		# Set the visibility.
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$functionName,$visibility);
	    }
	    # Create documentation.
	    system("mkdir -p doc/enumerations/definitions");
	    open(my $defHndl,">doc/enumerations/definitions/".$node->{'directive'}->{'name'}.".tex");
	    print $defHndl "\\subsection{\\large {\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."}}\\hypertarget{ht:AutoEnumerations".ucfirst($node->{'directive'}->{'name'})."}{}\\label{sec:AutoEnumerations".ucfirst($node->{'directive'}->{'name'})."}\\index{enumerations!".$node->{'directive'}->{'name'}."\@{\\normalfont \\ttfamily ".$node->{'directive'}->{'name'}."}}\n\n";
	    print $defHndl "\\begin{tabular}{rp{130mm}}\n";
	    print $defHndl "Description: & ".$node->{'directive'}->{'description'}." \\\\\n";
	    my $moduleNode = $node;
	    do {
		$moduleNode = $moduleNode->{'parent'};
	    } until ( ! $moduleNode ||  $moduleNode->{'type'} eq "module" );
	    if ( $moduleNode ) {
		(my $fileName = $moduleNode->{'parent'}->{'name'}) =~ s/\./_/g;
		print $defHndl "Provided by: & {\\normalfont \\ttfamily module} \\href{https://github.com/galacticusorg/galacticus/releases/download/masterRelease/Galacticus_Source.pdf\\#source.".$fileName.":".lc($moduleNode->{'name'})."}{\\normalfont \\ttfamily ".latex_encode($moduleNode->{'name'})."} \\\\\n";
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
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
