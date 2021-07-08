# Contains a Perl module which implements processing of "functionCall" directives in the Galacticus build system.

package Galacticus::Build::FunctionCall;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Scalar::Util 'reftype';
use Galacticus::Build::Hooks;
use Galacticus::Build::Dependencies;

# Insert hooks for our functions.
%Galacticus::Build::Hooks::moduleHooks = 
    (
     %Galacticus::Build::Hooks::moduleHooks,
     functionCall => {parse => \&Function_Calls_Parse_Directive, generate => \&Function_Calls_Generate_Output}
    );

sub Function_Calls_Parse_Directive {
    # Parse content for a "functionCall" directive.
    my $buildData = shift;

    # Assert that we have a prefix, currentDocument and directive.
    die("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: no currentDocument present"  )
	unless ( exists($buildData->{'currentDocument'}               ) );
    die("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: no unitName present"         )
	unless ( exists($buildData->{'currentDocument'}->{'unitName'} ) );
    die("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: no functionType present"     )
	unless ( exists($buildData->{'functionType'   }               ) );

    # For procedure pointer cases, assert that we have the pointer name.
    if ( $buildData->{'functionType'} eq "pointer" ) {
	die("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: no pointerName present")
	    unless ( exists($buildData->{'pointerName'}) );
    }

    # Extract the unit name.
    my $unitName = $buildData->{'currentDocument'}->{'unitName'};

    # If a globalized-version is to be called, append the suffix.
    $unitName .= "_"
	if ( exists($buildData->{'currentDocument'}->{'useGlobal'}) && $buildData->{'currentDocument'}->{'useGlobal'} eq "yes" );
    
    # Determine the type of function call.
    my $include = 1;
    if ( $buildData->{'functionType'} eq "pointer" ) {
	# Action is to assign a procedure pointer.
	$buildData->{'functionCall'}->{$unitName}->{'code'} = $buildData->{'pointerName'}." => ".$unitName."\n";
    } elsif ( $buildData->{'functionType'} eq "void" ) {
     	# Action is to call a void function.
	if ( exists($buildData->{'exclude'}) ) {
	    $include = 0
		if ( $buildData->{'exclude'} eq $buildData->{'currentDocument'}->{'unitName'} );
	}
	if ( $include == 1 ) {
	    $buildData->{'functionCall'}->{$unitName}->{'code'} = "call ".$unitName."(";
	    # If we have function arguments, append them to the call.
	    my $arguments = "";
	    if ( exists($buildData->{'functionArgs'}) ) {
		if ( defined(reftype($buildData->{'functionArgs'})) && reftype($buildData->{'functionArgs'}) eq "HASH" ) {
		    if ( exists($buildData->{'functionArgs'}->{$buildData->{'codeType'}}) ) {
			$arguments = $buildData->{'functionArgs'}->{$buildData->{'codeType'}};
		    } else {
			die ("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: function arguments not given for this language")
		    }
		} else {
		    die ("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: function arguments not given for this language")
			unless ( $buildData->{'codeType'} eq "fortran" );
		    $arguments = $buildData->{'functionArgs'};
		}
	    }
	    # Add arguments.
	    unless ( $arguments eq "" ) {
		$buildData->{'functionCall'}->{$unitName}->{'code'} .= ","
		    unless ( $buildData->{'functionCall'}->{$unitName}->{'code'} =~ m/\($/ );
		$buildData->{'functionCall'}->{$unitName}->{'code'} .= $arguments;
	    }
	    $buildData->{'functionCall'}->{$unitName}->{'code'} .= ")\n";	       
	}
    } elsif ( $buildData->{'functionType'} eq "function" ) {
     	# Action is to call a function that returns a result.
	if ( exists($buildData->{'exclude'}) ) {
	    $include = 0
		if ( $buildData->{'exclude'} eq $buildData->{'currentDocument'}->{'unitName'} );
	}
	if ( $include == 1 ) {
	    die('Galacticus::Build::FunctionCall: no return parameter specified')
		unless ( exists($buildData->{'returnParameter'}) );
	    $buildData->{'functionCall'}->{$unitName}->{'code'} = $buildData->{'returnParameter'}."=".$buildData->{'currentDocument'}->{'unitName'}."(";
	    # If we have function arguments, append them to the call.
	    my $arguments = "";
	    if ( exists($buildData->{'functionArgs'}) ) {
		if ( defined(reftype($buildData->{'functionArgs'})) && reftype($buildData->{'functionArgs'}) eq "HASH" ) {
		    if ( exists($buildData->{'functionArgs'}->{$buildData->{'codeType'}}) ) {
			$arguments = $buildData->{'functionArgs'}->{$buildData->{'codeType'}};
		    } else {
			die ("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: function arguments not given for this language")
		    }
		} else {
		    die ("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: function arguments not given for this language")
			unless ( $buildData->{'codeType'} eq "fortran" );
		    $arguments = $buildData->{'functionArgs'};
		}
	    }
	    # Add arguments.
	    unless ( $arguments eq "" ) {
		$buildData->{'functionCall'}->{$unitName}->{'code'} .= ","
		    unless ( $buildData->{'functionCall'}->{$unitName}->{'code'} =~ m/\($/ );
		$buildData->{'functionCall'}->{$unitName}->{'code'} .= $arguments;
	    }
	    $buildData->{'functionCall'}->{$unitName}->{'code'} .= ")\n";	       
	}
    }

    # Add in any dependency sort information.
    $buildData->{'functionCall'}->{$unitName}->{'after'   } = $buildData->{'currentDocument'}->{'after'   }
        if ( exists($buildData->{'currentDocument'}->{'after' }) );
    $buildData->{'functionCall'}->{$unitName}->{'before'  } = $buildData->{'currentDocument'}->{'before'  }
        if ( exists($buildData->{'currentDocument'}->{'before'}) );
    $buildData->{'functionCall'}->{$unitName}->{'sortName'} = $buildData->{'currentDocument'}->{'sortName'}
        if ( exists($buildData->{'currentDocument'}->{'sortName'}) );

    # If some action is specified, perform this action after the function call.
    $buildData->{'functionCall'}->{$unitName}->{'code'} .= $buildData->{'onReturn'}."\n"
	if ( exists($buildData->{'onReturn'}) && $include == 1 );

}

sub Function_Calls_Generate_Output {
    # Generate output for a "functionCall" directive.
    my $buildData = shift;

    # Assert that we have a file name and directive present.
    die("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: no fileName present"      )
	unless ( exists($buildData->{'fileName'}) );
    die("Galacticus::Build::FunctionCall::Function_Calls_Parse_Directive: no directive present")
	unless ( exists($buildData->{'directive'}) );

    # Add a header.
    $buildData->{'content'}  = "! Generated automatically by Galacticus::Build::FunctionCall\n";
    $buildData->{'content'} .= "!  From: ".$buildData->{'fileName'}."\n";

    # Perform the dependency sort.
    &Galacticus::Build::Dependencies::Dependency_Sort($buildData->{'functionCall'},$buildData);

    # Find the length of the longest unit name.
    my $longestNameLength = 0;
    foreach my $unitName ( @{$buildData->{'unitNames'}} ) {
	$longestNameLength = length($unitName)
	    if ( length($unitName) > $longestNameLength );
    }

    # Iterate over all labels, and add them to the content.
    foreach my $unitName ( @{$buildData->{'unitNames'}} ) {
	my $functionCall = $buildData->{'functionCall'}->{$unitName}->{'code'};
	# For calls to void functions, add padding to align the arguments.
	if ( $buildData->{'functionType'} eq "void" ) {
	    my $paddingLength = $longestNameLength-length($unitName);
	    my $padding       = " " x $paddingLength;
	    $functionCall     =~ s/^([^\(]+)/$1$padding/;
	}
	$buildData->{'content'} .= $functionCall;
    }

}

1;
