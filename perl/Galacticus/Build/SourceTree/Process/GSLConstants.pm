# Contains a Perl module which implements extraction of constants from GSL.

package Galacticus::Build::SourceTree::Process::GSLConstants;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Temp qw/ tempfile /;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'gslConstant'} = \&Process_GSLConstant;

sub Process_GSLConstant {
    # Get the tree.
    my $tree = shift();
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "gslConstant" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate, compile, and execute a minimal C code which simply outputs the relevant GSL constant.
	    my $type = exists($node->{'directive'}->{'type'}) ? $node->{'directive'}->{'type'} : "double precision";
	    my $code;
	    if ( $type eq "enum" ) {
		# An enum type.
		my @members = split(/\s*,\s*/,$node->{'directive'}->{'members'});
		(undef, my $tmpSourceFileName) = tempfile('tempXXXXX', SUFFIX => '.c', OPEN => 0, DIR => $ENV{'BUILDPATH'});
		(undef, my $tmpExecFileName  ) = tempfile('tempXXXXX', SUFFIX => '.x', OPEN => 0, DIR => $ENV{'BUILDPATH'});
		open(my $tmpSource,">".$tmpSourceFileName);
		print $tmpSource "\#include <stdio.h>\n";
		print $tmpSource "\#include <float.h>\n";
		print $tmpSource "\#include <gsl/".$node->{'directive'}->{'gslHeader'}.".h>\n";
		print $tmpSource "int main () {\n";
		print $tmpSource " printf(\"".join(" ",("%i") x scalar(@members))."\", ".$node->{'directive'}->{'members'}.");\n";
		print $tmpSource "}\n";	
		close($tmpSource);
		system($ENV{'CCOMPILER'}." -o ".$tmpExecFileName." ".$tmpSourceFileName." ".$ENV{'CFLAGS'});
		die("Galacticus::Build::SourceTree::Process::GSLConstants: failed to compile")
		    unless ( $? == 0 );
		open(my $gslPipe,$tmpExecFileName."|");
		my $gslEnum = <$gslPipe>;
		close($gslPipe);
		unlink($tmpSourceFileName,$tmpExecFileName);
		my @values = split(" ",$gslEnum);
		my @definitions;
		for(my $i=0;$i<scalar(@members);++$i) {
		    push(@definitions,$members[$i]."=".$values[$i]);
		}
		$code  = "enum, bind(c)\n";
		$code .= " enumerator :: ".join(", ",@definitions)."\n"; 
		$code .= "end enum\n";
	    } else {
		# A non-enum type.
		# For double precision constants we use a macro to define the precision with which the constant should be output to
		# preserve full precision (see https://stackoverflow.com/a/19897395).
		(undef, my $tmpSourceFileName) = tempfile('tempXXXXX', SUFFIX => '.c', OPEN => 0, DIR => $ENV{'BUILDPATH'});
		(undef, my $tmpExecFileName  ) = tempfile('tempXXXXX', SUFFIX => '.x', OPEN => 0, DIR => $ENV{'BUILDPATH'});
		open(my $tmpSource,">".$tmpSourceFileName);
		print $tmpSource "\#include <stdio.h>\n";
		print $tmpSource "\#include <float.h>\n";
		print $tmpSource "\#include <gsl/".$node->{'directive'}->{'gslHeader'}.".h>\n";
		if ( $type eq "double precision" ) {
		    print $tmpSource "\#ifdef LDBL_DECIMAL_DIG\n";
		    print $tmpSource "  \#define OP_LDBL_Digs (LDBL_DECIMAL_DIG)\n";
		    print $tmpSource "\#else  \n";
		    print $tmpSource "  \#ifdef DECIMAL_DIG\n";
		    print $tmpSource "    \#define OP_LDBL_Digs (DECIMAL_DIG)\n";
		    print $tmpSource "  \#else  \n";
		    print $tmpSource "    \#define OP_LDBL_Digs (LDBL_DIG + 3)\n";
		    print $tmpSource "  \#endif\n";
		    print $tmpSource "\#endif\n";
		    print $tmpSource "int main () {\n";
		    print $tmpSource " printf(\"%.*e\", OP_LDBL_Digs - 1, ".$node->{'directive'}->{'gslSymbol'}.");\n";
		    print $tmpSource "}\n";
		} elsif ( $type eq "integer" ) {
		    print $tmpSource "int main () {\n";
		    print $tmpSource " printf(\"%i\", ".$node->{'directive'}->{'gslSymbol'}.");\n";
		    print $tmpSource "}\n";	
		}
		close($tmpSource);
		system($ENV{'CCOMPILER'}." -o ".$tmpExecFileName." ".$tmpSourceFileName." ".$ENV{'CFLAGS'});
		die("Galacticus::Build::SourceTree::Process::GSLConstants: failed to compile")
		    unless ( $? == 0 );
		open(my $gslPipe,$tmpExecFileName."|");
		my $gslConstant = <$gslPipe>;
		close($gslPipe);
		unlink($tmpSourceFileName,$tmpExecFileName);
		$gslConstant =~ s/e/d/;
		$code = $type.", parameter, public :: ".$node->{'directive'}->{'variable'}."=".$gslConstant."\n";
	    }
	    my $newNode =
	    {
		type       => "code",
		content    => $code,
		firstChild => undef()
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
