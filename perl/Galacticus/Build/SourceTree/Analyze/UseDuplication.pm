# Contains a Perl module which analyzes pre-processed source code for duplication of "use" statements.

package Galacticus::Build::SourceTree::Analyze::UseDuplication;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use List::Uniq qw(uniq);
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::analyzeHooks{'useDuplication'} = \&Analyze_UseDuplication;

sub Analyze_UseDuplication {
    # Get the tree.
    my $tree  = shift();
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0;
    my $fileName;
    while ( $node ) {
	if ( $node->{'type'} eq "moduleUse" ) {
	    # Look for duplication of module uses in functions or subroutines only.
	    if ( $node->{'parent'}->{'type'} eq "function" || $node->{'parent'}->{'type'} eq "subroutine" ) {
		# Walk up the hierarchy finding already used modules.
		my @usedModules = ();
		my $nodeParent = $node->{'parent'}->{'parent'};
		while ( $nodeParent ) {
		    my $nodeChild = $nodeParent->{'firstChild'};
		    while ( $nodeChild ) {
			if ( $nodeChild->{'type'} eq "moduleUse" ) {
			    push(@usedModules,sort(keys(%{$nodeChild->{'moduleUse'}})));
			}
			$nodeChild = $nodeChild->{'sibling'};
		    }
		    $nodeParent = $nodeParent->{'parent'};
		}
		foreach my $moduleName ( sort(keys(%{$node->{'moduleUse'}})) ) {
		    if ( grep {$moduleName eq $_} @usedModules ) {
			print "Warning: module '".$moduleName."' used in ".$node->{'parent'}->{'name'}."() was already used in container [file: ".$fileName."]\n";
		    }
		}
	    }
	    # Determine symbols used from each module.
	    foreach my $module ( sort(keys(%{$node->{'moduleUse'}})) ) {
		my @symbolsExported = &Fortran::Utils::moduleSymbols(lc($module));
		my @symbolsUsed;
		# Search through code looking for use of these symbols.
		my $nodeCode  = $node->{'parent'};
		my $depthCode = 0;
		while ( $nodeCode->{'firstChild'} ) {
		    $nodeCode = $nodeCode->{'firstChild'};
		    ++$depthCode;
		}
		while ( $nodeCode && $depthCode >= 0 ) {
		    if ( $nodeCode->{'type'} eq "code" ) {
			open(my $code,"<",\$nodeCode->{'content'});
			do {
			    my $rawLine;
			    my $line;
			    my $comments;
			    &Fortran::Utils::Get_Fortran_Line($code,$rawLine,$line,$comments);
			    foreach my $symbol ( @symbolsExported ) {
				push(@symbolsUsed,$symbol)
				    if ( $line =~ m/\b$symbol\b/i );
			    }
			} until ( eof($code) );
			close($code);
		    }
		    $nodeCode = &Galacticus::Build::SourceTree::Walk_Tree($nodeCode,\$depthCode);
		}
		@symbolsUsed = uniq(sort(@symbolsUsed));
		# Detect cases where no symbol is used from a module.
		if ( scalar(@symbolsUsed) == 0 ) {
		    print "Warning: no symbols were used from module '".$module."' in ".$node->{'parent'}->{'type'}." '".$node->{'parent'}->{'name'}."' [file: ".$fileName."]\n";
		} else {
		    # If the module imported all symbols, report only those which were actually used.
		    if ( exists($node->{'moduleUse'}->{$module}->{'all'}) ) {
			print "Warning: ".$node->{'parent'}->{'type'}." '".$node->{'parent'}->{'name'}."' imports all symbols from '".$module."' but only uses\n";
			print join("\n",map {"           ".$_} @symbolsUsed)."\n";
		    } else {
			# Only specified symbols were imported. Check for those which were not needed.
			my @symbolsNonRequired;
			foreach my $symbolImported ( sort(keys(%{$node->{'moduleUse'}->{$module}->{'only'}})) ) {
			    push(@symbolsNonRequired,$symbolImported)
				unless ( grep {$_ eq $symbolImported} @symbolsUsed );
			}
			if ( scalar(@symbolsNonRequired) > 0 ) {
			    print "Warning: ".$node->{'parent'}->{'type'}." '".$node->{'parent'}->{'name'}."' imports these symbols from '".$module."' which are not used\n";
			    print join("\n",map {"           ".$_} @symbolsNonRequired)."\n";
			}
		    }
		}
	    }
	} elsif ( $node->{'type'} eq "file" ) {
	    $fileName = $node->{'name'};
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
