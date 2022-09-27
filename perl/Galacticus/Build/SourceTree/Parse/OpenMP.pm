# Contains a Perl module which implements parsing of OpenMP directives.

package Galacticus::Build::SourceTree::Parse::OpenMP;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Fortran::Utils;
use Clone 'clone';
use Regexp::Common;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks{'OpenMP'} = \&Parse_OpenMP;

sub Parse_OpenMP {
    # Get the tree.
    my $tree = shift();
    # Build a balanced parenthesis matcher.
    my $balanced = $RE{balanced}{-parens=>'()'};
    # List of options which attach to iterators.
    my @iteratorOptions = ( "schedule" );
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Find code blocks.
	if ( $node->{'type'} eq "code" ) {
	    # Read the code block.
	    my $rawCode;
	    my $moduleUses;
	    my $lineNumber       = exists($node->{'line'  }) ? $node->{'line'  } : 0        ;
	    my $source           = exists($node->{'source'}) ? $node->{'source'} : "unknown";
	    my @newNodesCombined;
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran::Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is an OpenMP parallel directive.
		if ( $processedLine =~ m/^\s*!\$omp\s+(end\s+)?parallel\s+(do|workshare)?\s?(.*)/ ) {
		    my $isCloser    = defined($1);
		    my $composite   =         $2 ;
		    my $options     =         $3 ;
		    $options        =~ s/^\s*//;
		    my @ompOptions;
		    while ( $options ne "" ) {
			if ( $options =~ m/^([a-z]+)(.*)/ ) {
			    my $name      = $1;
			    my $remainder = $2;
			    my @contents  = $remainder =~ /$balanced/;
			    my $content   = scalar(@contents) > 0 ? $contents[0] : "";
			    $options =~ s/^$name\s*//;
			    $options =~ s/^\Q$content//;
			    $options =~ s/^\s*//;
			    push(
				@ompOptions,
				{
				    name    => $name   ,
				    content => $content
				}
				);
			}
		    }
		    # Generate new code.
		    my $codeParallel;
		    my $codeComposite;
		    if ( defined($composite) ) {
			$codeParallel  = "!\$omp ".($isCloser ? "end " : "")."parallel";
			$codeComposite = "!\$omp ".($isCloser ? "end " : "").$composite;
			foreach my $ompOption ( @ompOptions ) {
			    my $codeOption = " ".$ompOption->{'name'}.$ompOption->{'content'};
			    if ( grep {$_ eq $ompOption->{'name'}} @iteratorOptions ) {
				$codeComposite .= $codeOption;
			    } else {
				$codeParallel  .= $codeOption;
			    }
			}
			$codeParallel  .= "\n";
			$codeComposite .= "\n";
		    } else {
			$codeParallel = $rawLine;
		    }
		    my @newNodes;
		    my @ompOptionsParallel;
		    foreach my $ompOption ( @ompOptions ) {
			push(@ompOptionsParallel,$ompOption)
			    unless ( grep {$_ eq $ompOption->{'name'}} @iteratorOptions );
		    }
		    my $nodeParallel =
		    {
			type     => "openMP"            ,
			name     => "parallel"          ,
			isCloser => $isCloser           ,
			options  => \@ompOptionsParallel
		    };	
		    $nodeParallel->{'firstChild'} =
		    {
		    	type       => "code"       ,
		    	content    => $codeParallel,
		    	parent     => $nodeParallel,
		    	sibling    => undef()      ,
		    	firstChild => undef()      ,
		    	source     => $source      ,
		    	line       => $lineNumber
		    };
		    push(@newNodes,$nodeParallel);
		    if ( defined($composite) ) {
			my @ompOptionsComposite;
			foreach my $ompOption ( @ompOptions ) {
			    push(@ompOptionsComposite,$ompOption)
				if ( grep {$_ eq $ompOption->{'name'}} @iteratorOptions );
			}
			my $nodeComposite =
			{
			    type     => "openMP"             ,
			    name     => $composite           ,
			    isCloser => $isCloser            ,
			    options  => \@ompOptionsComposite
			};	
			$nodeComposite->{'firstChild'} =
			{
			    type       => "code"         ,
			    content    => $codeComposite ,
			    parent     => $nodeComposite ,
			    sibling    => undef()        ,
			    firstChild => undef()        ,
			    source     => $source        ,
			    line       => $lineNumber
			};
			if ( $isCloser ) {
			    unshift(@newNodes,$nodeComposite);
			} else {
			    push   (@newNodes,$nodeComposite);
			}
		    }
		    if ( defined($rawCode) ) {
			my $rawNode =
			{
			    type    => "code",
			    content => $rawCode
			};
			undef($rawCode);
			unshift(@newNodes,$rawNode);
		    }
		    push(@newNodesCombined,@newNodes);
		} else {
		    $rawCode .= $rawLine;
		}
		$lineNumber += $rawLine =~ tr/\n//;
	    } until ( eof($code) );
	    close($code);
	    if ( defined($rawCode) ) {
		my $rawNode =
		{
		    type    => "code",
		    content => $rawCode
		};
		undef($rawCode);
		push(@newNodesCombined,$rawNode);
	    }
	    &Galacticus::Build::SourceTree::ReplaceNode($node,\@newNodesCombined); 
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }    
}

sub update {    
    # Update the contained code.
    my $node = shift();
    # Create the new content.
    $node->{'firstChild'}->{'content'} = "!\$omp ".($node->{'isCloser'} ? "end " : "").$node->{'name'}.join("",map {" ".$_->{'name'}.$_->{'content'}} @{$node->{'options'}})."\n";
}

sub copyin {
    # Add copyin directives for named variables.
    my $node          =   shift() ;
    my @variableNames = @{shift()};
    # If no copyin option exists, add one.
    push(
	@{$node->{'options'}},
	{
	    name    => "copyin",
	    content => "()"
	}
	)
	unless ( grep {$_->{'name'} eq "copyin"} @{$node->{'options'}} );
    # Extract a list of current names.
    my @option = grep {$_->{'name'} eq "copyin"} @{$node->{'options'}};
    my $content = $option[0]->{'content'};
    $content =~ s/^\((.*)\)$/$1/;
    my @names = split(/\s*,\s*/,$content);
    # Insert new names.
    foreach my $name ( @variableNames ) {
	push(@names,$name)
	    unless ( grep {$_ eq $name} @names );
    }
    # Update the option.
    $option[0]->{'content'} = "(".join(",",@names).")";
    # Update the code content.
    &update($node);
}

1;
