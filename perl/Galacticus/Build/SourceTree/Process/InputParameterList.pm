# Contains a Perl module which implements constructing a list of parameter names used in a unit.

package Galacticus::Build::SourceTree::Process::InputParameterList;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use XML::Simple;
use LaTeX::Encode;
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Directives;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'inputParameterList'} = \&Process_InputParameterList;

sub Process_InputParameterList {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get the file name.
    my $fileName;
    $fileName = $tree->{'name'}
        if ( $tree->{'type'} eq "file" );
    # Get code directive locations.
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
    # Initialize list of unlisted parameters.
    my @unlistedInputParameters;
    # Walk the tree, looking for input parameter list directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Look for inputParameterList directives and process them.
	if ( $node->{'type'} eq "inputParameterList" && ! $node->{'directive'}->{'processed'} ) {	    
	    # Record that this directive has been processed.
	    $node->{'directive'}->{'processed'} =  1;
	    # Determine the source name.
	    my $sourceMatch;
	    $sourceMatch = $node->{'directive'}->{'source'}
	        if ( exists($node->{'directive'}->{'source'}) );
	    # Step through sibling nodes looking for input parameter directives.
	    my @inputParameterNames;
	    my $sibling = $node->{'parent'}->{'firstChild'};
	    while ( $sibling ) {
		if ( $sibling->{'type'} eq "inputParameter" ) {
		    my $source = "globalParameters";
		    $source = $sibling->{'directive'}->{'source'}
		        if ( exists($sibling->{'directive'}->{'source'}) );
		    if ( ! $sourceMatch || $sourceMatch eq $source ) {
			if      ( exists($sibling->{'directive'}->{'name'}) ) {
			    # Single parameter defined by its name - simply push onto the list.
			    push(@inputParameterNames,$sibling->{'directive'}->{'name'});
			} elsif ( exists($sibling->{'directive'}->{'iterator'}) ) {
			    # A parameter whose name iterates over a set of possible names.
			    if ( $sibling->{'directive'}->{'iterator'} =~ m/\(\#([a-zA-Z0-9]+)\-\>([a-zA-Z0-9]+)\)/ ) {
				my $directiveName = $1;
				my $attributeName = $2;
				die('Process_InputParameterList(): locations not found for directives')
				    unless ( exists($directiveLocations->{$directiveName}) );
				foreach my $fileName ( &List::ExtraUtils::as_array($directiveLocations->{$directiveName}->{'file'}) ) {
				    foreach ( &Galacticus::Build::Directives::Extract_Directives($fileName,$directiveName) ) {
					(my $parameterName = $sibling->{'directive'}->{'iterator'}) =~ s/\(\#$directiveName\-\>$attributeName\)/$_->{$attributeName}/;
					push(@inputParameterNames,$parameterName);
				    }
				}
			    } else {
				die('Process_InputParameterList(): nothing to iterate over');
			    }
			}
		    }
		} elsif ( $sibling->{'type'} eq "objectBuilder" ) {
		    # Add methods read by objectBuilder directives.
		    push(@inputParameterNames,$sibling->{'directive'}->{'class'}."Method")
			if ( ! $sourceMatch || $sourceMatch eq $sibling->{'directive'}->{'source'} );
		}
		$sibling = $sibling->{'sibling'};
	    }
	    # Generate the variable declaration.
	    my $declaration =
	    {
		intrinsic  => "type",
		type       => "varying_string",
		attributes => [ "dimension(".scalar(@inputParameterNames).")" ],
		variables  => [ $node->{'directive'}->{'label'} ]
	    };
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},[$declaration]);
	    # Add module usage.
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},{moduleUse => {ISO_Varying_String => {all => 1}}});
	    # Generate the setting code.
	    if ( scalar(@inputParameterNames) > 0 ) {
		my $setter;
		for(my $i=1;$i<=scalar(@inputParameterNames);++$i) {
		    $setter .= $node->{'directive'}->{'label'}."(".$i.")='".$inputParameterNames[$i-1]."'\n";
		}
		# Insert new code.
		my $setterNode =
		{
		    type       => "code"           ,
		    content    => $setter          ,
		    sibling    => undef()          ,
		    parent     => $node->{'parent'},
		    firstChild => undef()
		};
		&Galacticus::Build::SourceTree::InsertAfterNode($node,[$setterNode]);
	    }
	}
	# Look for inputParameter directives for which no input parameter list is defined.
	if ( $node->{'type'} eq "inputParameter" && $fileName ) {
	    my $sibling = $node->{'parent'}->{'firstChild'};
	    my $inList  = 0;
	    while ( $sibling ) {
		if ( $sibling->{'type'} eq "inputParameterList" ) {
		    $inList = 1;
		    last;
		}
		$sibling = $sibling->{'sibling'};
	    }
	    unless ( $inList ) {
		if ( exists($node->{'directive'}->{'name'}) ) {
		    push(@unlistedInputParameters,         $node->{'directive'}->{'name' });
		} elsif ( exists($node->{'directive'}->{'regEx'}) ) {
		    push(@unlistedInputParameters,"regEx:".$node->{'directive'}->{'regEx'});
		}
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Output file of unlisted parameters.
    if ( @unlistedInputParameters ) {
	$fileName =~ s/\.F90$/.p/;
	open(my $parametersFile,">>".$ENV{'BUILDPATH'}."/".$fileName);
	print $parametersFile $_."\n"
	    foreach ( @unlistedInputParameters );
	close($parametersFile);
    }
    
}

1;
