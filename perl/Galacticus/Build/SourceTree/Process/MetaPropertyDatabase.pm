# Contains a Perl module which implements a database of meta-properties.

package Galacticus::Build::SourceTree::Process::MetaPropertyDatabase;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use List::ExtraUtils;
use List::Util 'max';
use Fortran::Utils;
use Galacticus::Build::Directives;
use Text::Template 'fill_in_string';
use Digest::MD5 qw(md5_hex);
use Galacticus::Build::SourceTree::Process::SourceIntrospection;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'metaPropertyDatabase'} = \&Process_MetaPropertyDatabase;

sub Process_MetaPropertyDatabase {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get code directive locations.
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
    my $stateStorables     = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml"    );
    # Walk the tree, looking for hook directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Handle metaPropertyDatabase directives.
	if ( $node->{'type'} eq "metaPropertyDatabase" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Iterate over all files that contain an "addMetaProperty" directive.
	    my @creators;
	    foreach my $fileName ( &List::ExtraUtils::as_array($directiveLocations->{'addMetaProperty'}->{'file'}) ) {
		# Extract all "addMetaProperty" directives in this file which create their meta-property, and skip if no such exist.
		my @addMetaProperties = grep {exists($_->{'isCreator'}) && $_->{'isCreator'} eq "yes"} &Galacticus::Build::Directives::Extract_Directives($fileName,'addMetaProperty');
		next
		    unless ( scalar(@addMetaProperties) > 0 );
		# Determine the functionClass defined in this file.
		my $functionClassName;
		my $implementationName;
		foreach my $directive ( &Galacticus::Build::Directives::Extract_Directives($fileName,'*',setRootElementType => 1) ) {
		    if ( grep {$_ eq $directive->{'rootElementType'}."Class"} keys(%{$stateStorables->{'functionClasses'}}) ) {
			$functionClassName  = $directive->{'rootElementType'};
			$implementationName = $directive->{'name'           };
			if ( $implementationName =~ m/^$functionClassName/ ) {
			    $implementationName =~ s/^$functionClassName//;
			    $implementationName = lcfirst($implementationName);
			} else {
			    die("Galacticus::Build::SourceTree::Process::MetaPropertyDatabase::Process_MetaPropertyDatabase(): functionClass implementation name has incorrect prefix");
			}
			last;
		    }
		}
		die("Galacticus::Build::SourceTree::Process::MetaPropertyDatabase::Process_MetaPropertyDatabase(): unable to determine functionClass implementation")
		    unless ( defined($functionClassName) && defined($implementationName) );
		foreach my $addMetaProperty ( @addMetaProperties ) {
		    $addMetaProperty->{'functionClass'     } = $functionClassName;
		    $addMetaProperty->{'implementationName'} = $implementationName;
		    $addMetaProperty->{'type'              } = "float"
			unless ( exists($addMetaProperty->{'type'}) );
		    $addMetaProperty->{'rank'              } = 0
			unless ( exists($addMetaProperty->{'rank'}) );
		    push(@creators,$addMetaProperty);
		}
	    }
	    my $noCreator;
	    $noCreator = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine metaPropertyNoCreator(component_,name_,type_,rank_)
 use :: Error, only : Error_Report
 implicit none
 character(len=* ), intent(in   ) :: component_, name_             , type_
 integer          , intent(in   ) :: rank_
 character(len=64)                :: className , implementationName, rankLabel

 write (rankLabel,'(i1)') rank_
CODE
	    my $join = "";
	    foreach my $creator ( @creators ) {
		unless ( $creator->{'name'} =~ m/^'/ ) {
		    $noCreator .= $join."if (component_ == '".$creator->{'component'}."' .and. name_ == '".$creator->{'name'}."' .and. type_ == '".$creator->{'type'}."' .and. rank_ == ".$creator->{'rank'}.") then\n";
		    $noCreator .= "          className='".$creator->{'functionClass'     }."'\n";
		    $noCreator .= " implementationName='".$creator->{'implementationName'}."'\n";
		    $join       = "else ";
		}
	    }
	    $code::location = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'});
	    $noCreator .= fill_in_string(<<'CODE', PACKAGE => 'code');
 else
  className         =""
  implementationName=""
  call Error_Report("no class creates the rank-"//trim(rankLabel)//" '"//trim(type_)//"' type meta-property '"//trim(name_)//"' in component '"//trim(component_)//"'"//{$location})
 end if
 call Error_Report("the rank-"//trim(rankLabel)//" '"//trim(type_)//"' type meta-property '"//trim(name_)//"' in component '"//trim(component_)//"' is required"//char(10)//"it is created by the '"//trim(implementationName)//"' implementation of the '"//trim(className)//"' class"//char(10)//"to create this meta-property include the following in your parameter file:"//char(10)//" <"//trim(className)//" value="""//trim(implementationName)//"""/>"//{$location})
 return
end subroutine metaPropertyNoCreator
CODE
	    # Insert the code.
	    my $newNode =
	    {
		type       => "code",
		content    => $noCreator,
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
