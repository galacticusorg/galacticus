# Contains a Perl module which implements profiling tools for OpenMP.

package Galacticus::Build::SourceTree::Process::ProfileOpenMP;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'profileOpenMP'} = \&Profile_OpenMP;

sub Profile_OpenMP {
    # Get the tree.
    my $tree = shift();
    # Check if debugging is required.
    my $profile = 0;
    if ( exists($ENV{'GALACTICUS_FCFLAGS'}) ) {
	$profile = 1
	    if ( grep {$_ eq "-DOMPPROFILE"} split(" ",$ENV{'GALACTICUS_FCFLAGS'}) );
    }
    return
	unless ( $profile );
    # Load critical section name enumeration.
    die("Galacticus::Build::SourceTree::Process::ProfileOpenMP::Profile_OpenMP(): critical section enumeration file does not exist")
	unless ( -e $ENV{'BUILDPATH'}."/openMPCriticalSections.xml" );
    my $xml        = new XML::Simple();
    my $descriptor = $xml->XMLin($ENV{'BUILDPATH'}."/openMPCriticalSections.xml");
    # Modules required.
    my $moduleUses =
    {
	type      => "moduleUse",
	moduleUse =>
	{
	    OpenMP_Utilities_Data => {intrinsic => 0, all => 1},
	    OMP_Lib               => {intrinsic => 0, all => 1}
	}
    };
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    my $newContent;
	    open(my $content,"<",\$node->{'content'});
	    while ( my $line = <$content> ) {
		if ( $line =~ m/^\s*\!\$omp\s+critical\s*\(([a-z0-9_]+)\)/i ) {
		    my $criticalSectionName = lc($1);
		    if ( exists($descriptor->{'critical'}->{$criticalSectionName}) ) {
			# Find parent function.
			my $nodeParent = $node;
			while ( $nodeParent->{'type'} ne "function" && $nodeParent->{'type'} ne "subroutine" ) {
			    $nodeParent = $nodeParent->{'parent'};
			}
			# Add modules required.
			&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($nodeParent,$moduleUses);
			# Generate code to time the wait at this critical section.
			$line = "ompProfileTimeWaitStart=OMP_Get_WTime()\n".$line."ompProfileTimeWaitEnd=OMP_Get_WTime()\nompProfileTimeWaitEnd=ompProfileTimeWaitEnd-ompProfileTimeWaitStart\ncriticalSectionWaitTime(".$descriptor->{'critical'}->{$criticalSectionName}->{'id'}.")=criticalSectionWaitTime(".$descriptor->{'critical'}->{$criticalSectionName}->{'id'}.")+ompProfileTimeWaitEnd\n";
		    }
		}
		$newContent .= $line;
	    }
	    close($content);
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
