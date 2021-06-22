#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use Galacticus::Options;
use GnuPlot::PrettyPlots;
use Class::Date qw(date now);

# Parse profiling information from a build log file.
# Andrew Benson (21-June-2021)

# Generates an HTML representation of the tasks being performed during each second of compilation.

# A color bar is shown next to each task indicating the period over which it was running.

# The color of the bar changes from red to green to indicate the number of simultaneous tasks running, such that red indicates
# regions where build parallelism is minimal.

# Get arguments.
die("Usage: buildProfiler.pl <buildLogFile> <profileFile> [options]")
    unless ( scalar(@ARGV) >= 2 );
my $buildLogFileName = $ARGV[0];
my $profileFileName  = $ARGV[1];
my %options = (
    durationMinimum => 1 # The shortest duration task which will be included in the profile.
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse profiling information from the file.
my @tasks;
my $timeEarliest;
my $timeLatest;
open(my $file,$buildLogFileName);
while ( my $line = <$file> ) {
    if ( $line =~ m/^\+\+Task: \{([\s\d\+\-:]+)\|([\s\d\+\-:]+)\} '(.*)/ ) {
	my $start     = $1;
	my $stop      = $2;
	my $command   = $3;
	my $startTime = date($start);
	my $stopTime  = date($stop );
	my $duration  = $stopTime-$startTime;
	# Process only tasks above the minimum duration.
	if ( $duration >= $options{'durationMinimum'} ) {
	    # Find earliest and latest times in the build.
	    if ( defined($timeEarliest) ) {
		$timeEarliest = $startTime
		    if ( $startTime < $timeEarliest );
	    } else {
		$timeEarliest = $startTime;
	    }
	    if ( defined($timeLatest) ) {
		$timeLatest = $stopTime
		    if ( $stopTime > $timeLatest );
	    } else {
		$timeLatest = $stopTime;
	    }
	    # Extract the command being run and convert to a more human-readable form where possible.
	    my @elements = split(" ",$command);
	    if ( $elements[0] eq "./scripts/build/preprocess.pl" ) {
		$elements[1] =~ s/source\///;
		$command = $elements[1]." (preprocess)";
	    } elsif ( $elements[0] eq "gfortran" ) {
		$elements[4] =~ s/\.\/work\/build\///;
		$command = $elements[4]." (compile)";
	    } elsif ( $elements[0] eq "./scripts/build/sourceDigests.pl" ) {
		$command = "source digests";
	    } elsif ( $elements[0] eq "./scripts/build/buildCode.pl" ) {
		$elements[2] =~ s/\.\/work\/build\///;
		$elements[2] =~ s/'//;
		$command = $elements[2]." (build)";
	    } elsif ( $elements[1] eq "-MRegexp::Common" ) {
		$elements[8] =~ s/work\/build\///;
		$command = $elements[8]." (cpp)";
	    } elsif ( $elements[0] =~ m/\.\/scripts\/build\/(.*)\.pl/ ) {
		$command = $1;
	    }
	    # Store this task to the list.
	    push(
		@tasks,
		{
		    description => $command,
		    startTime   => $startTime,
		    endTime     => $stopTime
		}
		);
	}
    }
}
close($file);

# Find the maximum time (in seconds) for the build.
my $timeMaximum = 0;
foreach my $task ( @tasks ) {
    $task->{'start'} = $task->{'startTime'}-$timeEarliest;
    $task->{'end'  } = $task->{'endTime'  }-$timeEarliest;
    $timeMaximum = $task->{'end'}
    if ( $task->{'end'} > $timeMaximum );    
}

# Find the number of threads running during each second of the build.
my $threadCountMaximum = 0;
my @threadCount;
for(my $i=0;$i<=$timeMaximum;++$i) {
    $threadCount[$i] = 0;
    foreach my $task ( @tasks ) {
	++$threadCount[$i]
	    if ( $i >= $task->{'start'} && $i <= $task->{'end'} );
    }
    $threadCountMaximum = $threadCount[$i]
	if ( $threadCount[$i] > $threadCountMaximum );
}

# Compute a cost for each task. This is the sum of the inverse of the number of threads running each second.
my $costMaximum = 0.0;
my $taskMaximum;
foreach my $task ( @tasks ) {
    $task->{'cost'} = 0.0;
    $task->{'isMaximumCost'} = 0;
    for(my $i=$task->{'start'};$i<=$task->{'end'};++$i) {
	$task->{'cost'} += 1.0/$threadCount[$i];
    }
    if ( $task->{'cost'} > $costMaximum ) {
	$costMaximum = $task->{'cost'};
	$taskMaximum = $task;
    }
}
$taskMaximum->{'isMaximumCost'} = 1;

# Create output.
open(my $profileFile,">",$profileFileName);
print $profileFile "<html>\n";
print $profileFile " <head>\n";
print $profileFile "  <style>\n";
print $profileFile "   body { font:16px Calibri;}\n";
print $profileFile "   td, th {\n";
print $profileFile "           margin: 0;\n";
print $profileFile "           border: 1px;\n";
print $profileFile "           height: 1em;\n";
print $profileFile "           white-space: nowrap;\n";
print $profileFile "           border-top-width: 0px;\n";
print $profileFile "           width: 1px;\n";
print $profileFile "          }\n";
print $profileFile "   div {\n";
print $profileFile "        width: 500px;\n";
print $profileFile "        overflow-x: scroll;\n";
print $profileFile "        margin-left: 22em;\n";
print $profileFile "        overflow-y: visible;\n";
print $profileFile "        padding-bottom: 1px;\n";
print $profileFile "       }\n";
print $profileFile "   .headcol {\n";
print $profileFile "             position: absolute;\n";
print $profileFile "             width: 0em;\n";
print $profileFile "             left: 0;\n";
print $profileFile "             top: auto;\n";
print $profileFile "             border-top-width: 1px;\n";
print $profileFile "             margin-top: -1px;\n";
print $profileFile "            }\n";
print $profileFile "   .headcolbold {\n";
print $profileFile "             position: absolute;\n";
print $profileFile "             width: 0em;\n";
print $profileFile "             left: 0;\n";
print $profileFile "             top: auto;\n";
print $profileFile "             border-top-width: 1px;\n";
print $profileFile "             margin-top: -1px;\n";
print $profileFile "             color: red;\n";
print $profileFile "            }\n";
print $profileFile "  </style>\n";
print $profileFile " </head>\n";
print $profileFile "<body>\n";
print $profileFile "Total compile time = ".$timeMaximum." seconds<p>\n";
print $profileFile "<div>\n";
print $profileFile "<table style=\"border-spacing: 0; border-collapse: separate; border-top: 1px\">\n";
foreach my $task ( @tasks ) {
    print $profileFile "<tr>";
    print $profileFile "<th class=\"headcol".($task->{'isMaximumCost'} ? "bold" : "")."\">".$task->{'description'}."</th>";
    for(my $i=0;$i<=$timeMaximum;++$i) {
	if ( $i < $task->{'start'}  ) {
	    print $profileFile "<td></td>";
	} elsif ( $i > $task->{'end'} ) {
	   # Do nothing here - missing columns at the end of the table row are ignored.
	} else {
	    my $fraction = $threadCount[$i]/$threadCountMaximum;
	    my $color    = &GnuPlot::PrettyPlots::gradientColor($fraction,$GnuPlot::PrettyPlots::colorGradients{'redGreen'});
	    print $profileFile "<td style=\"background: ".$color."\"></td>";
	}
    }
    print $profileFile "</tr>\n";
}
print $profileFile "</table>\n";
print $profileFile "</div>\n";
print $profileFile "</body>\n";
print $profileFile "</html>\n";
close($profileFile);

exit;
