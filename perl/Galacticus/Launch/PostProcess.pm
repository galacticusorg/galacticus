# Postprocess models.

package Galacticus::Launch::PostProcess;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Copy;
use File::Slurp;
use IO::Compress::Simple;
use System::Redirect;

sub Failed {
    # The run failed for some reason.
    my $job          = shift();
    my $launchScript = shift();
    # Move any core file to the output directory.
    opendir(gDir,".");
    while ( my $file = readdir(gDir) ) {
	move($file,$job->{'directory'}."/core")
	    if ( $file =~ m/core\.\d+/ );
    }
    closedir(gDir);
    # Report the model failure (by e-mail if we have an e-mail address to send a report to and if
    # so requested).
    my $message = "FAILED: A Galacticus model failed to finish:\n\n";
    $message   .= "  Host:\t".$ENV{"HOSTNAME"}."\n";
    $message   .= "  User:\t".$ENV{"USER"}."\n\n";
    $message   .= "Model output is in: ".$job->{'directory'}."\n\n";
    if (
	exists($launchScript->{'config'}->{'contact'}->{'email'}) 
	&&
	$launchScript->{'config'}->{'contact'}->{'email'} =~ m/\@/ 
	&&
	$launchScript->{'emailReport'} eq "yes" 
	) {
	require MIME::Lite;
	$message .= "Log file is attached.\n";
	my $msg = MIME::Lite->new(
	    From    => 'Galacticus',
	    To      => $launchScript->{'config'}->{'contact'}->{'email'},
	    Subject => 'Galacticus model failed',
	    Type    => 'TEXT',
	    Data    => $message
	    );
	$msg->attach(
	    Type     => "text/plain",
	    Path     => $job->{'directory'}."/galacticus.log",
	    Filename => "galacticus.log"
	    );
	$msg->send;
    } else {
	print $message;
	print "Log follows:\n";
	print read_file($job->{'directory'}."/galacticus.log");
    }
}

sub Analyze {
    # Model finished successfully.
    my $job          = shift();
    my $launchScript = shift();
    # Check if we need to merge models.
    push(@{$launchScript->{'mergeGroups'}->{$job->{'mergeGroup'}}},$job->{'directory'});
    return
	if ( scalar(@{$launchScript->{'mergeGroups'}->{$job->{'mergeGroup'}}}) < $launchScript->{'splitModels'} );
    if ( $launchScript->{'splitModels'} > 1 ) {
	# We must merge the models before continuing.
	system(
	    $ENV{'GALACTICUS_EXEC_PATH'}."/scripts/aux/Merge_Models.pl ".
	    join(" ",map {$_."/galacticus.hdf5"} @{$launchScript->{'mergeGroups'}->{$job->{'mergeGroup'}}})
	    ." ".
	    ${$launchScript->{'mergeGroups'}->{$job->{'mergeGroup'}}}[0]."/galacticusMerged.hdf5"
	    );
    }
    # Perform analysis.
    if ( defined($job->{'analysis'}) ) {
	my $analysisScript = $job->{'directory'}."/analysis.sh";
	open(my $analysisFile,">".$analysisScript);
	print $analysisFile $job->{'analysis'};
	close($analysisFile);
	&System::Redirect::tofile(
	    "chmod u=wrx ".$job->{'directory'}."/analysis.sh;".
	    $job->{'directory'}."/analysis.sh",
	    $job->{'directory'}."/analysis.out"
	    );
    }
}

sub CleanUp {
    # Clean up after a job is finished.
    my $job          = shift();
    my $launchScript = shift();
    # Compress output if requested.
    &IO::Compress::Simple::Compress_Directory($job->{'directory'})
	if ( $launchScript->{'compressModels'} eq "yes" );
}

1;
