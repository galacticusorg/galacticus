# Contains a Perl module which implements useful compression tools.

package Simple;
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);

my $status = 1;
$status;

sub Compress_Directory {
    # Compress all files in a directory.

    # Get the directory name.
    $directoryName = shift;

    opendir(dirHndl,$directoryName);
    while ( $fileName = readdir(dirHndl) ) {
	unless ( $fileName =~ m/^\./ || $fileName =~ m/\.gz$/ || $fileName =~ m/\.bz2/ ) {
	    if ( bzip2 "<".$directoryName."/".$fileName.">" => "<".$directoryName."/".$fileName.".bz2>" ) {
		unlink($directoryName."/".$fileName);
	    } else {
		print "Failed to compress ".$directoryName."/".$fileName."\n";
	    }
	}
    }
    closedir(dirHndl);
}
