#!/usr/bin/env perl
use strict;
use warnings;

# Check for broken links in Galacticus documentation.
# Andrew Benson (21-September-2020)

# Extract PDF destinations.
our $pdfDestinations;
foreach my $suffix ( "Usage", "Physics", "Development", "Source" ) {
    system("curl --silent -L --insecure --output Galacticus_".$suffix.".pdf https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_".$suffix.".pdf");
    system("./scripts/aux/pdfDestinationsExtract.py Galacticus_".$suffix.".pdf Galacticus_".$suffix.".dests");
    open(my $destFile,"./Galacticus_".$suffix.".dests");
    while ( my $dest = <$destFile> ) {
	chomp($dest);
	$pdfDestinations->{$suffix}->{$dest} = 1;
    }
    close($destFile);
}

# Initialize list of broken URLs.
our @brokenURLs;
our $urlCount = 0;

# Scan files in the doc folder.
my $docPath = "./doc";
opendir(my $docFolder,$docPath);
while ( my $fileName = readdir($docFolder) ) {
    &scanFile($fileName,$docPath)
	if ( $fileName =~ /^[A-Z].*\.tex$/ );
}
closedir($docFolder);

# Scan files in the source folder.
my $sourcePath = "./source";
opendir(my $sourceFolder,$sourcePath);
while ( my $fileName = readdir($sourceFolder) ) {
    &scanFile($fileName,$sourcePath)
	if ( $fileName =~ /\.(F90|Inc)$/ );
}
closedir($sourceFolder);

# Scan files in the wiki.
system("git clone https://github.com/galacticusorg/galacticus.wiki.git");
my $wikiPath = "./galacticus.wiki";
system("cd ".$wikiPath."; git pull");
opendir(my $wikiFolder,$wikiPath);
while ( my $fileName = readdir($wikiFolder) ) {
    &scanWiki($fileName,$wikiPath)
	if ( $fileName =~ /\.md$/ );
}
closedir($wikiFolder);

# Check status.
my $status = scalar(@brokenURLs) != 0;
exit $status;

sub scanFile {
    my $fileName   = shift();
    my $path       = shift();
    my $lineNumber = 0;
    open(my $file,$path."/".$fileName);
    while ( my $line = <$file> ) {
	++$lineNumber;
	while ( $line =~ m/\\href\{([^\}]+)\}/ ) {
	    my $url = $1;
	    $line =~ s/\\href\{([^\}]+)\}//;
            # Remove escapes.
	    $url =~ s/\\(.)/$1/g;
	    unless ( &checkLink($url) ) {
		print "Broken link: \"".$url."\" in ".$path."/".$fileName." line ".$lineNumber."\n";
		push(@brokenURLs,$url." in ".$path."/".$fileName." line ".$lineNumber);
	    }
	}
	while ( $line =~ m/\\refPhysics\{([^\}]+)\}/ ) {
	    my $ref = $1;
	    my $url = "https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#physics.".$ref;
	    $line =~ s/\\refPhysics\{([^\}]+)\}//;
	    unless ( &checkLink($url) ) {
		print "Broken refPhysics link: \"".$ref."\" in ".$path."/".$fileName." line ".$lineNumber."\n";
		push(@brokenURLs,$ref." (refPhysics) in ".$path."/".$fileName." line ".$lineNumber);
	    }
	}
	while ( $line =~ m/\\refClass\{([^\}]+)\}/ ) {
	    my $ref = $1;
	    my $url = "https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf\#class.".$ref;
	    $line =~ s/\\refClass\{([^\}]+)\}//;
	    unless ( &checkLink($url) ) {
		print "Broken refClass link: \"".$ref."\" in ".$path."/".$fileName." line ".$lineNumber."\n";
		push(@brokenURLs,$ref." (refClass) in ".$path."/".$fileName." line ".$lineNumber);
	    }
	}
    }
    close($file);
}

sub scanWiki {
    my $fileName   = shift();
    my $path       = shift();
    my $lineNumber = 0;
    open(my $file,$path."/".$fileName);
    while ( my $line = <$file> ) {
	++$lineNumber;
	while ( $line =~ m/\[[^\]]+\]\(([^\)]+)\)/ ) {
	    my $url = $1;
	    $line =~ s/\[[^\]]+\]\(([^\)]+)\)//;
	    unless ( &checkLink($url) ) {
		print "Broken link: \"".$url."\" in ".$path."/".$fileName." line ".$lineNumber." (see preceededing log)\n";
		push(@brokenURLs,$url." in ".$path."/".$fileName." line ".$lineNumber);
	    }
	}
    }
    close($file);
}

sub checkLink {
    my $url = shift();
    my $status;
    $url =~ s/\\#/#/g;
    if ( $url =~ m/^mailto:/ ) {
	# Ignore mailto URLs.
	$status = 1;
    } elsif ( $url =~ m/^#/ ) {
	# Ignore in-page anchors.
	$status = 1;
    } elsif ( $url =~ m/^https:\/\/github\.com\/galacticusorg\/galacticus\/releases\/download\/bleeding-edge\/Galacticus_(Usage|Physics|Development|Source)\.pdf\\??#(.+)/ ) {
	# Link to an anchor in Galacticus PDF documentation.
	my $suffix = $1;
	my $anchor = $2;
	$status = exists($pdfDestinations->{$suffix}) && exists($pdfDestinations->{$suffix}->{$anchor});
    } else {
	# An external link. Include a short sleep here to rate limit requests.
	## --cipher 'DEFAULT:!DH' - this reduces the default security level which otherwise prevents some URLs from being downloaded.
	## --range 0-0 - this causes no bytes to actually be downloaded - this is disabled on some sites as it seems to break them
	my $options = "--silent --insecure --location --output /dev/null --fail --cipher 'DEFAULT:!DH'";
	$options .= " --range 0-0"
	    unless ( $url =~ m/^https:\/\/www\.drdobbs\.com\// );
	$options .= " --user-agent \"Mozilla\""
	    if ( $url =~ m/sharepoint\.com/ );
	$options .= " --compressed"
	    if ( $url =~ m/docs\.github\.com/ );
	system("sleep 1; curl ".$options." \"".$url."\"");
	$status = $? == 0 ? 1 : 0;	
	$options =~ s/\-\-silent//;
	$options =~ s/\-\-fail//;
	system("sleep 1; curl ".$options." \"".$url."\"")
	    unless ( $status );
    }
    return $status;
}
