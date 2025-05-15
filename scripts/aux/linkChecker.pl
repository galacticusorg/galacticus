#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use System::Redirect;
use Regexp::Common;
use WWW::Curl::Easy;
use XML::Simple;
use JSON::PP qw(encode_json decode_json);

# Check for broken links in Galacticus documentation.
# Andrew Benson (21-September-2020)

# Read arguments.
die("Usage: linkChecker.pl <apiToken>")
    unless ( scalar(@ARGV) == 1 );
my $apiToken = $ARGV[0];

# Extract PDF destinations.
my $pdfDestinations;
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

# Initialize records of number of consecutive checks for which a URL has failed.
my $xml = new XML::Simple();
my $failures =  -e "linkCheckFailures.xml" ? $xml->XMLin("linkCheckFailures.xml") : {};
print "Current consecutive failure count (in order of increasing number of failures:\n";
foreach my $url ( sort { $failures->{'url'}->{$a}->{'consecutiveFailures'} <=> $failures->{'url'}->{$b}->{'consecutiveFailures'}} keys(%{$failures->{'url'}}) ) {
    print $failures->{'url'}->{$url}->{'consecutiveFailures'}.":\t".$url."\n";
}

# Initialize structure of links to check.
my $urls;

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

# Check the URLs.
my $status = &checkURLs($urls,$pdfDestinations,$apiToken,\$failures);

# Store records of consecutive failures.
open(my $record,">","linkCheckFailures.xml");
print $record $xml->XMLout($failures, RootName => "failures");
close($record);

# Finished.
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
	    push(@{$urls->{$url}},{file => $fileName, path => $path, lineNumber => $lineNumber});
	}
	while ( $line =~ m/\\refPhysics\{([^\}]+)\}/ ) {
	    my $ref = $1;
	    my $url = "https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#physics.".$ref;
	    $line =~ s/\\refPhysics\{([^\}]+)\}//;
	    push(@{$urls->{$url}},{type => "refPhysics", ref => $ref, file => $fileName, path => $path, lineNumber => $lineNumber});
	}
	while ( $line =~ m/\\refClass\{([^\}]+)\}/ ) {
	    my $ref = $1;
	    my $url = "https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf\#class.".$ref;
	    $line =~ s/\\refClass\{([^\}]+)\}//;
	    push(@{$urls->{$url}},{type => "refClass", ref => $ref, file => $fileName, path => $path, lineNumber => $lineNumber});
	}
    }
    close($file);
}

sub scanWiki {
    my $fileName   = shift();
    my $path       = shift();
    my $lineNumber = 0;
    my $bracket    = $RE{balanced}{-parens=>'[]'};
    my $parens     = $RE{balanced}{-parens=>'()'};
    open(my $file,$path."/".$fileName);
    while ( my $line = <$file> ) {
	++$lineNumber;
	while ( $line =~ m/$bracket$parens/ ) {
	    (my $url = $2) =~ s/^\((.*)\)$/$1/;
	    $line =~ s/\[[^\]]+\]\(([^\)]+)\)//;
	    push(@{$urls->{$url}},{file => $fileName, path => $path, lineNumber => $lineNumber});
	}
    }
    close($file);
}

sub checkURLs {
    # Check URLs are valid.
    my $urls            =   shift() ;
    my $pdfDestinations =   shift() ;
    my $apiToken        =   shift() ;
    my $failures        = ${shift()};
    my $status          = 0;
    my $bibCodes;
    foreach my $urlKey ( keys(%{$urls}) ) {
	(my $url = $urlKey) =~ s/\\#/#/g;
	# Ignore mailto URLs.
	next
	    if ( $url =~ m/^mailto:/ );
	# Ignore in-page anchors.
	next
	    if ( $url =~ m/^#/ );
	# Process links.
	if ( $url =~ m/^https:\/\/github\.com\/galacticusorg\/galacticus\/releases\/download\/bleeding-edge\/Galacticus_(Usage|Physics|Development|Source)\.pdf\\??#(.+)/ ) {
	    # Link to an anchor in Galacticus PDF documentation.
	    my $suffix = $1;
	    my $anchor = $2;
	    unless ( exists($pdfDestinations->{$suffix}) && exists($pdfDestinations->{$suffix}->{$anchor}) ) {
		$status = 1;
		print "Broken ".$urls->{$urlKey}->[0]->{'type'}."{".$urls->{$urlKey}->[0]->{'ref'}."} link in:\n";
		foreach my $source ( @{$urls->{$urlKey}} ) {
		    print "\t ".$source->{'path'}."/".$source->{'file'}." line ".$source->{'lineNumber'}."\n";
		}
	    }   
	} elsif ( $url =~ m/adsabs\.harvard\.edu\/abs\/([^\/]+)/ ) {
	    # NASA ADS link - accumulate the bibCode.
	    my $bibCode = $1;
	    # Translate escape codes.
	    $bibCode =~ s/%26/&/;
	    push(@{$bibCodes->{$bibCode}->{'sources'}},@{$urls->{$urlKey}});
	    ${$bibCodes->{$bibCode}->{'urls'}}{$urlKey} = 1;
	} else {
	    # An external link. Include a short sleep here to rate limit requests.
	    sleep(1);
	    ## --cipher 'DEFAULT:!DH' - this reduces the default security level which otherwise prevents some URLs from being downloaded.
	    ## --range 0-0 - this causes no bytes to actually be downloaded - this is disabled on some sites as it seems to break them.
	    my $options = "--max-time 60 --insecure --location --output /dev/null --fail-with-body --cipher 'DEFAULT:!DH'";
	    $options .= " --range 0-0"
		unless ( $url =~ m/^https:\/\/www\.drdobbs\.com\// );
	    $options .= " --user-agent \"Mozilla\""
		if ( $url =~ m/sharepoint\.com/ );
	    # docker.com issues a 403 unless we make cURL pretend to be wget...
	    $options .= " --user-agent \"Wget/1.21.2\""
		if ( $url =~ m/docker\.com/ );
	    $options .= " --compressed"
		if ( $url =~ m/docs\.github\.com/ );
	    $options .= " --http1.1"
		if ( $url =~ m/camb\.info/ );
	    $options .= " --retry 5"
		if ( $url =~ m/www\.gnu\.org/ );
	    my $sleepTime = 1;
	    &System::Redirect::tofile("curl ".$options." \"".$url."\"","curl.log");
	    my $error = $?;
	    if ( $error ) {
		# Check for known problems.
		open(my $logFile,"curl.log");
		while ( my $line = <$logFile> ) {
		    # Some servers do not correctly terminate their connections. Ignore such cases.
		    if ( $url =~ m/http:\/\/heasarc\.gsfc\.nasa\.gov\/xanadu\/xspec\// ) {
			if ( $line =~ m/error:0A000126:SSL routines::unexpected eof while reading, errno 0/ ) {
			    $error = 0;
			    last;
			}
		    }
		}
		close($logFile);
	    }
	    if ( $error ) {		
		$status = 1
		    if ( &recordFailure($url,$failures) );
		print "Broken link: \"".$url."\" (for past ".$failures->{'url'}->{$url}->{'consecutiveFailures'}." attempts) in:\n";
		foreach my $source ( @{$urls->{$urlKey}} ) {
		    print "\t".$source->{'path'}."/".$source->{'file'}." line ".$source->{'lineNumber'}."\n";
		}
		print "Log:\n";
		open(my $logFile,"curl.log");
		while ( my $line = <$logFile> ) {
		    print $line;
		}
		close($logFile);
	    } else {
		&recordSuccess($url,$failures);
	    }
	}
    }
    # Check NASA ADS bibcodes.
    my $curl = WWW::Curl::Easy->new();
    $curl->setopt(CURLOPT_HEADER,1);
    my $countRecords = scalar(keys(%{$bibCodes}));
    $curl->setopt(CURLOPT_URL, 'https://api.adsabs.harvard.edu/v1/search/bigquery?q=*:*&rows='.$countRecords.'&fl=bibcode,alternate_bibcode,title,author,year,pub,volume,page');
    $curl->setopt(CURLOPT_HTTPHEADER, ['Authorization: Bearer '.$apiToken,"Content-Type: big-query/csv"]);
    my $response_body;
    $curl->setopt(CURLOPT_WRITEDATA,\$response_body);
    $curl->setopt(CURLOPT_POST, 1);
    $curl->setopt(CURLOPT_POSTFIELDS, "bibcode\n".join("\n",sort(keys(%{$bibCodes}))));
    my $retcode = $curl->perform;
    my $records;
    if ($retcode == 0) {
	my $response_code = $curl->getinfo(CURLINFO_HTTP_CODE);
	if ( $response_code == 200 ) {
	    # Extract the JSON.
	    my $json;
	    my $startFound = 0;
	    open(my $response,"<",\$response_body);
	    while ( my $line = <$response> ) {
		$startFound = 1
		    if ( $line =~ m/^\{/ );
		$json .= $line
		    if ( $startFound );
	    }
	    close($response);
	    $records = decode_json($json);
	} else {
	    $status = 1;
	    print "Failed to retrieve record identifiers: ".$response_code.$response_body."\n";
	}
    } else {
	$status = 1;
	print "Failed to retrieve record identifiers: ".$retcode." ".$curl->strerror($retcode)." ".$curl->errbuf."\n";
    }
    # Parse records.
    foreach my $entry ( @{$records->{'response'}->{'docs'}} ) {
	my $found = 0;
	if ( exists($bibCodes->{$entry->{'bibcode'}}) ) {
	    $bibCodes->{$entry->{'bibcode'}}->{'found'} = 1;
	    $found = 1;
	}
	if ( exists($entry->{'alternate_bibcode'}) ) {
	    foreach my $altBibCode ( @{$entry->{'alternate_bibcode'}} ) {
		if ( exists($bibCodes->{$altBibCode}) ) {
		    $bibCodes->{$altBibCode}->{'found'} = 1;
		    $found = 1;
		}
	    }
	}
	unless ( $found ) {
	    $status = 1;
	    print "Received unrequested record for bibcode '".$entry->{'bibcode'}."'\n";
	}
    }
    # Look for any bibcodes that were not found.
    foreach my $bibCode ( keys(%{$bibCodes}) ) {
	if ( exists($bibCodes->{$bibCode}->{'found'}) ) {
	    &recordSuccess($bibCode,$failures);
	} else {
	    $status = 1
		if ( &recordFailure($bibCode,$failures) );
	    print "Broken link (for past ".$failures->{'url'}->{$bibCode}->{'consecutiveFailures'}." attempts): {bibCode: ".$bibCode."} \"".join("; ",keys(%{$bibCodes->{$bibCode}->{'urls'}}))."\" in:\n";
	    foreach my $source ( @{$bibCodes->{$bibCode}->{'sources'}} ) {
		print "\t".$source->{'path'}."/".$source->{'file'}." line ".$source->{'lineNumber'}."\n";
	    }
	}
    }
    # Return final status.
    return $status;
}

sub recordFailure {
    # Record a failure to retrieve a URL.
    my $url      = shift();
    my $failures = shift();
    # Increment the number of consecutive fails by 1.
    $failures->{'url'}->{$url}->{'consecutiveFailures'} += 1;
    # Return a status depending on the number of consecutive failures.
    my $status = $failures->{'url'}->{$url}->{'consecutiveFailures'} >= 3;
    return $status;
}

sub recordSuccess {
    # Record a failure to retrieve a URL.
    my $url      = shift();
    my $failures = shift();
    # Zero the number of consecutive failures for this URL.
    $failures->{'url'}->{$url}->{'consecutiveFailures'} = 0;    
}
