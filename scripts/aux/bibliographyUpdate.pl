#!/usr/bin/env perl
use strict;
use warnings;
use WWW::Curl::Easy;
use JSON::PP qw(encode_json decode_json);
use Text::BibTeX;
use Data::Dumper;

# Update bibliography entries
# Andrew Benson (18-March-2024)
## Entries that are arXiv records will be pulled from NASA ADS which will return the actual journal reference where available.

# Read arguments.
die("Usage: bibliographyUpdate.pl <apiToken>")
    unless ( scalar(@ARGV) == 1 );
my $apiToken = $ARGV[0];

# Construct a curl object.
my $curl = WWW::Curl::Easy->new();
$curl->setopt(CURLOPT_HEADER,1);
$curl->setopt(CURLOPT_HTTPHEADER, ['Authorization: Bearer '.$apiToken]);

# First read in all bibliography entries in their raw form. This will allow us to re-output the entirely unchanged if there is no
# need to update.
my @entries;
open(my $bibliographyCurrent,$ENV{'GALACTICUS_EXEC_PATH'}."/doc/Galacticus.bib");
while ( my $line = <$bibliographyCurrent> ) {
    ++$#entries
	if ( $line =~ m/^\@/ );
    $entries[$#entries] .= $line;
}
close($bibliographyCurrent);

# Parse the bibliography and process records.
open(my $bibliographyNew,">",$ENV{'GALACTICUS_EXEC_PATH'}."/doc/Galacticus.bib.new");
my $bibliography = Text::BibTeX::File->new($ENV{'GALACTICUS_EXEC_PATH'}."/doc/Galacticus.bib");
my $i            = -1;
my $countUpdated =  0;
while ( my $entry = Text::BibTeX::Entry->new($bibliography) ) {
    ++$i;
    unless ( $entry->parse_ok() ) {
	print $entries[$i];
	die("failed to parse bibliography entry");
    }
    # By default, assume we can simply reuse the raw text of the existing entry.
    my $entryText = $entries[$i];
    # Check for an arXiv entry.
    if ( $entry->exists('journal') ) {
	my $journal = $entry->get('journal');
	if ( $journal =~ m/arxiv/i ) {
	    my $identifier;
	    my $url;
	    $url = $entry->get('url')
		if ( $entry->exists('url') );
	    $url = $entry->get('adsurl')
		if ( $entry->exists('adsurl') );
	    if ( defined($url) ) {
		if ( $url =~ m/^https??:\/\/(ui\.)??adsabs\.harvard\.edu\/abs\/(.+)/ ) {
		    # URL is an ADS arXiv reference. Extract the identifier directly.
		    $identifier = $2;
		} elsif ( $url =~ m/^https??:\/\/arxiv\.org\/abs\/([\d\.]+)/ ) {
		    # URL is a direct arXiv reference. Construct the corresponding ADS identifier.
		    $identifier = $1;
		    if ( $entry->exists('year') ) {
			my $year = $entry->get('year');
			$identifier = $year."arXiv".$identifier;
		    } else {
			print "no 'year' field exists\n";
		    }
		    if ( $entry->exists('author') ) {
			my $author = $entry->get('author');
			$identifier .= substr($author,0,1);
		    } else {
			print "no 'author' field exists\n";
		    }
		} else {
		    print "URL '".$url."' is not recognized\n";
		}
	    } else {
		die("no URL found");
	    }
	    die('unable to update entry for arXiv record - key is "'.$entry->key().'"')
		unless ( defined($identifier) );
	    # We have an identifier for an arXiv reference. Pull the ADS record for this identifier.
	    $curl->setopt(CURLOPT_URL, 'https://api.adsabs.harvard.edu/v1/export/bibtex');
	    $curl->setopt(CURLOPT_HTTPHEADER, ['Authorization: Bearer '.$apiToken,"Content-Type: application/json"]);
	    my $response_body;
	    $curl->setopt(CURLOPT_WRITEDATA,\$response_body);
	    $curl->setopt(CURLOPT_POST, 1);
	    $curl->setopt(CURLOPT_POSTFIELDS, '{"bibcode": ["'.$identifier.'"]}');
	    my $retcode = $curl->perform;
	    if ($retcode == 0) {
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
		my $bibtexEntry = decode_json($json);
		if (
		    exists ($bibtexEntry->{'export'})
		    &&
		    defined($bibtexEntry->{'export'})
		    ) {
		    # If the returned record is something other than an arXiv reference, update to this new record.
		    unless ( 
			$bibtexEntry->{'export'} =~ m/^\s*adsurl\s*=\s*\{https:\/\/ui\.adsabs\.harvard\.edu\/abs\/\d+arXiv[\d\.]+[A-Z]\}/m
			||
			$bibtexEntry->{'export'} =~ m/^\s*adsurl\s*=\s*\{https:\/\/ui\.adsabs\.harvard\.edu\/abs\/\d+astro\.ph\.\d+[A-Z]\}/m
			) {
			# Replace the BibTeX key of this record with the original key - we do not want to have to re-write references
			# in all of our documentation.
			++$countUpdated;
			my $key = $entry->key();
			($entryText = $bibtexEntry->{'export'}) =~ s/^\@ARTICLE\{.*,/\@article\{$key,/;
		    }
		}
	    } else {
		die("Failed to retrieve record identifier '".$identifier."': ".$retcode." ".$curl->strerror($retcode)." ".$curl->errbuf);
	    }
	}
    }
    # Output the (possibly updated) record.
    print $bibliographyNew $entryText;
}
close($bibliographyNew);

# Replace the original bibliography.
system("mv ".$ENV{'GALACTICUS_EXEC_PATH'}."/doc/Galacticus.bib.new ".$ENV{'GALACTICUS_EXEC_PATH'}."/doc/Galacticus.bib");

# Report.
print "Updated ".$countUpdated." of ".scalar(@entries)." records\n";

exit;
