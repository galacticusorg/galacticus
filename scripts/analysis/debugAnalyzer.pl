#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

# Extract information from debug logs and export to CSV for further analysis.
# Andrew Benson (22-August-2023)

# Get the name of the debug log file.
die("Usage: debugAnalyzer.pl <debugFile>")
    unless ( scalar(@ARGV) == 1 );
my $debugFileName = $ARGV[0];

# Initialize step counter and data.
my $i = -1;
my $data;

# Open and parse the debug log file.
open(my $debugLog,$debugFileName);
while ( my $line = <$debugLog> ) {
    if ( $line =~ m/^\s*step:\s+([0-9\.eE\+\-]+)/ ) {
	++$i;
	$data->{'time'}->[$i] = $1;
    } elsif ( $line =~ m/^\s*value:\s+([a-zA-Z:]+)\s+([\d\.eE\+\-]+)/ ) {
	$data->{'properties'}->{$1}->{'value'}      ->[$i] = $2;
    } elsif ( $line =~ m/^\s*rate:\s+\(([a-zA-Z0-9_]+)\)\s+([a-zA-Z:\d\[\]]+)\s+([\d\.eE\+\-]+)/ ) {
	$data->{'properties'}->{$2}->{'rate' }->{$1}->[$i] = $3;
    } else {
	print $line;
	die("failed to parse line");
    }


    ## AJB HACK
    last
	if ( $i > 10 );

}
close($debugLog);
my $stepCount = $i+1;

# Output values and rates.
my @values = sort(keys(%{$data->{'properties'}}));
my @rates;
foreach my $property ( @values ) {
    push(@rates,map {{property => $property, function => $_}} sort(keys(%{$data->{'properties'}->{$property}->{'rate'}})));
}
open(my $valueLog,">","debugValues.csv");
print $valueLog "time , ".join(" , ",@values)."\n";
for(my $j=0;$j<$stepCount;++$j) {
    print $valueLog $data->{'time'}->[$j]." , ".join(" , ",map {defined($data->{'properties'}->{$_}->{'value'}->[$j]) ? $data->{'properties'}->{$_}->{'value'}->[$j] : "0.000000E+00"} @values)."\n";
}
close($valueLog);
open(my $ratesLog,">","debugRates.csv");
print $ratesLog "time , ".join(" , ",map {$_->{'property'}.":".$_->{'function'}} @rates)."\n";
for(my $j=0;$j<$stepCount;++$j) {
    print $ratesLog $data->{'time'}->[$j]." , ".join(" , ",map {my $rate = $data->{'properties'}->{$_->{'property'}}->{'rate'}->{$_->{'function'}}->[$j]; defined($rate) ? $rate : "0.000000E+00"} @rates)."\n";
}
close($ratesLog);
exit;
