#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use XML::Simple;
use POSIX qw(log10);
use Galacticus::Options;
use Scalar::Util qw(looks_like_number);
use Data::Dumper;

# Parse an MCMC config file and generate a LaTeX table describing the priors applied to each parameter.
# Andrew Benson (15-April-2024)

# Read command line arguments and options.
die("Usage: createPriorsTable.pl <configFileName> [options...]")
    unless ( scalar(@ARGV) > 1 );
my $configFileName = $ARGV[0];
my %options =
    (
     outputFile => "priors.tex",
     format     => "%+.2f"
     );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse parameter group labels.
my %groupLabels;
if ( exists($options{'groupBefore'}) ) {
    foreach my $option ( &List::ExtraUtils::as_array($options{'groupBefore'}) ) {
	if ( $option =~ m/^([a-zA-Z0-9:]+)\s+(.*)/ ) {
	    my $name = $1;
	    my $label = $2;
	    $groupLabels{$name} = $label;
	} else {
	    die('invalid parameter name');
	}
    }
}

# Read the config file
my $xml = new XML::Simple();
my $config = $xml->XMLin($configFileName);

# Build the table.
open(my $output,">",$options{'outputFile'});
print $output "\\begin{tabular}{ll}\n";
print $output "\\hline\n";
print $output "\\textbf{Parameter} & \\textbf{Prior} \\\\\n";
print $output "\\hline\n";
my $i = -1;
foreach my $modelParameter ( @{$config->{'posteriorSampleSimulation'}->{'modelParameter'}} ) {
    ++$i;
    my $name = exists($modelParameter->{'label'}) ? "\$".$modelParameter->{'label'}->{'value'}."\$" : $modelParameter->{'name'}->{'value'};
    my $priorName;
    my @priorParameters;
    if ( $modelParameter->{'distributionFunction1DPrior'}->{'value'} eq "uniform" ) {
	$priorName          = "Uniform";
	$priorParameters[0] = $modelParameter->{'distributionFunction1DPrior'}->{'limitLower'}->{'value'};
	$priorParameters[1] = $modelParameter->{'distributionFunction1DPrior'}->{'limitUpper'}->{'value'};
    } elsif ( $modelParameter->{'distributionFunction1DPrior'}->{'value'} eq "logUniform" ) {
	if ( $name =~ m/^\$(.+)\$$/ ) {
	    $name = "\$\\log_{10}(".$1.")\$";
	} else {
	    $name = "\$\\log_{10}(\\hbox{".$name."})\$";
	}
	$priorName          = "Uniform";
	$priorParameters[0] = log10($modelParameter->{'distributionFunction1DPrior'}->{'limitLower'}->{'value'});
	$priorParameters[1] = log10($modelParameter->{'distributionFunction1DPrior'}->{'limitUpper'}->{'value'});
    } elsif ( $modelParameter->{'distributionFunction1DPrior'}->{'value'} eq "normal" ) {
	$priorName          = "Normal";
	$priorParameters[0] = $modelParameter->{'distributionFunction1DPrior'}->{'mean'}->{'value'};
	$priorParameters[1] = $modelParameter->{'distributionFunction1DPrior'}->{'variance'}->{'value'};
	if ( exists($modelParameter->{'distributionFunction1DPrior'}->{'limitLower'}) ) {
	    push(@priorParameters,$modelParameter->{'distributionFunction1DPrior'}->{'limitLower'}->{'value'});
	} elsif ( exists($modelParameter->{'distributionFunction1DPrior'}->{'limitUpper'}) ) {
	    push(@priorParameters,"\$-\\infty\$");
	}
	if ( exists($modelParameter->{'distributionFunction1DPrior'}->{'limitUpper'}) ) {
	    push(@priorParameters,$modelParameter->{'distributionFunction1DPrior'}->{'limitUpper'}->{'value'});
	} elsif ( exists($modelParameter->{'distributionFunction1DPrior'}->{'limitLower'}) ) {
	    push(@priorParameters,"\$+\\infty\$");
	}
    } elsif ( $modelParameter->{'distributionFunction1DPrior'}->{'value'} eq "logNormal" ) {
	$priorName          = "Lognormal";
	$priorParameters[0] = $modelParameter->{'distributionFunction1DPrior'}->{'mean'}->{'value'};
	$priorParameters[1] = $modelParameter->{'distributionFunction1DPrior'}->{'variance'}->{'value'}
    } else {
	die("unknown prior");
    }
    # Write any parameter group label.
    if ( exists($groupLabels{$modelParameter->{'name'}->{'value'}}) ) {
	print $output "\\hline\n"
	    unless ( $i == 0 );
	print $output "\\multicolumn{2}{c}{\\emph{".$groupLabels{$modelParameter->{'name'}->{'value'}}."}} \\\\\n";
    }
    # Write the parameter prior.
    print $output $name." & \\hbox{".$priorName."}(".join(",",map {looks_like_number($_) ? sprintf($options{'format'},$_) : $_} @priorParameters).") \\\\\n";
}
print $output "\\hline\n";
print $output "\\end{tabular}\n";
close($output);

exit;
