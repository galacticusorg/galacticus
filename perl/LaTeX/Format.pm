package LaTeX::Format;
use POSIX;

# Formatting functions for LaTeX.
# Andrew Benson (28-July-2016)

sub Number {
    my $value   = shift();
    my $sigFigs = shift();
    my (%options) = @_
	if ( scalar(@_) > 0 );
    # Handle zero.
    my $valueFormatted;
    if ( $value == 0.0 ) {
	$valueFormatted = "0.0";
    } else {
	# Handle negative values.
	my $isNegative = $value < 0.0 ? "-" : "";
	$value = abs($value);
	# Determine order.
	my $o = floor(log($value)/log(10.0));
	my $order = floor(log10($value));
	# Scale value.
	$value /= 10.0**$order;
	# Format value.
	my $format = "%".$sigFigs.".".($sigFigs-1)."f";    
	# Created the formatted number.
	$valueFormatted = $isNegative.sprintf($format,$value)." \\times 10^{".$order."}";
    }
    # Add math mode delimiters?
    my $mathMode = exists($options{'mathMode'}) ? $options{'mathMode'} : 1;
    $valueFormatted = "\$".$valueFormatted."\$"
	if ( $mathMode );
    # Return the formatted number.
    return $valueFormatted;
}

1;
