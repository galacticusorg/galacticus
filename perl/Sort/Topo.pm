package Sort::Topo;
use strict;
use warnings;
use List::MoreUtils 'first_index';
use Data::Dumper;

sub sort {
    # Perform a topological sort. Based on the example from \href{https://rosettacode.org/wiki/Topological_sort}.
    ## "objects" is a list of the objects to be sorted.
    my @objects      = @{shift()};
    ## "dependencies" is a hash in which each key is an object name, and the value is an array of objects names upon which it depends.
    my %dependencies = %{shift()};
    # First build an array of dependencies.
    my @dependency;
    foreach my $object ( keys(%dependencies) ) {
	# Skip cases where the object is not present in the object list.
	next
	    unless ( grep {$_ eq $object} @objects );
	my $objectIndex = first_index {$_ eq $object} @objects;
	foreach my $dependent ( @{$dependencies{$object}} ) {
	    next
		unless ( grep {$_ eq $dependent} @objects );
	    my $dependentIndex = first_index {$_ eq $dependent} @objects;
	    push(@dependency,[$dependentIndex,$objectIndex]);
	}
    }
    # Count objects and dependencies, and initialize the order and position arrays.
    my $countDependencies = scalar(@dependency);
    my $countObjects      = scalar(@objects   );
    my @order             = 0..$countObjects-1;
    my @position          = 0..$countObjects-1;
    # Begin the sort.
    my $k = 0;
    my $j    ;
    do {
	$j = $k           ;
	$k = $countObjects;
	for(my $i=0;$i<$countDependencies;++$i) {
            my $dependencyLeft  = $dependency[$i][0];
            my $dependencyRight = $dependency[$i][1];
            my $positionLeft    = $position[$dependencyLeft ];
            my $positionRight   = $position[$dependencyRight];
	    next
		if ( $dependencyLeft == $dependencyRight || $positionLeft >= $k || $positionLeft < $j || $positionRight < $j);
	    --$k;
            $position[$order         [$k]]=  $positionLeft    ;
            $position[$dependencyLeft    ]=                $k ;
            $order   [  $positionLeft    ]=$order         [$k];
            $order   [                $k ]=$dependencyLeft    ;
	}
    } until ( $k <= $j );
    # Check for a circular dependency.
    my $countOrdered = $j;
    die("circular dependency")
	if ($countOrdered < $countObjects);
    # Construct the final ordered list of objects.
    return map {$objects[$_]} @order;
}

1;
