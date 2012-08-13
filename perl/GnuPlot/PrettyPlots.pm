# Contains a Perl module which implements simple, yet pretty plotting of PDL datasets using GnuPlot.

package PrettyPlots;
use strict;
use warnings;
use PDL;
use Switch;

# A set of color definitions in #RRGGBB format.
our %colors = (
    Snow                 => "#FFFAFA",
    Snow2                => "#EEE9E9",
    Snow3                => "#CDC9C9",
    Snow4                => "#8B8989",
    GhostWhite           => "#F8F8FF",
    WhiteSmoke           => "#F5F5F5",
    Gainsboro            => "#DCCDC",
    FloralWhite          => "#FFFAF0",
    OldLace              => "#FDF5E6",
    Linen                => "#FAF0E6",
    AntiqueWhite         => "#FAEBD7",
    AntiqueWhite2        => "#EEDFCC",
    AntiqueWhite3        => "#CDC0B0",
    AntiqueWhite4        => "#8B8378",
    PapayaWhip           => "#FFEFD5",
    BlanchedAlmond       => "#FFEBCD",
    Bisque               => "#FFE4C4",
    Bisque2              => "#EED5B7",
    Bisque3              => "#CDB79E",
    Bisque4              => "#8B7D6B",
    PeachPuff            => "#FFDAB9",
    PeachPuff2           => "#EECBAD",
    PeachPuff3           => "#CDAF95",
    PeachPuff4           => "#8B7765",
    NavajoWhite          => "#FFDEAD",
    Moccasin             => "#FFE4B5",
    Cornsilk             => "#FFF8DC",
    Cornsilk2            => "#EEE8DC",
    Cornsilk3            => "#CDC8B1",
    Cornsilk4            => "#8B8878",
    Ivory                => "#FFFFF0",
    Ivory2               => "#EEEEE0",
    Ivory3               => "#CDCDC1",
    Ivory4               => "#8B8B83",
    LemonChiffon         => "#FFFACD",
    Seashell             => "#FFF5EE",
    Seashell2            => "#EEE5DE",
    Seashell3            => "#CDC5BF",
    Seashell4            => "#8B8682",
    Honeydew             => "#F0FFF0",
    Honeydew2            => "#E0EEE0",
    Honeydew3            => "#C1CDC1",
    Honeydew4            => "#838B83",
    MintCream            => "#F5FFFA",
    Azure                => "#F0FFFF",
    AliceBlue            => "#F0F8FF",
    Lavender             => "#E6E6FA",
    LavenderBlush        => "#FFF0F5",
    MistyRose            => "#FFE4E1",
    White                => "#FFFFFF",
    # Grays
    Black                => "#000000",
    DarkSlateGray        => "#2F4F4F",
    DimGray              => "#696969",
    SlateGray            => "#708090",
    LightSlateGray       => "#778899",
    Gray                 => "#BEBEBE",
    LightGray            => "#D3D3D3",
    # Blues
    MidnightBlue         => "#191970",
    Navy                 => "#000080",
    CornflowerBlue       => "#6495ED",
    DarkSlateBlue        => "#483D8B",
    SlateBlue            => "#6A5ACD",
    MediumSlateBlue      => "#7B68EE",
    LightSlateBlue       => "#8470FF",
    MediumBlue           => "#0000CD",
    RoyalBlue            => "#4169E1",
    Blue                 => "#0000FF",
    DodgerBlue           => "#1E90FF",
    DeepSkyBlue          => "#00BFFF",
    SkyBlue              => "#87CEEB",
    LightSkyBlue         => "#87CEFA",
    SteelBlue            => "#4682B4",
    LightSteelBlue       => "#B0C4DE",
    LightBlue            => "#ADD8E6",
    PowderBlue           => "#B0E0E6",
    PaleTurquoise        => "#AFEEEE",
    DarkTurquoise        => "#00CED1",
    MediumTurquoise      => "#48D1CC",
    Turquoise            => "#40E0D0",
    Cyan                 => "#00FFFF",
    LightCyan            => "#E0FFFF",
    CadetBlue            => "#5F9EA0",
    # Greens
    MediumAquamarine     => "#66CDAA",
    Aquamarine           => "#7FFFD4",
    DarkGreen            => "#006400",
    DarkOliveGreen       => "#556B2F",
    DarkSeaGreen         => "#8FBC8F",
    SeaGreen             => "#2E8B57",
    MediumSeaGreen       => "#3CB371",
    LightSeaGreen        => "#20B2AA",
    PaleGreen            => "#98FB98",
    SpringGreen          => "#00FF7F",
    LawnGreen            => "#7CFC00",
    Chartreuse           => "#7FFF00",
    MediumSpringGreen    => "#00FA9A",
    GreenYellow          => "#ADFF2F",
    LimeGreen            => "#32CD32",
    YellowGreen          => "#9ACD32",
    ForestGreen          => "#228B22",
    OliveDrab            => "#6B8E23",
    DarkKhaki            => "#BDB76B",
    Khaki                => "#F0E68C",
    Yellow               => "#FFFF00",
    PaleGoldenrod        => "#EEE8AA",
    LightGoldenrodYellow => "#FAFAD2",
    LightYellow          => "#FFFFE0",
    Gold                 => "#FFD700",
    LightGoldenrod       => "#EEDD82",
    Goldenrod            => "#DAA520",
    DarkGoldenrod        => "#B8860B",
    # Browns
    RosyBrown            => "#BC8F8F",
    IndianRed            => "#CD5C5C",
    SaddleBrown          => "#8B4513",
    Sienna               => "#A0522D",
    Peru                 => "#CD853F",
    Burlywood            => "#DEB887",
    Beige                => "#F5F5DC",
    Wheat                => "#F5DEB3",
    SandyBrown           => "#F4A460",
    Tan                  => "#D2B48C",
    Chocolate            => "#D2691E",
    Firebrick            => "#B22222",
    Brown                => "#A52A2A",
    # Oranges
    DarkSalmon           => "#E9967A",
    Salmon               => "#FA8072",
    LightSalmon          => "#FFA07A",
    Orange               => "#FFA500",
    DarkOrange           => "#FF8C00",
    Coral                => "#FF7F50",
    LightCoral           => "#F08080",
    Tomato               => "#FF6347",
    OrangeRed            => "#FF4500",
    Red                  => "#FF0000",
    # Pinks/Violets
    HotPink              => "#FF69B4",
    DeepPink             => "#FF1493",
    Pink                 => "#FFC0CB",
    LightPink            => "#FFB6C1",
    PaleVioletRed        => "#DB7093",
    Maroon               => "#B03060",
    MediumVioletRed      => "#C71585",
    VioletRed            => "#D02090",
    Violet               => "#EE82EE",
    Plum                 => "#DDA0DD",
    Orchid               => "#DA70D6",
    MediumOrchid         => "#BA55D3",
    DarkOrchid           => "#9932CC",
    DarkViolet           => "#9400D3",
    BlueViolet           => "#8A2BE2",
    Purple               => "#A020F0",
    MediumPurple         => "#9370DB",
    Thistle              => "#D8BFD8"
    );

# Sets of color pairs suitable for plotting points with a light middle and darker border.
our %colorPairs = (
    redYellow      => [$colors{'Red'          },$colors{'Yellow'        }],
    blackGray      => [$colors{'SlateGray'    },$colors{'Black'         }],
    blueCyan       => [$colors{'Blue'         },$colors{'Cyan'          }],
    peachPuff      => [$colors{'Bisque3'      },$colors{'PeachPuff'     }],
    slateGray      => [$colors{'DarkSlateGray'},$colors{'SlateGray'     }],
    cornflowerBlue => [$colors{'DarkSlateBlue'},$colors{'CornflowerBlue'}],
    lightSkyBlue   => [$colors{'DodgerBlue'   },$colors{'LightSkyBlue'  }],
    mediumSeaGreen => [$colors{'SeaGreen'     },$colors{'MediumSeaGreen'}],
    yellowGreen    => [$colors{'OliveDrab'    },$colors{'YellowGreen'   }],
    lightGoldenrod => [$colors{'Goldenrod'    },$colors{'LightGoldenrod'}],
    indianRed      => [$colors{'Sienna'       },$colors{'IndianRed'     }],
    orange         => [$colors{'OrangeRed'    },$colors{'Orange'        }],
    plum           => [$colors{'VioletRed'    },$colors{'Plum'          }],
    thistle        => [$colors{'MediumPurple' },$colors{'Thistle'       }]
    );

# Sets of sequences of color pairs suitable for plotting multiple datasets.
our %colorPairSequences = (
    sequence1 => [
	 "peachPuff"  ,"slateGray"     ,"cornflowerBlue","lightSkyBlue","mediumSeaGreen"
	,"yellowGreen","lightGoldenrod","indianRed"     ,"orange"      ,"plum"
        ,"thistle"
    ],
    slideSequence => [
	 "yellowGreen", "thistle", "orange", "lightGoldenrod"
    ]
    );

# Subroutine to generate plotting commands for a specified dataset and accumulate them to a buffer for later writing to GnuPlot.
sub Prepare_Dataset {
    # Extract the plot structure and datasets to plot.
    my $plot = shift;
    my $x    = shift;
    my $y    = shift;
    # Extract any remaining options.
    my (%options) = @_ if ( $#_ >= 1 );
    
    # Determine the GnuPlot version.
    my ($versionMajor,$versionMinor,$versionPatchLevel);
    my $gnuplotVersion = `gnuplot -V`;
    chomp($gnuplotVersion);
    if ( $gnuplotVersion =~ m/gnuplot (\d+)\.(\d+) patchlevel (\d+)/ ) {
	$versionMajor      = $1;
	$versionMinor      = $2;
	$versionPatchLevel = $3;
    } else {
	$versionMajor      = -1;
	$versionMinor      = -1;
	$versionPatchLevel = -1;
    }

    # Determine the plot style, assuming points by default.
    my $style = "point";
    $style = $options{'style'} if ( exists($options{'style'}) );

    # Create attributes for line weight, assuming no specification if weight option is not present.
    my %lineWeight;
    $lineWeight{'lower'} = "";
    $lineWeight{'upper'} = "";
    $lineWeight{'lower'} = " lw ".${$options{'weight'}}[0] if ( exists($options{'weight'}) );	
    $lineWeight{'upper'} = " lw ".${$options{'weight'}}[1] if ( exists($options{'weight'}) );	

    # Create attribute for line color, assuming no specification if color option is not present.
    my %lineColor;
    $lineColor{'lower'} = "";
    $lineColor{'upper'} = "";
    $lineColor{'lower'} = " lc rgbcolor \"".${$options{'color'}}[0]."\"" if ( exists($options{'color'}) );
    $lineColor{'upper'} = " lc rgbcolor \"".${$options{'color'}}[1]."\"" if ( exists($options{'color'}) );

    # Create attribute for line type, assuming no specification if type option is not present.
    my %lineType;
    $lineType{'lower'} = "";
    $lineType{'upper'} = "";
    $lineType{'lower'} = " lt ".$options{'linePattern'} if ( exists($options{'linePattern'}) );
    $lineType{'upper'} = " lt ".$options{'linePattern'} if ( exists($options{'linePattern'}) );

    # Create attribute for point type, assuming no specification if symbol option is not present.
    my %pointType;
    $pointType{'lower'} = "";
    $pointType{'upper'} = "";
    $pointType{'lower'} = " pt \"".${$options{'symbol'}}[0]."\"" if ( exists($options{'symbol'}) );
    $pointType{'upper'} = " pt \"".${$options{'symbol'}}[1]."\"" if ( exists($options{'symbol'}) );

    # Create attribute for point size, assuming no specification if pointSize option is not present.
    my %pointSize;
    $pointSize{'lower'} = "";
    $pointSize{'upper'} = "";
    $pointSize{'lower'} = " ps ".$options{'pointSize'} if ( exists($options{'pointSize'}) );
    $pointSize{'upper'} = " ps ".$options{'pointSize'} if ( exists($options{'pointSize'}) );

    # Create a title attribute, defaulting to no title if none is specified.
    my $title = " notitle";
    $title = " title \"".$options{'title'}."\"" if ( exists($options{'title'}) );

    # Define the dummy point and end points.
    my $dummyPoint;
    if (
	$versionMajor > 4
	||
	( $versionMajor == 4 &&
	  ( $versionMinor > 4 
	    ||
	    ( $versionMinor == 4 && $versionPatchLevel >= 2 )
	  )
	)
       )
    {
	$dummyPoint = "0 0\n";
    } else {
	$dummyPoint = "inf inf\n";
    }
    my $endPoint   = "e\n";

    # We have three plotting phases:
    #  keyLower - the lower layer of the key.
    #  keyUpper - the upper layer of the key.
    #  data     - actually plots the data.
    my %phaseRules = (
	keyLower => {
	    level => "lower"
	},
	keyUpper => {
	    level => "upper"
	},
	data => {
	}
	);

    # Loop over phases.
    foreach my $phase ( keys(%phaseRules) ) {
	# Initialize phase prefix if necessary.
	${$plot}->{$phase}->{'prefix'} = "plot" unless ( exists(${$plot}->{$phase}->{'prefix'}) );
	# Branch depending on style of dataset, points, boxes or lines.
	switch ( $style ) {
	    case ("line") {
		# Check if we are asked to plot just a single level.
		if ( exists($phaseRules{$phase}->{'level'}) ) {
		    # Plot just a single level, no real data.
		    ${$plot}->{$phase}->{'command'} .= ${$plot}->{$phase}->{'prefix'}." '-' with lines".$title
			.$lineType  {$phaseRules{$phase}->{'level'}}
		        .$lineColor {$phaseRules{$phase}->{'level'}}
		        .$lineWeight{$phaseRules{$phase}->{'level'}};
		    ${$plot}->{$phase}->{'data'   } .= $dummyPoint;
		    ${$plot}->{$phase}->{'data'   } .= $endPoint;
		    ${$plot}->{$phase}->{'prefix'} = ",";
		} else {
		    # Plot the actual data.
		    foreach my $level ( 'lower', 'upper' ) {
			${$plot}->{$phase}->{'data'} .= "plot '-' with lines notitle".$lineType{$level}.$lineColor{$level}.$lineWeight{$level}."\n";
			for(my $iPoint=0;$iPoint<nelem($x);++$iPoint) {
			    ${$plot}->{$phase}->{'data'} .= $x->index($iPoint)." ".$y->index($iPoint)."\n";
			}
			${$plot}->{$phase}->{'data'} .= $endPoint;
		    }
		}
	    }
	    case ("filledCurve") {
		# Draw a filled curve - using the "y2" option as the second set of y points.
		die ("GnuPlot::PrettyPlots - filledCurve requires a 'y2' vector") unless (exists($options{'y2'}));
		if ( exists($phaseRules{$phase}->{'level'}) ) {
		    # Plot just a single level, no real data.
		    my $level = "upper";
		    ${$plot}->{$phase}->{'command'} .= ${$plot}->{$phase}->{'prefix'}." '-' with filledcurve".$title
			.$lineType  {$level}
		        .$lineColor {$level}
		        .$lineWeight{$level}
		        ." fill noborder";
		    ${$plot}->{$phase}->{'data'   } .= $dummyPoint;
		    ${$plot}->{$phase}->{'data'   } .= $endPoint;
		    ${$plot}->{$phase}->{'prefix'} = ",";
		} else {
		    my $level = "upper";
		    ${$plot}->{$phase}->{'data'} .= "set style fill solid 1.0 noborder\n";
		    ${$plot}->{$phase}->{'data'} .= "plot '-' with filledcurve notitle".$lineType{$level}.$lineColor{$level}.$lineWeight{$level}." fill border\n";
		    ${$plot}->{$phase}->{'data'} .= $x->index(0)." ".$y->index(0)." ".$y->index(0)."\n";
		    for(my $iPoint=0;$iPoint<nelem($x);++$iPoint) {
			${$plot}->{$phase}->{'data'} .= $x->index($iPoint)." ".$y->index($iPoint)." ".$options{'y2'}->index($iPoint)."\n";
		    }
		    ${$plot}->{$phase}->{'data'} .= $x->index(nelem($x)-1)." ".$y->index(nelem($x)-1)." ".$y->index(nelem($x)-1)."\n";
		    ${$plot}->{$phase}->{'data'} .= $endPoint;
		}
	    }
	    case ("boxes") {
		# Check if we are asked to plot just a single level.
		if ( exists($phaseRules{$phase}->{'level'}) ) {
		    # Plot just a single level, no real data.
		    ${$plot}->{$phase}->{'command'} .= ${$plot}->{$phase}->{'prefix'}." '-' with boxes".$title
			.$lineType  {$phaseRules{$phase}->{'level'}}
		        .$lineColor {$phaseRules{$phase}->{'level'}}
		        .$lineWeight{$phaseRules{$phase}->{'level'}};
		    ${$plot}->{$phase}->{'data'   } .= $dummyPoint;
		    ${$plot}->{$phase}->{'data'   } .= $endPoint;
		    ${$plot}->{$phase}->{'prefix'} = ",";
		} else {
		    # Plot the actual data.
		    foreach my $level ( 'lower', 'upper' ) {
			${$plot}->{$phase}->{'data'} .= "set boxwidth 0.9 relative\n";
			${$plot}->{$phase}->{'data'} .= "set style fill solid 1.0\n";
			${$plot}->{$phase}->{'data'} .= "plot '-' with boxes notitle".$lineType{$level}.$lineColor{$level}.$lineWeight{$level}."\n";
			for(my $iPoint=0;$iPoint<nelem($x);++$iPoint) {
			    ${$plot}->{$phase}->{'data'} .= $x->index($iPoint)." ".$y->index($iPoint)."\n";
			}
			${$plot}->{$phase}->{'data'} .= $endPoint;
		    }
		}
	    }
	    case ( "point" ) {
		# Check if we are asked to plot just a single level.
		if ( exists($phaseRules{$phase}->{'level'}) ) {
		    # Plot just a single level, no real data.
		    ${$plot}->{$phase}->{'command'} .= ${$plot}->{$phase}->{'prefix'}." '-'".$title
			.$pointType{$phaseRules{$phase}->{'level'}}.$pointSize{$phaseRules{$phase}->{'level'}}.$lineColor{$phaseRules{$phase}->{'level'}}.$lineWeight{$phaseRules{$phase}->{'level'}};
		    ${$plot}->{$phase}->{'data'   } .= $dummyPoint;
		    ${$plot}->{$phase}->{'data'   } .= $endPoint;
		    ${$plot}->{$phase}->{'prefix'} = ",";
		} else {
		    # Plot the actual data.
		    for(my $iPoint=0;$iPoint<nelem($x);++$iPoint) {
			# Determine if vertical errors are present.
			my $showVerticalErrors = 0;
			$showVerticalErrors = 1 if (
			    ( exists($options{'errorDown' }) && $options{'errorDown' }->index($iPoint) > 0.0 )
			    ||
			    ( exists($options{'errorUp'   }) && $options{'errorUp'   }->index($iPoint) > 0.0 )
			    );
			my $showHorizontalErrors = 0;
			$showHorizontalErrors = 1 if (
			    ( exists($options{'errorLeft' }) && $options{'errorLeft' }->index($iPoint) > 0.0 )
			    ||
			    ( exists($options{'errorRight'}) && $options{'errorRight'}->index($iPoint) > 0.0 )
			    );
			my $errorCommand;
			if ( $showVerticalErrors == 1 ) {
			    if ( $showHorizontalErrors == 1 ) {
				$errorCommand = " with xyerrorbars";
			    } else {
				$errorCommand = " with errorbars";
			    }
			} else {
			    if ( $showHorizontalErrors == 1 ) {
				$errorCommand = " with xerrorbars";
			    } else {
				$errorCommand = "";
			    }
			}
			
			# Add error bar data.
			my $errors = "";
			if ( $showHorizontalErrors == 1 ) {
			    if ( exists($options{'errorLeft'}) && $options{'errorLeft'}->index($iPoint) > 0.0 ) {
				# Add a standard leftware error bar.
				my $errorPosition = $x->index($iPoint)-$options{'errorLeft'}->index($iPoint);
				$errors .= " ".$errorPosition;
			    } else  {
				# No leftward error bar.
				$errors .= " ".$x->index($iPoint);
			    }
			    if ( exists($options{'errorRight'}) && $options{'errorRight'}->index($iPoint) > 0.0 ) {
				# Add a standard rightward error bar.
				my $errorPosition = $x->index($iPoint)+$options{'errorRight'}->index($iPoint);
				$errors .= " ".$errorPosition;
			    } else  {
				# No rightward error bar.
				$errors .= " ".$x->index($iPoint);
			    }
			}
			if ( $showVerticalErrors == 1 ) {
			    if ( exists($options{'errorDown'}) && $options{'errorDown'}->index($iPoint) > 0.0 ) {
				# Add a standard downward error bar.
				my $errorPosition = $y->index($iPoint)-$options{'errorDown'}->index($iPoint);
				$errors .= " ".$errorPosition;
			    } else  {
				# No downward error bar.
				$errors .= " ".$y->index($iPoint);
			    }
			    if ( exists($options{'errorUp'}) && $options{'errorUp'}->index($iPoint) > 0.0 ) {
				# Add a standard upward error bar.
				my $errorPosition = $y->index($iPoint)+$options{'errorUp'}->index($iPoint);
				$errors .= " ".$errorPosition;
			    } else  {
				# No upward error bar.
				$errors .= " ".$y->index($iPoint);
			    }
			}

			# Add arrows.
			my $arrows      = "";
			my $clearArrows = "";
			if ( exists($options{'errorLeft'}) && $options{'errorLeft'}->index($iPoint) < 0.0 ) {
			    my $tipX = $x->index($iPoint)+$options{'errorLeft'}->index($iPoint);
			    $arrows .= "set arrow from first ".$x->index($iPoint).",".$y->index($iPoint)." to first "
				.$tipX.",".$y->index($iPoint)." filled back ".$lineColor{'lower'}.$lineWeight{'upper'}."\n";
			    $clearArrows = "unset arrow\n";
			}
			if ( exists($options{'errorRight'}) && $options{'errorRight'}->index($iPoint) < 0.0 ) {
			    my $tipX = $x->index($iPoint)-$options{'errorRight'}->index($iPoint);
			    $arrows .= "set arrow from first ".$x->index($iPoint).",".$y->index($iPoint)." to first "
				.$tipX.",".$y->index($iPoint)." filled back ".$lineColor{'lower'}.$lineWeight{'upper'}."\n";
			    $clearArrows = "unset arrow\n";
			}
			if ( exists($options{'errorDown'}) && $options{'errorDown'}->index($iPoint) < 0.0 ) {
			    my $tipY = $y->index($iPoint)+$options{'errorDown'}->index($iPoint);
			    $arrows .= "set arrow from first ".$x->index($iPoint).",".$y->index($iPoint)." to first "
				.$x->index($iPoint).",".$tipY." filled back ".$lineColor{'lower'}.$lineWeight{'upper'}."\n";
			    $clearArrows = "unset arrow\n";
			}
			if ( exists($options{'errorUp'}) && $options{'errorUp'}->index($iPoint) < 0.0 ) {
			    my $tipY = $y->index($iPoint)-$options{'errorUp'}->index($iPoint);
			    $arrows .= "set arrow from first ".$x->index($iPoint).",".$y->index($iPoint)." to first "
				.$x->index($iPoint).",".$tipY." filled back ".$lineColor{'lower'}.$lineWeight{'upper'}."\n";
			    $clearArrows = "unset arrow\n";
			}

			# Output the point.
			foreach my $level ( 'lower', 'upper' ) {
			    ${$plot}->{$phase}->{'data'} .= $arrows if ( $level eq "lower" );
			    ${$plot}->{$phase}->{'data'} .= "plot '-' notitle".$errorCommand.$pointType{$level}.$pointSize{$level}.$lineColor{$level}.$lineWeight{$level}."\n";
			    ${$plot}->{$phase}->{'data'} .= $x->index($iPoint)." ".$y->index($iPoint).$errors."\n";
			    ${$plot}->{$phase}->{'data'} .= $endPoint;
			    ${$plot}->{$phase}->{'data'} .= $clearArrows if ( $level eq "lower" );
			    # Clear any error bars after lower layer plot.
			    $errorCommand = "";
			    $errors       = "";	
			}
		    }
		}
	    }
	}
    }
}

# Subroutine to output plotting commands to GnuPlot for a set of accumulated datasets.
sub Plot_Datasets {
    my $gnuPlot = shift;
    my $plot    = shift;
    print $gnuPlot "set multiplot\n";
    foreach my $phase ( 'keyLower', 'keyUpper' ) {
	print $gnuPlot ${$plot}->{$phase}->{'command'}."\n";
	print $gnuPlot ${$plot}->{$phase}->{'data'   };
	# Switch off borders and tics after the first plot.
	print $gnuPlot "unset border; unset xtics; unset ytics; unset x2tics; unset y2tics\n" if ( $phase eq "keyLower" );
    }
    print $gnuPlot ${$plot}->{'data'}->{'data'};
    print $gnuPlot "unset multiplot\n";
    undef(${$plot});
}
