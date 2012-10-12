# Contains a Perl module which implements processing GnuPlot output through LaTeX into either PDF or PNG formats.

package LaTeX;
use strict;
use warnings;
use File::Copy;
use Digest::MD5  qw(md5_hex);
require System::Redirect;

sub GnuPlot2PNG {
    # Generate a PNG image of the plot with a transparent background.

    # Get the name of the GnuPlot-generated EPS file.
    my $gnuplotEpsFile = shift;
    my (%options) = @_ if ( $#_ >= 1 );
    
    # Get the root name.
    (my $gnuplotRoot = $gnuplotEpsFile) =~ s/\.eps//;
    
    # Set background color, assuming white by default.
    my @rgbFractional = [  1,  1,  1];
    my @rgb           = [255,255,255];
    if ( exists($options{'backgroundColor'}) ) {
	$options{'backgroundColor'} =~ s/^\#//;
	@rgbFractional = map {hex($_)/255.0 } unpack 'a2a2a2', $options{'backgroundColor'};
	@rgb           = map {hex($_)       } unpack 'a2a2a2', $options{'backgroundColor'};
    }

    # Edit the LaTeX input to reverse black and white for labels.
    open(iHndl,$gnuplotRoot.".tex");
    open(oHndl,">".$gnuplotRoot.".tex.swapped");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/LTw/ ) {$line =~ s/white/black/};
	if ( $line =~ m/LTb/ ) {$line =~ s/black/white/};
	if ( $line =~ m/LTa/ ) {$line =~ s/black/white/};
	if ( $line =~ m/\\begin{picture}/ ) {
	    print oHndl "\\definecolor{myColor}{rgb}{".join(",",@rgbFractional)."}\n";
	    print oHndl "\\colorbox{myColor}{\n";
	}
	print oHndl $line;
	if ( $line =~ m/\\end{picture}/ ) {print oHndl "}\n"};
    }
    close(oHndl);
    close(iHndl);
    unlink($gnuplotRoot.".tex");
    move($gnuplotRoot.".tex.swapped",$gnuplotRoot.".tex");

    # Edit the EPS file to reverse black and white for axes.
    open(iHndl,$gnuplotRoot.".eps");
    open(oHndl,">".$gnuplotRoot.".eps.swapped");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/LCw/ ) {$line =~ s/1 1 1/0 0 0/};
	if ( $line =~ m/LCb/ ) {$line =~ s/0 0 0/1 1 1/};
	if ( $line =~ m/LCa/ ) {$line =~ s/0 0 0/1 1 1/};
	print oHndl $line;
    }
    close(oHndl);
    close(iHndl);
    unlink($gnuplotRoot.".eps");
    move($gnuplotRoot.".eps.swapped",$gnuplotRoot.".eps");

    # Ensure PNG file does not exist.
    unlink($gnuplotRoot.".png");

    # Convert to PDF.
    &GnuPlot2PDF($gnuplotEpsFile);

    # Convert to PNG.
    my $backgroundColor = "'rgb(".join(",",@rgb).")'";
    system("convert -shave 1x1 -scale 3072 -density 300 -antialias -background transparent -transparent ".$backgroundColor." ".$gnuplotRoot.".pdf ".$gnuplotRoot.".png");

    # Remove the PDF file.
    unlink($gnuplotRoot.".pdf");
}

sub GnuPlot2ODG {
    # Generate an ODG file of the plot.

    # Get the name of the GnuPlot-generated EPS file.
    my $gnuplotEpsFile = shift;
    my (%options) = @_;
    
    # Get the root name.
    (my $gnuplotRoot = $gnuplotEpsFile) =~ s/\.eps//;
    
    # Edit the LaTeX input to reverse black and white for labels.
    open(iHndl,$gnuplotRoot.".tex");
    open(oHndl,">".$gnuplotRoot.".tex.swapped");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/LTw/ ) {$line =~ s/white/black/};
	if ( $line =~ m/LTb/ ) {$line =~ s/black/white/};
	if ( $line =~ m/LTa/ ) {$line =~ s/black/white/};
	print oHndl $line;
    }
    close(oHndl);
    close(iHndl);
    unlink($gnuplotRoot.".tex");
    move($gnuplotRoot.".tex.swapped",$gnuplotRoot.".tex");

    # Edit the EPS file to reverse black and white for axes.
    open(iHndl,$gnuplotRoot.".eps");
    open(oHndl,">".$gnuplotRoot.".eps.swapped");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/LCw/ ) {$line =~ s/1 1 1/0 0 0/};
	if ( $line =~ m/LCb/ ) {$line =~ s/0 0 0/1 1 1/};
	if ( $line =~ m/LCa/ ) {$line =~ s/0 0 0/1 1 1/};
	print oHndl $line;
    }
    close(oHndl);
    close(iHndl);
    unlink($gnuplotRoot.".eps");
    move($gnuplotRoot.".eps.swapped",$gnuplotRoot.".eps");

    # Convert to PDF.
    &GnuPlot2PDF($gnuplotEpsFile,%options);

    # Convert to SVG.
    system("pdf2svg ".$gnuplotRoot.".pdf ".$gnuplotRoot."_percent.svg");

    # Convert SVG RGB colors from percentages to 0-255 scale.
    open(iHndl,$gnuplotRoot."_percent.svg");
    open(oHndl,">".$gnuplotRoot.".svg");
    while ( my $line = <iHndl> ) {
	while ( $line =~ m/rgb\(([\d\.]+)\%,([\d\.]+)\%,([\d\.]+)\%\)/ ) {
	    my $r = int($1*255.0/100.0);
	    my $g = int($2*255.0/100.0);
	    my $b = int($3*255.0/100.0);
	    $line =~ s/rgb\([\d\.]+\%,[\d\.]+\%,[\d\.]+\%\)/rgb($r,$g,$b)/;
	}
	print oHndl $line;
    }
    close(oHndl);
    close(iHndl);

    # Convert to ODG.
    my @svg2officeLocations = ( "/usr/bin", "/usr/local/bin", $ENV{'HOME'}."/bin" );
    push(
	@svg2officeLocations,
	$ENV{'GALACTICUS_ROOT_V091'}."../Tools/bin"
	)
	if ( exists($ENV{'GALACTICUS_ROOT_V091'}) );
    my $svg2office;
    foreach my $location ( @svg2officeLocations ) {
	$svg2office = $location."/svg2office-1.2.2.jar"
	    if ( -e $location."/svg2office-1.2.2.jar" );
    }
    die("GnuPlot::LaTeX.pm: svg2office-1.2.2.jar not found")
	unless ( defined($svg2office) );
    system("java -jar ".$svg2office." ".$gnuplotRoot.".svg");

    # Remove the PDF and SVG files.
    unlink($gnuplotRoot.".pdf",$gnuplotRoot.".svg",$gnuplotRoot."_percent.svg");
}

sub GnuPlot2PDF {
    # Get the name of the GnuPlot-generated EPS file.
    my $gnuplotEpsFile = shift;
    # Extract any remaining options.
    my (%options) = @_ if ( $#_ >= 0 );
    
    # Get the root name.
    (my $gnuplotRoot = $gnuplotEpsFile) =~ s/\.eps//;

    # Get any folder name.
    my $folderName = "./";
    my $gnuplotBase = $gnuplotRoot;
    if ( $gnuplotRoot =~ m/^(.*\/)(.*)$/ ) {
	$folderName = $1;
	$gnuplotBase = $2;
    }

    # Construct the name of the corresponding LaTeX files.
    my $gnuplotLatexFile = $gnuplotRoot.".tex";
    my $gnuplotAuxFile   = $gnuplotRoot.".aux";

    # Construct the name of the corresponding pdf file.
    my $gnuplotPdfFile = $gnuplotRoot.".pdf";

    # Remove duplicated labelling.
    # Initialize hash that will store MD5s of previously seen labels.
    my %labels;
    # Populate with MD5s for empty front and back text labels.
    $labels{'079019eb66826e1085e8f82634835781'} = 1;
    $labels{'47add47b938f7b849abb82e200954083'} = 1;
    # Process the file.
    open(iHndl,$gnuplotRoot.".tex");
    open(oHndl,">".$gnuplotRoot.".tex.swapped");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/^\s*\\gplgaddtomacro\\gpl(front|back)text\{\%\s*$/ ) {
	    my $buffer = $line;
	    do {
		$line    = <iHndl>;
		$buffer .= $line;
	    } until ( $line =~ m/^\s*\}\%\s*$/ );
	    my $md5 = md5_hex($buffer);
	    $line = "";
	    $line = $buffer
		unless ( exists($labels{$md5}) );
	    $labels{$md5} = 1;
	}
	$line =~ s/includegraphics\{$folderName/includegraphics\{/;
	print oHndl $line;
    }
    close(oHndl);
    close(iHndl);
    unlink($gnuplotRoot.".tex");
    move($gnuplotRoot.".tex.swapped",$gnuplotRoot.".tex");

    # Convert the plot body from EPS to PDF.
    &SystemRedirect::tofile("epstopdf ".$gnuplotEpsFile,"/dev/null");

    # Get the dimensions of the plot body.
    my $width;
    my $height;
    open(pHndl,"pdfinfo ".$gnuplotPdfFile."|");
    while ( my $line = <pHndl> ) {
	if ( $line =~ m/Page size:\s*(\d+)\s*x\s*(\d+)/ ) {
	    $width  = $1+100;
	    $height = $2+100;
	}
    }
    close(pHndl);

    # Create a wrapper file for the LaTeX.
    my $fontSize = "10";
    $fontSize = $options{'fontSize'}
        if ( exists($options{'fontSize'}) );
    my $wrapper = "gnuplotWrapper".$$;
    open(wHndl,">".$folderName.$wrapper.".tex");
    print wHndl "\\documentclass[".$fontSize."pt]{article}\n";
    print wHndl "\\usepackage[margin=0pt";
    print wHndl ",paperwidth=".$width."pt,paperheight=".$height."pt"
	if ( defined($width) );
    print wHndl "]{geometry}\n";
    print wHndl "\\usepackage{graphicx}\n\\usepackage{nopageno}\n\\usepackage{txfonts}\n\\usepackage[usenames]{color}\n\\begin{document}\n\\include{".$gnuplotBase."}\n\\end{document}\n";
    close(wHndl);
    &SystemRedirect::tofile("cd ".$folderName."; pdflatex ".$wrapper."; pdfcrop ".$wrapper.".pdf","/dev/null");
    move($folderName.$wrapper."-crop.pdf",$gnuplotPdfFile);
    unlink(
    	$folderName.$wrapper.".pdf",
    	$folderName.$wrapper.".tex",
    	$folderName.$wrapper.".log",
    	$folderName.$wrapper.".aux",
    	$gnuplotEpsFile,
    	$gnuplotLatexFile,
    	$gnuplotAuxFile
    	)
	unless ( exists($options{'cleanUp'}) && $options{'cleanUp'} == 0 );
}

1;
