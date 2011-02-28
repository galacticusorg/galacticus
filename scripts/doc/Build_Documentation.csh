#!/bin/csh

# Build the Galacticus documentation.
# Andrew Benson (20-February-2011)

# Extract source code data.
scripts/doc/Extract_Data.pl source doc/data
if ( $? != 0 ) then
 echo Failed to extract source code data
 exit 1
endif

# Analyze source code.
scripts/doc/Code_Analyzer.pl source doc/source_documentation.tex
if ( $? != 0 ) then
 echo Failed to analyze source code
 exit 1
endif

# Move to the documentation folder.
cd doc

# Demangle the bibliography.
Bibliography_Demangle.pl
if ( $? != 0 ) then
 echo Failed to demangle bibliography
 exit 1
endif

# Compile the manual.
@ iPass = 1
while( $iPass <= 3 )
 # Run pdflatex.
 pdflatex Galacticus | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
 if ( $? != 0 ) then
  echo pdflatex failed
  exit 1
 endif

 # Run bibtex.
 bibtex Galacticus
 if ( $? != 0 ) then
  echo bibtex failed
  exit 1
 endif

 # Run makeindex.
 makeindex -s Galacticus.ist Galacticus
 if ( $? != 0 ) then
  echo makeindex failed
  exit 1
 endif

 @ iPass++
end

exit 0
