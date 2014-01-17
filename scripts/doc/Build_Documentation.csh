#!/usr/bin/env tcsh

# Build the Galacticus documentation.
# Andrew Benson (20-February-2011)

# Ensure that nodeComponent and treeNode objects are built, along with any functions.
make all
if ( $? != 0 ) then
 echo Failed to build all executables
 exit 1
endif

# Extract source code data.
scripts/doc/Extract_Data.pl source doc/data
if ( $? != 0 ) then
 echo Failed to extract source code data
 exit 1
endif

# Extract contributor data.
scripts/doc/Extract_Contributors.pl source doc/contributions.tex
if ( $? != 0 ) then
 echo Failed to extract contributor data
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

# Order method descriptions.
ls methods/*.tex | sort | awk '{print "\\input{"substr($1,1,length($1)-4)"}"}' > autoMethods.tex

# Compile the manual.
@ iPass = 1
while( $iPass <= 6 )
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
 makeindex Galacticus
 if ( $? != 0 ) then
  echo makeindex failed for main index
  exit 1
 endif

 # Run makeindex for code index.
 makeindex -s Galacticus.isty Galacticus.cdx -o Galacticus.cnd
 if ( $? != 0 ) then
  echo makeindex failed for code index
  exit 1
 endif

 # Run makeglossaries.
 makeglossaries Galacticus
 if ( $? != 0 ) then
  echo make glossaries failed
  exit 1
 endif

 @ iPass++
end

exit 0
