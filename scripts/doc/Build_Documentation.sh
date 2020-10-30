#!/bin/sh

# Build the Galacticus documentation.
# Andrew Benson (20-February-2011)

# Set defaults.
PPN=1

# Get options.
while getopts p: option
do
case "${option}"
in
p) PPN=${OPTARG};;
esac
done

# Clear out old build files.
rm -f                                                                                                                                     \
   doc/physics/*.tex doc/inputParameters/*.tex doc/enumerations/definitions/*.tex doc/enumerations/specifiers/*.tex doc/contributions.tex \
   doc/source_documentation.tex doc/dataEnumerationSpecifiers.tex doc/dataEnumerations.tex doc/dataMethods.tex

# Ensure that nodeComponent and treeNode objects are built, along with any functions.
rm -rf work/build
make -j$PPN all
if [ $? -ne 0 ]; then
 echo Failed to build all executables
 exit 1
fi

# Extract source code data.
scripts/doc/extractData.pl source doc/data
if [ $? -ne 0 ]; then
 echo Failed to extract source code data
 exit 1
fi

# Extract contributor data.
scripts/doc/Extract_Contributors.pl . doc/contributions.tex
if [ $? -ne 0 ]; then
 echo Failed to extract contributor data
 exit 1
fi

# Analyze source code.
scripts/doc/Code_Analyzer.pl source doc/source_documentation.tex
if [ $? -ne 0 ]; then
 echo Failed to analyze source code
 exit 1
fi

# Move to the documentation folder.
cd doc

# Demangle the bibliography.
Bibliography_Demangle.pl
if [ $? -ne 0 ]; then
 echo Failed to demangle bibliography
 exit 1
fi

# Order physics descriptions.
ls physics/*.tex | sort | awk '{print "\\input{"substr($1,1,length($1)-4)"}"}' > autoPhysics.tex

# Order enumeration definitions.
ls enumerations/definitions/*.tex | sort | awk '{print "\\input{"substr($1,1,length($1)-4)"}"}' > autoEnumerationDefinitions.tex

# Order enumeration specifiers.
ls enumerations/specifiers/*.tex | sort | awk '{print "\\input{"substr($1,1,length($1)-4)"}"}' > autoEnumerationSpecifiers.tex

# Iterate over manuals.
for type in "Usage" "Physics" "Development" "Source"; do

    # Compile the manuals.
    iPass=1
    while [ $iPass -le 6 ]; do
	# Run pdflatex.
	if [ $iPass -le 5 ]; then
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d >& /dev/null
	else
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
	fi
	if [ $? -ne 0 ]; then
	    echo pdflatex failed
	    exit 1
	fi

	# Run bibtex.
	if [ $iPass -le 5 ]; then
	    bibtex Galacticus_$type >& /dev/null
	else
	    bibtex Galacticus_$type
	fi
	if [ $? -ne 0 ]; then
	    echo bibtex failed
	    exit 1
	fi

	# Run makeindex.
	if [ $iPass -le 5 ]; then
	    makeindex Galacticus_$type >& /dev/null
	else
	    makeindex Galacticus_$type
	fi
	if [ $? -ne 0 ]; then
	    echo makeindex failed for main index
	    exit 1
	fi

	# Run makeglossaries.
	if [ $iPass -le 5 ]; then
	    makeglossaries Galacticus_$type >& /dev/null
	else
	    makeglossaries Galacticus_$type
	fi
	if [ $? -ne 0 ]; then
	    echo make glossaries failed
	    exit 1
	fi
	
	iPass=$((iPass+1))
    done

done

exit 0
