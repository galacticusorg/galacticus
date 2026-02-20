#!/bin/sh

# Build the Galacticus documentation.
# Andrew Benson (20-February-2011)

# Set defaults.
PPN=1
FORCE=yes
SUFFIX=
DIR=./work/build/
CLEAN=no
OPTIMIZE=no

# Get options.
while getopts ":p:f:d:s:c:o:" option; do
case "${option}"
in
f) FORCE=${OPTARG};;
p) PPN=${OPTARG};;
d) DIR=${OPTARG};;
s) SUFFIX=${OPTARG};;
c) CLEAN=${OPTARG};;
o) OPTIMIZE=${OPTARG};;
\?) echo "Invalid option: $OPTARG";;
:) echo "Invalid option: $OPTARG requires an argument";;
esac
done

# Export build directory.
export BUILDPATH=$DIR

# Clear out old build files.
if [ "$FORCE" = "yes" ]; then
    rm -f                                                                                                                                     \
       doc/physics/*.tex doc/inputParameters/*.tex doc/enumerations/definitions/*.tex doc/enumerations/specifiers/*.tex doc/contributions.tex \
       doc/source_documentation.tex doc/dataEnumerationSpecifiers.tex doc/dataEnumerations.tex doc/dataMethods.tex
    rm -rf $DIR
fi
# Ensure that nodeComponent and treeNode objects are built, along with any functions.
make -j$PPN GALACTICUS_BUILD_DOCS=yes SUFFIX=$SUFFIX BUILDPATH=$DIR all
if [ $? -ne 0 ]; then
 echo Failed to build all executables
 exit 1
fi

# Optimize storage if requested.
if [ "$OPTIMIZE" = "yes" ]; then
    rm -f $DIR/*.o $DIR/*.md5.blob
    rm -f *.exe$SUFFIX
fi

# Extract source code data.
./scripts/doc/extractData.pl source doc/data
if [ $? -ne 0 ]; then
 echo Failed to extract source code data
 exit 1
fi

# Extract contributor data.
./scripts/doc/Extract_Contributors.pl . doc/contributions.tex
if [ $? -ne 0 ]; then
 echo Failed to extract contributor data
 exit 1
fi

# Analyze source code.
./scripts/doc/Code_Analyzer.pl source doc/source_documentation.tex
if [ $? -ne 0 ]; then
 echo Failed to analyze source code
 exit 1
fi

# Extract constants data.
./scripts/doc/constants.pl $DIR doc/constants.tex
if [ $? -ne 0 ]; then
 echo Failed to extract constants data
 exit 1
fi

# Move to the documentation folder.
cd doc

# Demangle the bibliography.
./Bibliography_Demangle.pl
if [ $? -ne 0 ]; then
 echo Failed to demangle bibliography
 exit 1
fi

# Order physics descriptions.
ls physics/*.tex | sort | awk '{print "\\input{"substr($1,1,length($1)-4)"}"}' > autoPhysics.tex

# Order enumeration definitions.
ls enumerations/definitions/*.tex | sort | awk '{print "\\input{"substr($1,1,length($1)-4)"}"}' > autoEnumerationDefinitions.tex

# Iterate over manuals.
for type in "Usage" "Physics" "Development" "Source"; do

    # Compile the manuals.
    iPass=1
    while [ $iPass -le 6 ]; do
	# Run pdflatex.
	if [ $iPass -le 5 ]; then
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d > /dev/null 2>&1
	else
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
	fi
	if [ $? -ne 0 ]; then
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
	    echo pdflatex failed
	    exit 1
	fi

	# Run bibtex.
	if [ $iPass -le 5 ]; then
	    bibtex Galacticus_$type > /dev/null 2>&1
	else
	    bibtex Galacticus_$type
	fi
	if [ $? -ne 0 ]; then
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
	    bibtex Galacticus_$type
	    echo bibtex failed
	    exit 1
	fi

	# Run makeindex.
	if [ $iPass -le 5 ]; then
	    makeindex Galacticus_$type > /dev/null 2>&1
	else
	    makeindex Galacticus_$type
	fi
	if [ $? -ne 0 ]; then
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
	    makeindex Galacticus_$type
	    echo makeindex failed for main index
	    exit 1
	fi

	# Run makeglossaries.
	if [ $iPass -le 5 ]; then
	    makeglossaries Galacticus_$type > /dev/null 2>&1
	else
	    makeglossaries Galacticus_$type
	fi
	if [ $? -ne 0 ]; then
	    pdflatex Galacticus_$type | grep -v -i -e overfull -e underfull | sed -r /'^$'/d | sed -r /'\[[0-9]*\]'/d
	    makeglossaries Galacticus_$type
	    echo make glossaries failed
	    exit 1
	fi
	
	iPass=$((iPass+1))
    done

done

# Clean build files if requested.
if [ "$CLEAN" = "yes" ]; then
    cd ..
    rm -rf $DIR
    rm -f *.exe$SUFFIX
fi

exit 0
