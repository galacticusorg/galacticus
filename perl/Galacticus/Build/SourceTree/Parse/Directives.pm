# Contains a Perl module which implements parsing of directives in the Galacticus preprocessor system.

package Galacticus::Build::SourceTree::Parse::Directives;
use strict;
use warnings;
use utf8;
use Data::Dumper;
use XML::Simple;
use XML::LibXML;
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks      {'directives'} = \&Parse_Directives;
$Galacticus::Build::SourceTree::Hooks::postprocessHooks{'directives'} = \&PostProcess_Directives;

sub Parse_Directives {
    # Get the tree.
    my $tree = shift();
    # Initialize state storables database.
    our $stateStorables;
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get state storables database if we do not have it.
    $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml")
	if ( ! $stateStorables && exists($ENV{'BUILDPATH'}) && -e $ENV{'BUILDPATH'}."/stateStorables.xml" );
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Find code blocks.
	if ( $node->{'type'} eq "code" ) {
	    # Initialize a set of new nodes.
	    my @newNodes;
	    # Read the code block, accumulating directives as we go.
	    my $rawCode;
	    my $rawOpener;
	    my $rawDirective;
	    my $strippedDirective;
	    my $inDirective      = 0;
	    my $inXML            = 0;
	    my $directiveRoot;
	    my $lineNumber       = exists($node->{'line'  }) ? $node->{'line'  } : 0        ;
	    my $source           = exists($node->{'source'}) ? $node->{'source'} : "unknown";
	    my $rawCodeLine      = $lineNumber;
	    my $rawDirectiveLine = $lineNumber;
	    open(my $code,"<",\$node->{'content'});
	    while ( my $line = <$code> ) {
		# Detect the end of an XML section and change state.
		$inXML = 0
		    if ( $line =~ m/^\s*!!\]/ );
		# Process XML blocks.
		my $isDirective  = 0;
		my $endDirective = 0;
		# Strip the initial "!<" and any following whitespace (but not tabs - this allows us to use tabs for formatting
		# purposes).
		(my $strippedLine = $line) =~ s/^\s*\!<\s*//;
		if ( $inXML ) {
		    # Determine if line is a directive line.
		    $isDirective    = 1
			if ( $strippedLine =~ m/^\s*\<([^\s\>\/]+)/ || $inDirective == 1 );
		    $directiveRoot = $1
			if ( $isDirective == 1 && $inDirective == 0 );		
		    # Catch the end of directives.
		    $endDirective = 1
			if ( $isDirective == 1 && $strippedLine =~ m/\s*\<\/$directiveRoot\>/ );
		    $endDirective = 1
			if ( $isDirective == 1 && $inDirective == 0 && ( $strippedLine =~ m/\s*\<$directiveRoot\s.*\/\>/ || $strippedLine =~ m/\s*\<$directiveRoot\/\>/ ) );
		    # Record whether we are currently in or out of a directive.
		    $inDirective = 1
			if ( $isDirective == 1 );
		}
		# Accumulate raw text.
		if ( $inDirective ) {
		    # Process non-breaking spaces as a special case.
		    $strippedLine =~ s/&nbsp;/ /g;
		    # Accumulate the line.
		    $rawDirective      .= $line;
		    $strippedDirective .= $strippedLine;
		} elsif ( $line !~ m/^\s*!!(\[|\])/ ) {
		    $rawCode           .= $line;
		} elsif ( $line =~ m/^\s*!!\[/ ) {
		    $rawOpener          = $line;
		} else {
		    $rawCodeLine      = $lineNumber+1;
		    $rawDirectiveLine = $lineNumber+1;
		}
		# Process code and directive blocks as necessary.
		if ( ( $inDirective == 1 || eof($code) ) && $rawCode      ) {
		    # Create a new node.
		    my $newNode =
		    {
			type       => "code"            ,
			content    => $rawCode          ,
			firstChild => undef(),
			source     => $source,
			line       => $rawCodeLine
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    push(
			@newNodes,
			$newNode
			);
		    # Reset the raw code text.
		    undef($rawCode);
		    $rawCodeLine      = $lineNumber;
		    $rawDirectiveLine = $lineNumber;
		}
		if ( ( $inDirective == 0 || eof($code) || $endDirective ) && $rawDirective ) {
		    # Attempt to parse the directive XML.
		    my $directive = eval{$xml->XMLin($strippedDirective, keepRoot => 1)};
		    if ( $@ ) {
			print "Parse_Directives: failed parsing with message:\n".$@."\n";
			print $strippedDirective;
			die();
		    }
		    my $directiveName = (keys %{$directive})[0];
		    # Validate the directive if possible.
		    my $schema;
		    if ( $stateStorables && exists($stateStorables->{'functionClasses'}{$directiveName."Class"}) ) {
			# functionClass instances - create a custom schema for this instance.
			$functionClass::name = $directiveName;
			my $schemaDocument = fill_in_string(<<'SCHEMA', PACKAGE => 'functionClass');
<?xml version="1.0"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xs:element name="{$name}">
    <xs:complexType>
      <xs:sequence>
       <xs:element name="description"       type="xs:string" minOccurs="1" maxOccurs="1"/>
       <xs:element name="descriptorSpecial" type="xs:string" minOccurs="0" maxOccurs="1"/>
       <xs:element name="linkedList"                         minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:attribute name="type"       use="required"/>
         <xs:attribute name="variable"   use="required"/>
         <xs:attribute name="next"       use="required"/>
         <xs:attribute name="object"     use="required"/>
         <xs:attribute name="objectType" use="required"/>
         <xs:attribute name="module"     use="optional"/>
        </xs:complexType>
       </xs:element>
       <xs:element name="deepCopy"                     minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="ignore"        minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="functionClass" minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="increment"     minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="atomic"    use="optional" >
             <xs:simpleType>
              <xs:restriction base="xs:string">
               <xs:enumeration value="no" />
               <xs:enumeration value="yes"/>
              </xs:restriction>
             </xs:simpleType>
            </xs:attribute>
           </xs:complexType>
          </xs:element>
          <xs:element name="setTo"         minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="value"     use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
        </xs:complexType>
       </xs:element>
       <xs:element name="assignment"                   minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="functionClass" minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
         <xs:attribute name="forceArrayAssign" use="optional"/>
        </xs:complexType>
       </xs:element>
       <xs:element name="stateStorable"                minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="functionClass" minOccurs="0" maxOccurs="1"         >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="restoreTo"     minOccurs="0" maxOccurs="unbounded" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="state"     use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="exclude"       minOccurs="0" maxOccurs="1"         >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
        </xs:complexType>
       </xs:element>
       <xs:element name="stateStore"                   minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="stateStore" minOccurs="1" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="store"     use="required"/>
            <xs:attribute name="restore"   use="required"/>
            <xs:attribute name="module"    use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
        </xs:complexType>
       </xs:element>
       <xs:element name="runTimeFileDependencies"      minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:attribute name="paths" use="required"/>
        </xs:complexType>
       </xs:element>
      </xs:sequence>
      <xs:attribute name="name"      use="required"/>
      <xs:attribute name="recursive" use="optional" >
       <xs:simpleType>
        <xs:restriction base="xs:string">
         <xs:enumeration value="no" />
         <xs:enumeration value="yes"/>
        </xs:restriction>
       </xs:simpleType>
      </xs:attribute>
      <xs:attribute name="abstract"  use="optional" >
       <xs:simpleType>
        <xs:restriction base="xs:string">
         <xs:enumeration value="no" />
         <xs:enumeration value="yes"/>
        </xs:restriction>
       </xs:simpleType>
      </xs:attribute>
    </xs:complexType>
  </xs:element>
</xs:schema>
SCHEMA
			$schema = XML::LibXML::Schema->new( string => $schemaDocument );
		    } elsif ( $stateStorables && grep {$_ eq $directiveName} @{$stateStorables->{'eventHookStatics'}} ) {
			# eventHookStatic instances - create a custom schema for this instance.
			$eventHookStatic::name = $directiveName;
			my $schemaDocument = fill_in_string(<<'SCHEMA', PACKAGE => 'eventHookStatic');
<?xml version="1.0"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xs:simpleType name="yesNo">
    <xs:restriction base="xs:string">
      <xs:enumeration value="yes"/>
      <xs:enumeration value="no" />
    </xs:restriction>
  </xs:simpleType> 
  <xs:element name="{$name}">
    <xs:complexType>
      <xs:attribute name="function"  use="required"             />
      <xs:attribute name="after"     use="optional"             />
      <xs:attribute name="before"    use="optional"             />
      <xs:attribute name="useGlobal" use="optional" type="yesNo"/>
    </xs:complexType>
  </xs:element>
</xs:schema>
SCHEMA
			$schema = XML::LibXML::Schema->new( string => $schemaDocument );
		    } elsif ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/schema/".$directiveName.".xsd" ) {
			$schema = XML::LibXML::Schema->new( location =>  $ENV{'GALACTICUS_EXEC_PATH'}."/schema/".$directiveName.".xsd");
		    }
		    if ( $schema ) {
			my $document = XML::LibXML->load_xml( string => $strippedDirective );
			eval { $schema->validate( $document ) };
			if ( $@ ) {
			    my $nodeParent = $node;
			    while ( $nodeParent->{'type'} ne "file" ){
				$nodeParent = $nodeParent->{'parent'};
			    }
			    die "Galacticus::Build::SourceTree::Parse::Directives(): validation of directive '".$directiveName."' failed in file ".$nodeParent->{'name'}." at line ".$node->{'line'}.":\n".$@;
			}
		    }
		    # Create a new node.
		    my $newNode =
		    {
			type       => $directiveName              ,
			directive  => $directive->{$directiveName},
			line       => $node     ->{'line'        },
			processed  => 0                           ,
			source     => $source                     ,
			line       => $rawDirectiveLine
		    };
		    (my $rawCloser = $rawOpener) =~ s/\[/\]/;
		    $rawDirective = $rawOpener.$rawDirective.$rawCloser;
		    $newNode->{'firstChild'} =
		    {
			type       => "code"           ,
			content    => $rawDirective    ,
			parent     => $newNode         ,
			sibling    => undef()          ,
			firstChild => undef()          ,
			source     => $source          ,
			line       => $rawDirectiveLine
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    # If the end of the code has been reached and we're in a code block, pop that code block from the children
		    # array before we push our directive node.
		    my $codeNode = pop(@newNodes)
			if ( eof($code) && $isDirective == 0 );
		    push(
			@newNodes,
			$newNode
			);
		    push(@newNodes,$codeNode)
			if ( eof($code) && $isDirective == 0 );
		    # Reset the raw directive text.
		    $inDirective = 0;
		    undef($rawDirective     );
		    undef($strippedDirective);
		    $rawCodeLine      = $lineNumber;
		    $rawDirectiveLine = $lineNumber;
		}
		# Detect the start of an XML section and change state.
		$inXML = 1
		    if ( $line =~ m/^\s*!!\[/ );
		# Increment line number count.
		++$lineNumber;
	    }
	    close($code);
	    # If we have a single code block, nothing needs to change.
	    unless ( scalar(@newNodes) == 1 && $newNodes[0]->{'type'} eq "code" ) {
		# New nodes created, insert them, replacing the old node.
		&Galacticus::Build::SourceTree::ReplaceNode($node,\@newNodes);
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }    
}

sub PostProcess_Directives {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Find directives.
	if ( exists($node->{'directive'}) ) {
	    die("directive '".$node->{'type'}."' was not processed at line ".$node->{'line'}." in ".$tree->{'name'})
		unless ( exists($node->{'directive'}->{'processed'}) && $node->{'directive'}->{'processed'} );
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }    
}

1;
