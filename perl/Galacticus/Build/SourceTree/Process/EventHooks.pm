# Contains a Perl module which implements processing of event directives.

package Galacticus::Build::SourceTree::Process::EventHooks;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Directives;
use Text::Template 'fill_in_string';
use Digest::MD5 qw(md5_hex);

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'eventHooks'} = \&Process_EventHooks;

sub Process_EventHooks {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get code directive locations.
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
    # Walk the tree, looking for hook directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Handle eventHookManger directives by building eventHook objects for all events.
	if ( $node->{'type'} eq "eventHookManager" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Find all hook directives.
	    my @hooks = map {&Galacticus::Build::Directives::Extract_Directives($_,'eventHook')} &List::ExtraUtils::as_array($directiveLocations->{'eventHook'}->{'file'});
	    # Create an object for each event hook.
	    foreach my $hook ( @hooks ) {
		# Determine the interface type for this hook.
		$code::interfaceType = &interfaceTypeGet($hook);
		my $hookObject;
		unless ( $code::interfaceType eq "Unspecified" ) {
		    # Parse the interface definition.
		    $code::declarations = $hook->{'interface'};
		    @code::arguments    = ();
		    open(my $declarations,"<",\$hook->{'interface'});
		    while ( ! eof($declarations) ) {
			&Fortran::Utils::Get_Fortran_Line($declarations,my $rawLine, my $processedLine, my $bufferedComments);
			foreach my $declarator ( keys(%Fortran::Utils::intrinsicDeclarations) ) {
			    if ( my @matches = ( $processedLine =~ $Fortran::Utils::intrinsicDeclarations{$declarator}->{'regEx'} ) ) {
				push(@code::arguments,&Fortran::Utils::Extract_Variables($matches[$Fortran::Utils::intrinsicDeclarations{$declarator}->{'variables'}],keepQualifiers => 0));
			    }
			}
		    }
		    close($declarations);
		    # Build the required types and functions.
		    $hookObject = fill_in_string(<<'CODE', PACKAGE => 'code');

type, extends(hook) :: hook{$interfaceType}
   procedure(interface{$interfaceType}), pointer, nopass :: function_ => null()
end type hook{$interfaceType}
 
type, extends(eventHook) :: eventHook{$interfaceType}
  private
 contains
  procedure :: attach => eventHook{$interfaceType}Attach
end type eventHook{$interfaceType}

abstract interface
 subroutine interface{$interfaceType}(self{scalar(@arguments) > 0 ? ",".join(",",@arguments) : ""})
  class(*), intent(inout) :: self
{$declarations}
 end subroutine interface{$interfaceType}
end interface
CODE
		    my $attacher = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventHook{$interfaceType}Attach(self,object_,function_)
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)          :: self
  class    (*                        ), intent(in   ), target  :: object_
  procedure(interface{$interfaceType})                         :: function_
  class    (hook                     )               , pointer :: hook_

  if (associated(self%first_)) then
     hook_ => self%first_
     do while (associated(hook_%next))
        hook_ => hook_%next
     end do
     allocate(hook{$interfaceType} :: hook_%next )
     hook_ => hook_%next
  else
     allocate(hook{$interfaceType} :: self%first_)
     hook_ => self%first_
  end if
  select type (hook_)
  type is (hook{$interfaceType})
     hook_%object_   => object_
     hook_%function_ => function_
  end select
  self%count_=self%count_+1
  return
end subroutine eventHook{$interfaceType}Attach
CODE
		    my $newNode   =
		    {
			type       => "code",
			content    => $attacher,
			firstChild => undef(),
			source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
			line       => 1
		    };
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$newNode]);
		    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},"hook".$code::interfaceType,"public");
		}
		 $hookObject .= "type(eventHook".$code::interfaceType."), public :: ".$hook->{'name'}."Event\n";
		my $newNode   =
		{
		    type       => "code",
		    content    => $hookObject,
		    firstChild => undef(),
		    source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
		    line       => 1
		};
		&Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    }
	}
	# Handle eventHook directives by creating code to call any hooked functions.
	if ( $node->{'type'} eq "eventHook" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} =  1;
	    # Insert the module.
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    Events_Hooks =>
		    {
			intrinsic => 0,
			all       => 1
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Insert required variables.
	    my @declarations =
		(
		 {
		     intrinsic  => "class"       ,
		     type       => "hook"       ,
		     variables  => [ "hook_"   ],
		     attributes => [ "pointer" ]
		 }
		);
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
	    # Create the code.
	    $code::interfaceType = &interfaceTypeGet($node->{'directive'});
	    $code::callWith      = exists($node->{'directive'}->{'callWith'}) ? ",".$node->{'directive'}->{'callWith'} : "";
	    $code::eventName     = $node->{'directive'}->{'name'    };
	    my $eventHookCode    = fill_in_string(<<'CODE', PACKAGE => 'code');
hook_ => {$eventName}Event%first()
do while (associated(hook_))
   select type (hook_)
   type is (hook{$interfaceType})
     call hook_%function_(hook_%object_{$callWith})
   end select
   hook_ => hook_%next
end do
CODE
	    # Insert the code.
	    my $newNode =
	    {
		type       => "code",
		content    => $eventHookCode,
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub interfaceTypeGet {
    my $hook = shift();
    my $interfaceType;
    if ( exists($hook->{'interface'}) ) {
	$interfaceType = md5_hex($hook->{'name'});
   } else {
	$interfaceType = "Unspecified";
    }
    return $interfaceType;
}

1;
