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
  !@ <objectMethods>
  !@   <object>eventHook{$interfaceType}</object>
  !@   <objectMethod>
  !@     <method>attach</method>
  !@     <type>\void</type>
  !@     <arguments>\textcolor\{red\}\{\textless class(*)\textgreater\} *object\_\argin, \textcolor\{red\}\{\textless procedure()\textgreater\} *function\_\argin</arguments>
  !@     <description>Attach a hook to the event.</description>
  !@   </objectMethod>
  !@   <objectMethod>
  !@     <method>isAttached</method>
  !@     <type>\logicalzero</type>
  !@     <arguments>\textcolor\{red\}\{\textless class(*)\textgreater\} *object\_\argin, \textcolor\{red\}\{\textless procedure()\textgreater\} *function\_\argin</arguments>
  !@     <description>Return true if the object is attached to this event.</description>
  !@   </objectMethod>
  !@   <objectMethod>
  !@     <method>detach</method>
  !@     <type>\void</type>
  !@     <arguments>\textcolor\{red\}\{\textless class(*)\textgreater\} *object\_\argin, \textcolor\{red\}\{\textless procedure()\textgreater\} *function\_\argin</arguments>
  !@     <description>Detach a hook from the event.</description>
  !@   </objectMethod>
  !@ </objectMethods>
  procedure :: attach     => eventHook{$interfaceType}Attach
  procedure :: isAttached => eventHook{$interfaceType}IsAttached
  procedure :: detach     => eventHook{$interfaceType}Detach
end type eventHook{$interfaceType}

abstract interface
 subroutine interface{$interfaceType}(self{scalar(@arguments) > 0 ? ",".join(",",@arguments) : ""})
  class(*), intent(inout) :: self
{$declarations}
 end subroutine interface{$interfaceType}
end interface
CODE
		    $code::location = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'});
		    my $attacher = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventHook{$interfaceType}Attach(self,object_,function_,openMPThreadBinding)
  !$ use OMP_Lib
  use Galacticus_Error, only : Galacticus_Error_Report
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)           :: self
  class    (*                        ), intent(in   ), target   :: object_
  integer                             , intent(in   ), optional :: openMPThreadBinding
  procedure(interface{$interfaceType})                          :: function_
  class    (hook                     )               , pointer  :: hook_
  integer                                                       :: i                  , openMPThreadBinding_

  !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{$location})
  !$ call OMP_Set_Lock(self%lock_)
  if (present(openMPThreadBinding)) then
     openMPThreadBinding_=openMPThreadBinding
  else
     openMPThreadBinding_=openMPThreadBindingNone
  end if
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
     hook_%object_             => object_
     hook_%function_           => function_
     hook_%openMPThreadBinding =  openMPThreadBinding_
     if (hook_%openMPThreadBinding == openMPThreadBindingAtLevel .or. hook_%openMPThreadBinding == openMPThreadBindingAllLevels) then
        hook_%openMPLevel=OMP_Get_Level()
        allocate(hook_%openMPThread(0:hook_%openMPLevel))
        do i=0,hook_%openMPLevel
           hook_%openMPThread(i)=OMP_Get_Ancestor_Thread_Num(i)
        end do
     end if
  end select
  self%count_=self%count_+1
  !$ call OMP_Unset_Lock(self%lock_)
  return
end subroutine eventHook{$interfaceType}Attach
CODE
		    my $attacherNode   =
		    {
			type       => "code",
			content    => $attacher,
			firstChild => undef(),
			source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
			line       => 1
		    };
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$attacherNode]);
		    my $detacher = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventHook{$interfaceType}Detach(self,object_,function_)
  use Galacticus_Error, only : Galacticus_Error_Report
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)          :: self
  class    (*                        ), intent(in   ), target  :: object_
  procedure(interface{$interfaceType})                         :: function_
  class    (hook                     )               , pointer :: hook_    , hookPrevious_

  !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{$location})
  !$ call OMP_Set_Lock(self%lock_)
  if (associated(self%first_)) then
     hookPrevious_ => null()
     hook_         => self%first_
     do while (associated(hook_))
        select type (hook_)
        type is (hook{$interfaceType})
           if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
             self%count_=self%count_-1
             if (associated(hookPrevious_)) then
               hookPrevious_%next   => hook_%next
             else
               self         %first_ => hook_%next
             end if
             deallocate(hook_)
             !$ call OMP_Unset_Lock(self%lock_)
             return
           end if
        end select
        hookPrevious_ => hook_
        hook_         => hook_%next
     end do
  end if
  !$ call OMP_Unset_Lock(self%lock_)
  call Galacticus_Error_Report('object/function not attached to this event'//{$location})
  return
end subroutine eventHook{$interfaceType}Detach
CODE
		    my $detacherNode   =
		    {
			type       => "code",
			content    => $detacher,
			firstChild => undef(),
			source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
			line       => 1
		    };
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$detacherNode]);
		    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},"hook".$code::interfaceType,"public");
		    my $isAttacher = fill_in_string(<<'CODE', PACKAGE => 'code');
logical function eventHook{$interfaceType}IsAttached(self,object_,function_)
  use Galacticus_Error, only : Galacticus_Error_Report
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)          :: self
  class    (*                        ), intent(in   ), target  :: object_
  procedure(interface{$interfaceType})                         :: function_
  class    (hook                     )               , pointer :: hook_

  !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{$location})
  !$ call OMP_Set_Lock(self%lock_)
  if (associated(self%first_)) then
     hook_ => self%first_
     do while (associated(hook_))
        select type (hook_)
        type is (hook{$interfaceType})
           if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
             eventHook{$interfaceType}IsAttached=.true.
             !$ call OMP_Unset_Lock(self%lock_)
             return
           end if
        end select
        hook_ => hook_%next
     end do
  end if
  eventHook{$interfaceType}IsAttached=.false.
  !$ call OMP_Unset_Lock(self%lock_)
  return
end function eventHook{$interfaceType}IsAttached
CODE
		    my $isAttacherNode   =
		    {
			type       => "code",
			content    => $isAttacher,
			firstChild => undef(),
			source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
			line       => 1
		    };
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$isAttacherNode]);
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
            # Build a function to initialize all event hooks.
            my $initializor = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksInitialize()
   implicit none
CODE
            foreach my $hook ( @hooks ) {
		$initializor .= "   call ".$hook->{'name'}."Event%initialize()\n";
            }
            $initializor .= fill_in_string(<<'CODE', PACKAGE => 'code');
   return
end subroutine eventsHooksInitialize
CODE
            my $initializorNode   =
            {
		type       => "code",
		content    => $initializor,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$initializorNode]);
            &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},"eventsHooksInitialize","public");
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
		    Events_Hooks     =>
		    {
			intrinsic => 0,
			all       => 1
		    },
                    Galacticus_Error =>
		    {
			intrinsic => 0,
			all       => 1
		    },
		    OMP_Lib         =>
		    {
			intrinsic => 0,
			all       => 1,
			openMP    => 1
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Insert required variables.
	    my @declarationsPotential =
		(
		 {
		     intrinsic  => "class"      ,
		     type       => "hook"       ,
		     variables  => [ "hook_"   ],
		     attributes => [ "pointer" ]
		 },
                 {
		     intrinsic  => "logical"            ,
		     variables  => [ "functionActive_" ]
		 },
                 {
		     intrinsic  => "integer"      ,
		     variables  => [ "ompLevel_" ]
		 }
		);
            my @declarations;
            foreach my $declaration ( @declarationsPotential ) {
		push(@declarations,$declaration)
                    unless ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node->{'parent'},$declaration->{'variables'}->[0]) );
            }
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
	    # Create the code.
	    $code::interfaceType = &interfaceTypeGet($node->{'directive'});
	    $code::callWith      = exists($node->{'directive'}->{'callWith'}) ? ",".$node->{'directive'}->{'callWith'} : "";
	    $code::eventName     = $node->{'directive'}->{'name'    };
            $code::location      = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'});
	    my $eventHookCode    = fill_in_string(<<'CODE', PACKAGE => 'code');
call {$eventName}Event%lock()
hook_ => {$eventName}Event%first()
do while (associated(hook_))
   select type (hook_)
   type is (hook{$interfaceType})
      select case (hook_%openMPThreadBinding)
      case (openMPThreadBindingNone)
         ! Not bound to any OpenMP thread, so always call.
         functionActive_=.true.
      case (openMPThreadBindingAtLevel)
         ! Binds at the OpenMP level - check levels match, and that this hooked object matches the OpenMP thread number across all levels.
         if (hook_%openMPLevel == OMP_Get_Level()) then
            functionActive_=.true.
            do ompLevel_=0,hook_%openMPLevel
               if (hook_%openMPThread(ompLevel_) /= OMP_Get_Ancestor_Thread_Num(ompLevel_)) then
                  functionActive_=.false.
                  exit
               end if
            end do
         else
            functionActive_=.false.
         end if
      case (openMPThreadBindingAllLevels)
         ! Binds at all levels at or above the level of the hooked object - check this condition is met, and that the hooked object matches the OpenMP thread number across all levels.
         if (hook_%openMPLevel <= OMP_Get_Level()) then
            functionActive_=.true.
            do ompLevel_=0,OMP_Get_Level()
               if (hook_%openMPThread(min(ompLevel_,hook_%openMPLevel)) /= OMP_Get_Ancestor_Thread_Num(ompLevel_)) then
                  functionActive_=.false.
                  exit
               end if
            end do
         else
            functionActive_=.false.
         end if
      case default
         functionActive_=.false.
         call Galacticus_Error_Report('unknown OpenMP binding'//{$location})
      end select
      if (functionActive_) call hook_%function_(hook_%object_{$callWith})
   end select
   hook_ => hook_%next
end do
call {$eventName}Event%unlock()
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
