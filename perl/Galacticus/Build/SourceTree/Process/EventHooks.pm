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
use List::Util 'max';
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
		    # Look for symbols to import.
		    $code::imports = "";
		    if ( exists($hook->{'import'}) ) {
			my $usesNode =
			{
			    type      => "moduleUse"
			};
			my @imports;
			foreach my $module ( &List::ExtraUtils::as_array($hook->{'import'}->{'module'}) ) {
			    my @symbolNames = split(/\s*,\s*/,$module->{'symbols'});
			    push(@imports,@symbolNames);
			    my %symbols;
			    foreach my $symbol ( @symbolNames ) {
				$symbols{$symbol} = 1;
			    }
			    $usesNode->{'moduleUse'}->{$module->{'name'}} =
			    {
				intrinsic => 0,
				only      => \%symbols
			    };
			}
			$code::imports = "  import ".join(", ",@imports)." \n";
			# Insert the modules.
			&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
		    }
		    # Build the required types and functions.
		    $hookObject = fill_in_string(<<'CODE', PACKAGE => 'code');

type, extends(hook) :: hook{$interfaceType}
   procedure(interface{$interfaceType}), pointer, nopass :: function_ => null()
end type hook{$interfaceType}
 
type, extends(eventHook) :: eventHook{$interfaceType}
  private
 contains
  !# <methods>
  !#   <method method="attach"     description="Attach a hook to the event."                         />
  !#   <method method="isAttached" description="Return true if the object is attached to this event."/>
  !#   <method method="detach"     description="Detach a hook from the event."                       />
  !# </methods>
  procedure :: attach     => eventHook{$interfaceType}Attach
  procedure :: isAttached => eventHook{$interfaceType}IsAttached
  procedure :: detach     => eventHook{$interfaceType}Detach
end type eventHook{$interfaceType}

abstract interface
 subroutine interface{$interfaceType}(self{scalar(@arguments) > 0 ? ",".join(",",@arguments) : ""})
{$imports}
  class(*), intent(inout) :: self
{$declarations}
 end subroutine interface{$interfaceType}
end interface
CODE
		    $code::location = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'});
		    my $attacher = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventHook{$interfaceType}Attach(self,object_,function_,openMPThreadBinding,label,dependencies)
  use    :: Galacticus_Error  , only : Galacticus_Error_Report
  use    :: ISO_Varying_String, only : assignment(=)
  !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)                         :: self
  class    (*                        ), intent(in   ), target                 :: object_
  integer                             , intent(in   ), optional               :: openMPThreadBinding
  character(len=*                    ), intent(in   ), optional               :: label
  class    (dependency               ), intent(in   ), optional, dimension(:) :: dependencies
  procedure(interface{$interfaceType})                                        :: function_
  class    (hook                     )               , pointer                :: hook_
  integer                                                                     :: i                  , openMPThreadBinding_

  !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{$location})
  call self%lock()
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
     if (present(label)) then
        hook_%label=label
     else
        hook_%label=""
     end if
     !$ if (hook_%openMPThreadBinding == openMPThreadBindingAtLevel .or. hook_%openMPThreadBinding == openMPThreadBindingAllLevels) then
     !$    hook_%openMPLevel=OMP_Get_Level()
     !$    allocate(hook_%openMPThread(0:hook_%openMPLevel))
     !$    do i=0,hook_%openMPLevel
     !$       hook_%openMPThread(i)=OMP_Get_Ancestor_Thread_Num(i)
     !$    end do
     !$ end if
  end select
  ! Increment the count of hooks into this event and resolve dependencies.
  self%count_=self%count_+1
  call self%resolveDependencies(hook_,dependencies)
  call self%unlock             (                  )
  return
end subroutine eventHook{$interfaceType}Attach
CODE
		    my $attacherTree  = &Galacticus::Build::SourceTree::ParseCode($attacher,"null()");
		    my @attacherNodes = &Galacticus::Build::SourceTree::Children($attacherTree);
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@attacherNodes);
		    my $detacher = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventHook{$interfaceType}Detach(self,object_,function_)
  use Galacticus_Error, only : Galacticus_Error_Report
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)          :: self
  class    (*                        ), intent(in   ), target  :: object_
  procedure(interface{$interfaceType})                         :: function_
  class    (hook                     )               , pointer :: hook_    , hookPrevious_

  !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{$location})
  call self%lock()
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
             call self%unlock()
             return
           end if
        end select
        hookPrevious_ => hook_
        hook_         => hook_%next
     end do
  end if
  call self%unlock()
  call Galacticus_Error_Report('object/function not attached to this event'//{$location})
  return
end subroutine eventHook{$interfaceType}Detach
CODE
		    my $detacherTree  = &Galacticus::Build::SourceTree::ParseCode($detacher,"null()");
		    my @detacherNodes = &Galacticus::Build::SourceTree::Children($detacherTree);
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@detacherNodes);
		    my $isAttacher = fill_in_string(<<'CODE', PACKAGE => 'code');
logical function eventHook{$interfaceType}IsAttached(self,object_,function_)
  use Galacticus_Error, only : Galacticus_Error_Report
  implicit none
  class    (eventHook{$interfaceType}), intent(inout)          :: self
  class    (*                        ), intent(in   ), target  :: object_
  procedure(interface{$interfaceType})                         :: function_
  class    (hook                     )               , pointer :: hook_

  !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{$location})
  call self%lock(writeLock=.false.)
  if (associated(self%first_)) then
     hook_ => self%first_
     do while (associated(hook_))
        select type (hook_)
        type is (hook{$interfaceType})
           if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
             eventHook{$interfaceType}IsAttached=.true.
             call self%unlock(writeLock=.false.)
             return
           end if
        end select
        hook_ => hook_%next
     end do
  end if
  eventHook{$interfaceType}IsAttached=.false.
  call self%unlock(writeLock=.false.)
  return
end function eventHook{$interfaceType}IsAttached
CODE
		    my $isAttacherTree  = &Galacticus::Build::SourceTree::ParseCode($isAttacher,"null()");
		    my @isAttacherNodes = &Galacticus::Build::SourceTree::Children($isAttacherTree);
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@isAttacherNodes);
		    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},"hook".$code::interfaceType,"public");
	        }
		$hookObject .= "type(eventHook".$code::interfaceType."), public :: ".$hook->{'name'}."Event\n";
		my $hookTree = &Galacticus::Build::SourceTree::ParseCode($hookObject,"null()");
		my @hookNodes = &Galacticus::Build::SourceTree::Children($hookTree);
		&Galacticus::Build::SourceTree::InsertAfterNode($node,\@hookNodes);
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
            # Build a function to output meta-data on event hook lock wait times.
            $code::hookCount             = scalar(                           @hooks);
            $code::hookNameLengthMaximum = max   (map {length($_->{'name'})} @hooks);
            my $waitTimeWriter           = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksWaitTimes()
#ifdef OMPPROFILE
    use :: Galacticus_HDF5   , only : galacticusOutputFile
    use :: IO_HDF5           , only : hdf5Access          , hdf5Object
    use :: ISO_Varying_String, only : varying_string      , var_str
#endif
    implicit none
#ifdef OMPPROFILE
    type            (hdf5Object                  )                          :: waitTimeGroup         , waitTimeDataset        , metaDataGroup
    character       (len={$hookNameLengthMaximum}), dimension({$hookCount}) :: eventHookNames
    double precision                              , dimension({$hookCount}) :: eventHookReadWaitTimes, eventHookWriteWaitTimes

CODE
my $i = 0;
foreach my $hook ( @hooks ) {
    ++$i;
    $waitTimeWriter .= "   eventHookNames         (".$i.")='".$hook->{'name'}."'\n";
    $waitTimeWriter .= "   eventHookReadWaitTimes (".$i.")=" .$hook->{'name'}."Event%waitTimeRead\n";
    $waitTimeWriter .= "   eventHookWriteWaitTimes(".$i.")=" .$hook->{'name'}."Event%waitTimeWrite\n";
}
$waitTimeWriter .= fill_in_string(<<'CODE', PACKAGE => 'code');
    ! Open output group.
    !$ call hdf5Access%set()
    metaDataGroup=galacticusOutputFile%openGroup('metaData','Galacticus meta data.'           )
    waitTimeGroup=metaDataGroup       %openGroup('openMP'  ,'Meta-data on OpenMP performance.')
    ! Write wait time data.
    call waitTimeGroup%writeDataset(eventHookNames         ,"eventHookNames"         ,"Names of event hooks"                                                              )
    call waitTimeGroup%writeDataset(eventHookReadWaitTimes ,"eventHookReadWaitTimes" ,"Total time spent waiting to read-lock event hooks" ,datasetReturned=waitTimeDataset)
    call waitTimeDataset%writeAttribute(1.0d0,"unitsInSI")
    call waitTimeDataset%close()
    call waitTimeGroup%writeDataset(eventHookWriteWaitTimes,"eventHookWriteWaitTimes","Total time spent waiting to write-lock event hooks",datasetReturned=waitTimeDataset)
    call waitTimeDataset%writeAttribute(1.0d0,"unitsInSI")
    call waitTimeDataset%close()
    ! Close output groups.
    call waitTimeGroup%close()
    call metaDataGroup%close()
    !$ call hdf5Access%unset()
#endif
   return
end subroutine eventsHooksWaitTimes
CODE
my $waitTimeWriterNode   =
            {
		type       => "code",
		content    => $waitTimeWriter,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$waitTimeWriterNode]);
            &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},"eventsHooksWaitTimes","public");




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
                    OMP_Lib          =>
		    {
			intrinsic => 0,
			openMP    => 1,
			all       => 1
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
		     variables  => [ "ompLevel_", "ompLevelCurrent_" ]
		 },
                 {
		     intrinsic  => "logical"      ,
		     variables  => [ "ompAncestorGot_" ]
		 },
                 {
                     intrinsic  => "integer"                        ,
                     attributes => [ "dimension(:)", "allocatable" ],
		     variables  => [ "ompAncestorThreadNum_" ]
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
call {$eventName}Event%lock(writeLock=.false.)
!$ ompAncestorGot_=.false.
!$ ompLevelCurrent_=-1
hook_ => {$eventName}Event%first()
do while (associated(hook_))
   select type (hook_)
   type is (hook{$interfaceType})
      select case (hook_%openMPThreadBinding)
      case (openMPThreadBindingAtLevel,openMPThreadBindingAllLevels)
         !$ if (.not.ompAncestorGot_) then
         !$    ompLevelCurrent_=OMP_Get_Level()
         !$    allocate(ompAncestorThreadNum_(0:ompLevelCurrent_))
         !$    do ompLevel_=0,ompLevelCurrent_
         !$       ompAncestorThreadNum_(ompLevel_)=OMP_Get_Ancestor_Thread_Num(ompLevel_)
         !$    end do
         !$    ompAncestorGot_=.true.
         !$ end if
      end select
      !$ select case (hook_%openMPThreadBinding)
      !$ case (openMPThreadBindingNone)
         ! Not bound to any OpenMP thread, so always call.
         functionActive_=.true.
      !$ case (openMPThreadBindingAtLevel)
      !$    ! Binds at the OpenMP level - check levels match, and that this hooked object matches the OpenMP thread number across all levels.
      !$    if (hook_%openMPLevel == ompLevelCurrent_) then
      !$       functionActive_=.true.
      !$       do ompLevel_=0,hook_%openMPLevel
      !$          if (hook_%openMPThread(ompLevel_) /= ompAncestorThreadNum_(ompLevel_)) then
      !$             functionActive_=.false.
      !$             exit
      !$          end if
      !$       end do
      !$    else
      !$       functionActive_=.false.
      !$    end if
      !$ case (openMPThreadBindingAllLevels)
      !$    ! Binds at all levels at or above the level of the hooked object - check this condition is met, and that the hooked object matches the OpenMP thread number across all levels.
      !$    if (hook_%openMPLevel <= ompLevelCurrent_) then
      !$       functionActive_=.true.
      !$       do ompLevel_=0,ompLevelCurrent_
      !$          if (hook_%openMPThread(min(ompLevel_,hook_%openMPLevel)) /= ompAncestorThreadNum_(ompLevel_)) then
      !$             functionActive_=.false.
      !$             exit
      !$          end if
      !$       end do
      !$    else
      !$       functionActive_=.false.
      !$    end if
      !$ case default
      !$    functionActive_=.false.
      !$    call Galacticus_Error_Report('unknown OpenMP binding'//{$location})
      !$ end select
      if (functionActive_) call hook_%function_(hook_%object_{$callWith})
   end select
   hook_ => hook_%next
end do
call {$eventName}Event%unlock(writeLock=.false.)
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
