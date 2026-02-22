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
		# Skip hooks that are duplicates.
		next
		    if ( exists($hook->{'isDuplicate'}) && $hook->{'isDuplicate'} eq "yes" );
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
			foreach my $declarator ( sort(keys(%Fortran::Utils::intrinsicDeclarations)) ) {
			    if ( my @matches = ( $processedLine =~ $Fortran::Utils::intrinsicDeclarations{$declarator}->{'regEx'} ) ) {
				push(@code::arguments,&Fortran::Utils::Extract_Variables($matches[$Fortran::Utils::intrinsicDeclarations{$declarator}->{'variables'}],keepQualifiers => 0));
				last;
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
			my @modules;
			if ( exists($hook->{'import'}->{'module'}->{'name'}) ) {
			    @modules = &List::ExtraUtils::as_array($hook->{'import'}->{'module'}                 );
			} else {
			    @modules = &List::ExtraUtils::hashList($hook->{'import'}->{'module'}, keyAs => "name");
			}
			foreach my $module ( @modules ) {
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
  !![
  <methods>
    <method method="attach"     description="Attach a hook to the event."                         />
    <method method="isAttached" description="Return true if the object is attached to this event."/>
    <method method="detach"     description="Detach a hook from the event."                       />
  </methods>
  !!]
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
    !!\{
    Attach an object to an event hook.
    !!\}
    use    :: Display           , only : displayMessage             , verbosityLevelInfo
    use    :: Error             , only : Error_Report
    use    :: ISO_Varying_String, only : varying_string             , var_str           , assignment(=), operator(//)
    use    :: String_Handling   , only : operator(//)
    !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
    implicit none
    class     (eventHook{$interfaceType}         ), intent(inout)                            :: self
    class     (*                                 ), intent(in   ), target                    :: object_
    type      (enumerationOpenMPThreadBindingType), intent(in   ), optional                  :: openMPThreadBinding
    character (len=*                             ), intent(in   ), optional                  :: label
    class     (dependency                        ), intent(in   ), optional   , dimension(:) :: dependencies
    procedure (interface{$interfaceType}         )                                           :: function_
    type      (hookList                          )               , allocatable, dimension(:) :: hooksTmp
    type      (hook{$interfaceType}              )                            , pointer      :: hook_
    type      (varying_string                    )                                           :: threadLabel      , message
    integer                                                                                  :: ompLevelEffective
    !$ integer                                                                               :: i
    !![
    <optionalArgument name="openMPThreadBinding" defaultsTo="openMPThreadBindingNone" />
    !!]

    ! Validate the thread binding model.
    if (self%isGlobal) then
       if (openMPThreadBinding_ /= openMPThreadBindingNone) call Error_Report("global event hooks permit only 'openMPThreadBindingNone'"         //{$location})
    else
       if (openMPThreadBinding_ == openMPThreadBindingNone) call Error_Report("threadprivate event hooks do not permit 'openMPThreadBindingNone'"//{$location})
    end if
    ! Check if atLevel attachment should be promoted.
    if (atLevelToAllLevels_ .and. openMPThreadBinding_ == openMPThreadBindingAtLevel) &
         openMPThreadBinding_=openMPThreadBindingAllLevels
    !$ if (self%isGlobal) call self%lock()
    ! Resize the array of hooks.
    if (allocated(self%hooks_)) then
       call move_alloc(self%hooks_,hooksTmp)
       allocate(self%hooks_(self%count_+1))
       self%hooks_(1:self%count_)=hooksTmp
       deallocate(hooksTmp)
    else
       allocate(self%hooks_(1))
    end if
    ! Create the new hook.
    allocate(hook_)
    hook_%object_             => object_
    hook_%function_           => function_
    hook_%openMPThreadBinding =  openMPThreadBinding_
    if (present(label)) then
       hook_%label=label
    else
       hook_%label=""
    end if
    !$omp atomic
    eventID            =eventID+1
    hook_      %eventID=eventID
    threadLabel        =""
    !$ threadLabel=" from thread "
    !$ ompLevelEffective=OMP_Get_Level()
    !$ if (futureThread_ /= -1) ompLevelEffective=ompLevelEffective+1
    !$ hook_%openMPLevel=ompLevelEffective
    !$ allocate(hook_%openMPThread(0:hook_%openMPLevel))
    !$ do i=0,hook_%openMPLevel
    !$    if (i == hook_%openMPLevel .and. futureThread_ /= -1) then
    !$      hook_%openMPThread(i)=futureThread_
    !$    else
    !$      hook_%openMPThread(i)=OMP_Get_Ancestor_Thread_Num(i)
    !$    end if
    !$    if (i > 0) threadLabel=threadLabel//" -> "
    !$    threadLabel=threadLabel//hook_%openMPThread(i)
    !$ end do
    ! Insert the hook into the list.
    self%hooks_(self%count_+1)%hook_ => hook_
    ! Increment the count of hooks into this event and resolve dependencies.
    self%count_=self%count_+1
    call self%resolveDependencies(hook_,dependencies)
    ! Report
    message=var_str("attaching '")//trim(hook_%label)//"' ["//hook_%eventID//"] to event"//trim(self%label)//threadLabel//" [count="//self%count_//"]"
    call displayMessage(message,verbosityLevelInfo)
    !$ if (self%isGlobal) call self%unlock()
    return
  end subroutine eventHook{$interfaceType}Attach
CODE
		    my $attacherTree  = &Galacticus::Build::SourceTree::ParseCode($attacher,"null()");
		    my @attacherNodes = &Galacticus::Build::SourceTree::Children($attacherTree);
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@attacherNodes);
		    my $detacher = fill_in_string(<<'CODE', PACKAGE => 'code');
  subroutine eventHook{$interfaceType}Detach(self,object_,function_)
    !!\{
    Attach an object to an event hook.
    !!\}
    use    :: Display           , only : displayMessage             , verbosityLevelInfo
    use    :: Error             , only : Error_Report
    use    :: ISO_Varying_String, only : varying_string             , var_str           , assignment(=), operator(//)
    use    :: String_Handling   , only : operator(//)
    !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
    implicit none
    class    (eventHook{$interfaceType}), intent(inout)               :: self
    class    (*                        ), intent(in   ), target       :: object_
    procedure(                         )                              :: function_
    type     (hookList                 ), allocatable  , dimension(:) :: hooksTmp
    type     (varying_string           )                              :: threadLabel, message
    integer                                                           :: i          , j
    
    !$ if (self%isGlobal) call self%lock()
    if (allocated(self%hooks_)) then
       do i=1,self%count_
          select type (hook_ => self%hooks_(i)%hook_)
          type is (hook{$interfaceType})
             if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
                ! Report
                threadLabel   =""
                !$ threadLabel=" from thread "
                !$ do j=0,OMP_Get_Level()
                !$    if (j > 0) threadLabel=threadLabel//" -> "
                !$    threadLabel=threadLabel//OMP_Get_Ancestor_Thread_Num(j)
                !$ end do
                message=var_str("detaching '")//trim(self%hooks_(i)%hook_%label)//"' ["//self%hooks_(i)%hook_%eventID//"] from event"//trim(self%label)//threadLabel//" [count="//self%count_//"]"
                call displayMessage(message,verbosityLevelInfo)
                deallocate(self%hooks_(i)%hook_)
                if (self%count_ > 1) then
                   call move_alloc(self%hooks_,hooksTmp)
                   allocate(self%hooks_(self%count_-1))
                   if (i >           1) self%hooks_(1:          i-1)=hooksTmp(1  :          i-1)
                   if (i < self%count_) self%hooks_(i:self%count_-1)=hooksTmp(i+1:self%count_  )
                   deallocate(hooksTmp)
                else
                   deallocate(self%hooks_)
                end if
                self%count_=self%count_-1
                !$ if (self%isGlobal) call self%unlock()
                return
             end if
          end select
       end do
    end if
    call Error_Report('object/function not attached to this event'//{$location})
    !$ if (self%isGlobal) call self%unlock()
    return
  end subroutine eventHook{$interfaceType}Detach
CODE
		    my $detacherTree  = &Galacticus::Build::SourceTree::ParseCode($detacher,"null()");
		    my @detacherNodes = &Galacticus::Build::SourceTree::Children($detacherTree);
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@detacherNodes);
		    my $isAttacher = fill_in_string(<<'CODE', PACKAGE => 'code');
  logical function eventHook{$interfaceType}IsAttached(self,object_,function_)
    !!\{
    Return true if an object is attached to an event hook.
    !!\}
    use :: Error, only : Error_Report
    implicit none
    class    (eventHook{$interfaceType}), intent(inout)          :: self
    class    (*                        ), intent(in   ), target  :: object_
    procedure(                         )                         :: function_
    integer                                                      :: i
    
    !$ if (self%isGlobal) call self%lock(writeLock=.false.)
    if (allocated(self%hooks_)) then
       do i=1,self%count_
          select type (hook_ => self%hooks_(i)%hook_)
          type is (hook{$interfaceType})
             if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
                eventHook{$interfaceType}IsAttached=.true.
                !$ if (self%isGlobal) call self%unlock(writeLock=.false.)
                return
             end if
          end select
       end do
    end if
    eventHook{$interfaceType}IsAttached=.false.
    !$ if (self%isGlobal) call self%unlock(writeLock=.false.)
    return
  end function eventHook{$interfaceType}IsAttached
CODE
		    my $isAttacherTree  = &Galacticus::Build::SourceTree::ParseCode($isAttacher,"null()");
		    my @isAttacherNodes = &Galacticus::Build::SourceTree::Children($isAttacherTree);
		    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},\@isAttacherNodes);
		    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},"hook".$code::interfaceType,"public");
	        }
		$hookObject .= "type(eventHook".$code::interfaceType."), public  :: ".$hook->{'name'}."Event                 , ".$hook->{'name'}."EventGlobal\n";
		$hookObject .= "type(eventHook".$code::interfaceType.")          :: ".$hook->{'name'}."Event_\n";
		$hookObject .= "type(eventHookList                    ), pointer :: ".$hook->{'name'}."EventBackups => null()\n";
		$hookObject .= "!\$omp threadprivate (".$hook->{'name'}."Event,".$hook->{'name'}."EventBackups)\n";
		my $hookTree = &Galacticus::Build::SourceTree::ParseCode($hookObject,"null()");
		my @hookNodes = &Galacticus::Build::SourceTree::Children($hookTree);
		&Galacticus::Build::SourceTree::InsertAfterNode($node,\@hookNodes);
	    }
            # Build a function to perform copy out of the current event lists before entering a new OpenMP parallel region.
            my $copyOut = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksFilterCopyOut_()
   implicit none
   call copyLock%set()
CODE
            foreach my $hook ( @hooks ) {
		$copyOut .= "   ".$hook->{'name'}."Event_=".$hook->{'name'}."Event\n";
            }
            $copyOut .= fill_in_string(<<'CODE', PACKAGE => 'code');
   return
end subroutine eventsHooksFilterCopyOut_
CODE
            my $copyOutNode   =
            {
		type       => "code",
		content    => $copyOut,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$copyOutNode]);
            # Build a function to perform copy in of the current event lists on entering a new OpenMP parallel region.
            my $copyIn = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksFilterCopyIn_()
   use    :: Display           , only : displayMessage             , verbosityLevelInfo
   use    :: ISO_Varying_String, only : var_str                    , operator(//)      , varying_string, assignment(=)
   use    :: String_Handling   , only : operator(//)
   !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
   implicit none
   type   (eventHookList ), pointer :: eventHookBackup
   type   (varying_string)          :: threadLabel    , message
   integer                          :: i

   threadLabel=""
   !$ threadLabel=" from thread "
   !$ do i=0,OMP_Get_Level()
   !$    if (i > 0) threadLabel=threadLabel//" -> "
   !$    threadLabel=threadLabel//OMP_Get_Ancestor_Thread_Num(i)
   !$ end do
CODE
            foreach my $hook ( @hooks ) {
		$code::name = $hook->{'name'};
		$copyIn .= fill_in_string(<<'CODE', PACKAGE => 'code');
   allocate(eventHookBackup)
   {$name}Event               =  {$name}Event_
   eventHookBackup%eventHook_ =  {$name}Event_
   if (associated({$name}EventBackups)) eventHookBackup%next => {$name}EventBackups
   {$name}EventBackups        => eventHookBackup
   message=var_str("{$name}: storing ")//eventHookBackup%eventHook_%count_//" hooks"//threadLabel
   call displayMessage(message,verbosityLevelInfo)
   nullify(eventHookBackup)
CODE
            }
            $copyIn .= fill_in_string(<<'CODE', PACKAGE => 'code');
   return
end subroutine eventsHooksFilterCopyIn_
CODE
            my $copyInNode   =
            {
		type       => "code",
		content    => $copyIn,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$copyInNode]);
            # Build a function to perform restore of the current event lists before leaving a OpenMP parallel region.
            my $restore = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksFilterRestore_()
   use    :: Display           , only : displayMessage             , verbosityLevelInfo
   use    :: Error             , only : Error_Report
   use    :: ISO_Varying_String, only : var_str                    , operator(//)      , varying_string, assignment(=)
   use    :: String_Handling   , only : operator(//)
   !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
   implicit none
   type   (eventHookList ), pointer :: eventHookBackup
   type   (varying_string)          :: threadLabel    , message
   integer                          :: i

   threadLabel=""
   !$ threadLabel=" from thread "
   !$ do i=0,OMP_Get_Level()
   !$    if (i > 0) threadLabel=threadLabel//" -> "
   !$    threadLabel=threadLabel//OMP_Get_Ancestor_Thread_Num(i)
   !$ end do
CODE
            foreach my $hook ( @hooks ) {
		$code::name          = $hook->{'name'};
		$code::interfaceType = &interfaceTypeGet($hook);
		$code::location      = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'});
		$restore .= fill_in_string(<<'CODE', PACKAGE => 'code');
   eventHookBackup     => {$name}EventBackups
   {$name}EventBackups => {$name}EventBackups%next
   select type (eventHook_ => eventHookBackup%eventHook_)
   type is (eventHook{$interfaceType})
      if (allocated({$name}Event%hooks_)) deallocate({$name}Event%hooks_)
      {$name}Event        =  eventHook_
   class default
      call Error_Report('eventHook has incorrect class'//{$location})
   end select
   message=var_str("{$name}: restoring ")//eventHookBackup%eventHook_%count_//" hooks"//threadLabel
   call displayMessage(message,verbosityLevelInfo)
   deallocate(eventHookBackup)
   nullify   (eventHookBackup)
CODE
            }
            $restore .= fill_in_string(<<'CODE', PACKAGE => 'code');
   return
end subroutine eventsHooksFilterRestore_
CODE
            my $restoreNode   =
            {
		type       => "code",
		content    => $restore,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$restoreNode]);
	    # Build a function to finalize copy of the current event lists on entering a new OpenMP parallel region.
            my $copyDone = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksFilterCopyDone_()
   implicit none
   call copyLock%unset()
   return
end subroutine eventsHooksFilterCopyDone_
CODE
            my $copyDoneNode   =
            {
		type       => "code",
		content    => $copyDone,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$copyDoneNode]);
            # Build a function to filter the list of hooks on entering a new OpenMP parallel region.
            my $filter = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksFilterFunction_()
   implicit none
CODE
            foreach my $hook ( @hooks ) {
		$filter .= "   call ".$hook->{'name'}."Event%filter()\n";
            }
            $filter .= fill_in_string(<<'CODE', PACKAGE => 'code');
   return
end subroutine eventsHooksFilterFunction_
CODE
            my $filterNode   =
            {
		type       => "code",
		content    => $filter,
		firstChild => undef(),
                source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
                line       => 1
	    };
            &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$filterNode]);
	    # Build a function to initialize all event hooks.
	    my $initializor = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine eventsHooksInitialize()
   use :: Events_Filters, only : eventsHooksFilterFunction, eventsHooksFilterCopyOut, eventsHooksFilterCopyIn, eventsHooksFilterCopyDone, &
        &                        eventsHooksFilterRestore
   implicit none

   eventsHooksFilterFunction => eventsHooksFilterFunction_
   eventsHooksFilterCopyOut  => eventsHooksFilterCopyOut_
   eventsHooksFilterCopyIn   => eventsHooksFilterCopyIn_
   eventsHooksFilterCopyDone => eventsHooksFilterCopyDone_
   eventsHooksFilterRestore  => eventsHooksFilterRestore_
   copyLock=ompLock()
CODE
	    foreach my $hook ( @hooks ) {
		$initializor .= "   ".$hook->{'name'}."EventGlobal%isGlobal=.true.\n";
		$initializor .= "   ".$hook->{'name'}."EventGlobal%lock_   =ompReadWriteLock()\n";
		$initializor .= "   ".$hook->{'name'}."Event      %label   ='".$hook->{'name'}."'\n";
		$initializor .= "   ".$hook->{'name'}."EventGlobal%label   ='".$hook->{'name'}." (global)'\n";
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
    use :: Output_HDF5       , only : outputFile
    use :: IO_HDF5           , only : hdf5Object
    use :: HDF5_Access       , only : hdf5Access
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
    metaDataGroup=outputFile%openGroup('metaData','Galacticus meta data.'           )
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
                    Error =>
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
	    my @declarations =
		(
                 {
                     intrinsic     => "integer",
		     variables     => [ $node->{'directive'}->{'name'}."Iterator" ]
		 }
		);
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@declarations);
	    # Create the code.
	    $code::interfaceType = &interfaceTypeGet($node->{'directive'});
	    $code::callWith      = exists($node->{'directive'}->{'callWith'}) ? ",".$node->{'directive'}->{'callWith'} : "";
	    $code::eventName     = $node->{'directive'}->{'name'};
	    my $eventHookCode    = fill_in_string(<<'CODE', PACKAGE => 'code');
if ({$eventName}EventGlobal%count() > 0) then
  do {$eventName}Iterator=1,{$eventName}EventGlobal%count()
     select type (hook_ => {$eventName}EventGlobal%hooks_({$eventName}Iterator)%hook_)
     type is (hook{$interfaceType})
       call hook_%function_(hook_%object_{$callWith})
     end select
  end do
end if
if ({$eventName}Event%count() > 0) then
  do {$eventName}Iterator=1,{$eventName}Event%count()
     select type (hook_ => {$eventName}Event%hooks_({$eventName}Iterator)%hook_)
     type is (hook{$interfaceType})
       call hook_%function_(hook_%object_{$callWith})
     end select
  end do
end if
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
	# Handle OpenMP parallel sections by adding copyin of our hooks, followed by per-thread filtering.
	if ( $node->{'type'} eq "openMP" && $node->{'name'} eq "parallel" && ! $node->{'isCloser'} && ! exists($node->{'eventFilterInserted'}) ) {
	    $node->{'eventFilterInserted'} =  1;
	    # Find all hook directives.
	    my @hooks = map {&Galacticus::Build::Directives::Extract_Directives($_,'eventHook')} &List::ExtraUtils::as_array($directiveLocations->{'eventHook'}->{'file'});
	    # Insert the required module uses.
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse => 
		{
		    Events_Filters =>
		    {
			intrinsic => 0,
			only => {
			    eventsHooksFilterFunction => 1,
			    eventsHooksFilterCopyOut  => 1,
			    eventsHooksFilterCopyIn   => 1,
			    eventsHooksFilterCopyDone => 1
			}
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Insert a call to our filter function.
	    my $copyOutNode =
	    {
		type       => "code",
		content    => "call eventsHooksFilterCopyOut()\n",
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertBeforeNode($node,[$copyOutNode]);
	    # Insert a call to our filter function.
	    my $filterCode = fill_in_string(<<'CODE', PACKAGE => 'code');
call eventsHooksFilterCopyIn()
!$omp barrier
!$omp single
call eventsHooksFilterCopyDone()
!$omp end single
call eventsHooksFilterFunction()
CODE
	    my $filterNode =
	    {
		type       => "code",
		content    => $filterCode,
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$filterNode]);

	}
	# Handle OpenMP end parallel sections by adding restore of our hooks.
	if ( $node->{'type'} eq "openMP" && $node->{'name'} eq "parallel" && $node->{'isCloser'} && ! exists($node->{'eventFilterInserted'}) ) {
	    $node->{'eventFilterInserted'} =  1;
	    # Find all hook directives.
	    my @hooks = map {&Galacticus::Build::Directives::Extract_Directives($_,'eventHook')} &List::ExtraUtils::as_array($directiveLocations->{'eventHook'}->{'file'});
	    # Insert the required module uses.
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse => 
		{
		    Events_Filters =>
		    {
			intrinsic => 0,
			only => {
			    eventsHooksFilterRestore => 1
			}
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Insert a call to our restore function.
	    my $restoreNode =
	    {
		type       => "code",
		content    => "call eventsHooksFilterRestore()\n",
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::EventHooks::Process_EventHooks()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertBeforeNode($node,[$restoreNode]);
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
