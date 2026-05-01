# Processes `eventHookManager`, `eventHook`, and `openMP parallel` nodes
# to synthesize the full event-hooks infrastructure: the hook / eventHook
# types, attach / detach / isAttached methods, filter copy-in/copy-out/
# restore subroutines, initializer, wait-time writer, per-call-site
# dispatch blocks, and OpenMP-parallel wrappers.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/EventHooks.pm

import hashlib
import io
import os
import re
import xml.etree.ElementTree as ET


from Galacticus.Build.FortranUtils                           import get_fortran_line
from Fortran.Utils                                 import (
    INTRINSIC_DECLARATIONS, extract_variables,
)
from List.ExtraUtils                               import as_array
from XML.Utils                                     import xml_to_dict
from Galacticus.Build.Directives                   import extract_directives
from Galacticus.Build.SourceTree                   import (
    walk_tree, parse_code, children, set_visibility,
    insert_after_node, insert_before_node, insert_post_contains,
)
from Galacticus.Build.SourceTree.Process           import register_process, process_tree
from Galacticus.Build.SourceTree.Parse.Declarations import add_declarations
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

_DIRECTIVE_LOCATIONS = None


def _load_directive_locations():
    global _DIRECTIVE_LOCATIONS
    if _DIRECTIVE_LOCATIONS is not None:
        return _DIRECTIVE_LOCATIONS
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError("process_event_hooks: BUILDPATH is not set")
    path = os.path.join(build_path, 'directiveLocations.xml')
    _DIRECTIVE_LOCATIONS = xml_to_dict(ET.parse(path).getroot())
    return _DIRECTIVE_LOCATIONS


def _interface_type_get(hook):
    """Return the md5_hex of the hook name when an `<interface>` is declared,
    otherwise the string "Unspecified".

    Mirrors interfaceTypeGet() at EventHooks.pm:722-731.
    """
    if isinstance(hook, dict) and 'interface' in hook:
        return hashlib.md5(hook['name'].encode('utf-8')).hexdigest()
    return 'Unspecified'


def _collect_hooks(directive_locations):
    """Return every `<eventHook>` directive across every file listed under
    directiveLocations/eventHook/file, skipping duplicates.
    """
    block = directive_locations.get('eventHook') or {}
    files = list(as_array(block.get('file')))
    hooks = []
    for f in files:
        for h in extract_directives(f, 'eventHook'):
            if h.get('isDuplicate') == 'yes':
                continue
            hooks.append(h)
    return hooks


def _substitute(template, subs):
    """Replace every `{$key}` occurrence in `template` with `subs[key]`.

    A simple targeted substitutor: keeps Fortran `{…}` sequences and LaTeX
    `!!\\{` / `!!\\}` blocks untouched.  Matches the Perl
    `fill_in_string(<<'CODE', PACKAGE => 'code')` idiom without pulling in a
    templating engine.
    """
    out = template
    for key, value in subs.items():
        out = out.replace('{$' + key + '}', str(value))
    return out


def _code_node(content, source='EventHooks', line=1):
    return {
        'type':       'code',
        'content':    content,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     source,
        'line':       line,
    }


# ---------------------------------------------------------------------------
# Call-site dispatch block (`<eventHook name=… />` at a regular call site)
# ---------------------------------------------------------------------------

_CALL_SITE_TEMPLATE = """if ({$eventName}EventGlobal%count() > 0) then
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
"""


def _process_event_hook_call_site(node):
    """Handle an `<eventHook name=…/>` directive at a call site.

    Mirrors EventHooks.pm:562-631.
    """
    directive = node.setdefault('directive', {})
    directive['processed'] = True

    add_uses(node['parent'], {
        'moduleUse': {
            'Events_Hooks': {'intrinsic': False, 'all': True},
            'Error':        {'intrinsic': False, 'all': True},
            'OMP_Lib':      {'intrinsic': False, 'openMP': True, 'all': True},
        },
        'moduleOrder': ['Events_Hooks', 'Error', 'OMP_Lib'],
    })

    add_declarations(node['parent'], [{
        'intrinsic':     'integer',
        'type':          None,
        'openMP':        False,
        'attributes':    [],
        'variables':     [directive['name'] + 'Iterator'],
        'variableNames': [directive['name'] + 'Iterator'],
    }])

    code = _substitute(_CALL_SITE_TEMPLATE, {
        'eventName':     directive['name'],
        'interfaceType': _interface_type_get(directive),
        'callWith':      ',' + directive['callWith'] if 'callWith' in directive else '',
    })
    insert_after_node(node, [_code_node(
        code,
        source='Galacticus.Build.SourceTree.Process.EventHooks'
               '.process_event_hooks()',
    )])


# ---------------------------------------------------------------------------
# OpenMP-parallel wrappers (copy-out / copy-in / restore around !$omp parallel)
# ---------------------------------------------------------------------------

_OMP_FILTER_AFTER = """call eventsHooksFilterCopyIn()
!$omp barrier
!$omp single
call eventsHooksFilterCopyDone()
!$omp end single
call eventsHooksFilterFunction()
"""


def _process_openmp_opener(node):
    """Wrap a non-closer `!$omp parallel` with copy-out + copy-in filter calls.

    Mirrors EventHooks.pm:633-685.
    """
    node['eventFilterInserted'] = 1
    add_uses(node['parent'], {
        'moduleUse': {
            'Events_Filters': {
                'intrinsic': False,
                'only': {
                    'eventsHooksFilterFunction': True,
                    'eventsHooksFilterCopyOut':  True,
                    'eventsHooksFilterCopyIn':   True,
                    'eventsHooksFilterCopyDone': True,
                },
            },
        },
        'moduleOrder': ['Events_Filters'],
    })
    insert_before_node(node, [_code_node(
        "call eventsHooksFilterCopyOut()\n",
        source='Galacticus.Build.SourceTree.Process.EventHooks'
               '.process_event_hooks()',
    )])
    insert_after_node(node, [_code_node(
        _OMP_FILTER_AFTER,
        source='Galacticus.Build.SourceTree.Process.EventHooks'
               '.process_event_hooks()',
    )])


def _process_openmp_closer(node):
    """Inject a `call eventsHooksFilterRestore()` before `!$omp end parallel`.

    Mirrors EventHooks.pm:687-717.
    """
    node['eventFilterInserted'] = 1
    add_uses(node['parent'], {
        'moduleUse': {
            'Events_Filters': {
                'intrinsic': False,
                'only': {'eventsHooksFilterRestore': True},
            },
        },
        'moduleOrder': ['Events_Filters'],
    })
    insert_before_node(node, [_code_node(
        "call eventsHooksFilterRestore()\n",
        source='Galacticus.Build.SourceTree.Process.EventHooks'
               '.process_event_hooks()',
    )])


# ---------------------------------------------------------------------------
# eventHookManager — the big one.  Per hook with a declared <interface>, we
# synthesize a dedicated `hook{InterfaceType}` + `eventHook{InterfaceType}`
# pair with attach/detach/isAttached.  Then, across ALL hooks, we synthesize
# the module-level filter copy-out / copy-in / restore / copy-done subroutines,
# the initializer, and the (OMPPROFILE-gated) wait-time writer.
# Mirrors EventHooks.pm:33-561.
# ---------------------------------------------------------------------------

_HOOK_TYPE_TEMPLATE = """
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
 subroutine interface{$interfaceType}(self{$argumentList})
{$imports}
  class(*), intent(inout) :: self
{$declarations}
 end subroutine interface{$interfaceType}
end interface
"""


_ATTACH_TEMPLATE = """  subroutine eventHook{$interfaceType}Attach(self,object_,function_,openMPThreadBinding,label,dependencies)
    !!{
    Attach an object to an event hook.
    !!}
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
"""


_DETACH_TEMPLATE = """  subroutine eventHook{$interfaceType}Detach(self,object_,function_)
    !!{
    Attach an object to an event hook.
    !!}
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
"""


_ISATTACHED_TEMPLATE = """  logical function eventHook{$interfaceType}IsAttached(self,object_,function_)
    !!{
    Return true if an object is attached to an event hook.
    !!}
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
"""


def _parse_interface_arguments(interface_text):
    """Extract argument names from the Fortran declarations under a hook's
    `<interface>` tag.

    Mirrors EventHooks.pm:46-59 in spirit: iterate each Fortran line, match
    against the intrinsic-type patterns, and extract the variable list after
    the `::` separator.

    Note that an attribute list may contain `dimension(:)` (with a colon
    inside the parens), so the line cannot be split on the first `::` found
    by `[^:]*::` — that pattern would lock onto the colon inside `(:)` and
    then fail to match the real `::` that follows.  Use a leading-intrinsic
    check, then split on the *last* `::` to isolate the trailing variable
    list.
    """
    arguments = []
    intrinsic_re = re.compile(
        r'^\s*(?:integer|real|double\s+precision|double\s+complex|'
        r'logical|character|type|class)\b',
        re.IGNORECASE,
    )

    fh = io.StringIO(interface_text)
    while True:
        raw_line, processed_line, _ = get_fortran_line(fh)
        if not raw_line and not processed_line:
            break
        if not intrinsic_re.match(processed_line):
            continue
        sep = processed_line.rfind('::')
        if sep < 0:
            continue
        var_list = processed_line[sep + 2:].strip()
        if var_list:
            arguments.extend(extract_variables(var_list, keep_qualifiers=False))
    return arguments


def _render_imports_and_uses(hook, manager_parent):
    """Render the `  import a, b, c\\n` line for the abstract interface body,
    and inject the matching `use … only : …` on the manager's parent module.

    Returns the import-line text (empty string when no `<import>` declared).
    """
    if 'import' not in hook:
        return ''
    imports_block = hook['import']
    module_block  = imports_block.get('module') if isinstance(imports_block, dict) else None
    if module_block is None:
        return ''
    # Perl supports both `{module: {name, symbols}}` (single) and
    # `{module: {fooMod: {symbols}, barMod: {symbols}}}` (multi — keyAttr-style
    # keyed by name).  Our xml_to_dict emits a list or a single dict, not the
    # keyAttr grouping, so we normalise via `as_array`.
    if isinstance(module_block, dict) and 'name' in module_block:
        modules = [module_block]
    elif isinstance(module_block, dict):
        # dict keyed by name (one-level Perl-style grouping).
        modules = [
            dict(attrs, name=name) if isinstance(attrs, dict) else {'name': name}
            for name, attrs in module_block.items()
        ]
    else:
        modules = list(as_array(module_block))
    all_symbols = []
    module_uses = {}
    for mod in modules:
        if not isinstance(mod, dict) or 'name' not in mod:
            continue
        symbols = [s.strip() for s in re.split(
            r'\s*,\s*', mod.get('symbols', '')) if s.strip()]
        all_symbols.extend(symbols)
        module_uses[mod['name']] = {
            'intrinsic': False,
            'only':      {s: True for s in symbols},
        }
    if module_uses:
        add_uses(manager_parent, {
            'moduleUse':   module_uses,
            'moduleOrder': list(module_uses.keys()),
        })
    if not all_symbols:
        return ''
    return "  import " + ", ".join(all_symbols) + " \n"


def _emit_typed_hook_block(hook, parent_node, manager_parent):
    """Emit the per-interface-type hook block: hook type, eventHook type,
    abstract interface, and the attach/detach/isAttached subroutines.
    """
    interface_type = _interface_type_get(hook)
    loc_expr = location(parent_node, parent_node.get('line', 0))
    declarations = hook.get('interface', '')
    arguments = _parse_interface_arguments(declarations)
    argument_list = ',' + ','.join(arguments) if arguments else ''
    imports_text = _render_imports_and_uses(hook, manager_parent)

    subs = {
        'interfaceType': interface_type,
        'argumentList':  argument_list,
        'imports':       imports_text,
        'declarations':  declarations,
        'location':      loc_expr,
    }

    # The type block itself stays near the directive site (Perl uses
    # InsertAfterNode for it alongside the module-level event variables).
    type_block = _substitute(_HOOK_TYPE_TEMPLATE, subs)

    # The subroutine bodies go after the `contains` marker of the manager's
    # parent module.  Run process_tree on each parsed sub_tree so directives
    # embedded in the templates (e.g. `<optionalArgument>` in `_ATTACH_TEMPLATE`)
    # are expanded and marked processed — the dispatchers for those directive
    # types have already run by the time eventHooks executes.
    for template in (_ATTACH_TEMPLATE, _DETACH_TEMPLATE, _ISATTACHED_TEMPLATE):
        sub_tree = parse_code(
            _substitute(template, subs),
            name='EventHooks',
        )
        process_tree(sub_tree)
        kids = children(sub_tree)
        for k in kids:
            k['parent'] = None
        insert_post_contains(manager_parent, kids)

    set_visibility(manager_parent, 'hook' + interface_type, 'public')
    return type_block


# ---------------------------------------------------------------------------
# Module-level subroutines synthesized once per manager directive.
# ---------------------------------------------------------------------------

_COPY_OUT_HEAD = """subroutine eventsHooksFilterCopyOut_()
   implicit none
   call copyLock%set()
"""
_COPY_OUT_TAIL = """   return
end subroutine eventsHooksFilterCopyOut_
"""


_COPY_IN_HEAD = """subroutine eventsHooksFilterCopyIn_()
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
"""
_COPY_IN_PER_HOOK = """   allocate(eventHookBackup)
   allocate(eventHookBackup%eventHook_,mold={$name}Event_)
   {$name}Event               =  {$name}Event_
   eventHookBackup%eventHook_ =  {$name}Event_
   if (associated({$name}EventBackups)) eventHookBackup%next => {$name}EventBackups
   {$name}EventBackups        => eventHookBackup
   message=var_str("{$name}: storing ")//eventHookBackup%eventHook_%count_//" hooks"//threadLabel
   call displayMessage(message,verbosityLevelInfo)
   nullify(eventHookBackup)
"""
_COPY_IN_TAIL = """   return
end subroutine eventsHooksFilterCopyIn_
"""


_RESTORE_HEAD = """subroutine eventsHooksFilterRestore_()
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
"""
_RESTORE_PER_HOOK = """   eventHookBackup     => {$name}EventBackups
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
"""
_RESTORE_TAIL = """   return
end subroutine eventsHooksFilterRestore_
"""


_COPY_DONE = """subroutine eventsHooksFilterCopyDone_()
   implicit none
   call copyLock%unset()
   return
end subroutine eventsHooksFilterCopyDone_
"""


_FILTER_HEAD = """subroutine eventsHooksFilterFunction_()
   implicit none
"""
_FILTER_TAIL = """   return
end subroutine eventsHooksFilterFunction_
"""


_INITIALIZE_HEAD = """subroutine eventsHooksInitialize()
   use :: Events_Filters, only : eventsHooksFilterFunction, eventsHooksFilterCopyOut, eventsHooksFilterCopyIn, eventsHooksFilterCopyDone, &
        &                        eventsHooksFilterRestore
   implicit none

   eventsHooksFilterFunction => eventsHooksFilterFunction_
   eventsHooksFilterCopyOut  => eventsHooksFilterCopyOut_
   eventsHooksFilterCopyIn   => eventsHooksFilterCopyIn_
   eventsHooksFilterCopyDone => eventsHooksFilterCopyDone_
   eventsHooksFilterRestore  => eventsHooksFilterRestore_
   copyLock=ompLock()
"""
_INITIALIZE_TAIL = """   return
end subroutine eventsHooksInitialize
"""


_WAIT_TIMES_HEAD = """subroutine eventsHooksWaitTimes()
#ifdef OMPPROFILE
    use :: Output_HDF5       , only : outputFile
    use :: IO_HDF5           , only : hdf5Object
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : varying_string      , var_str
    use :: Units_MetaData    , only : unitType
#endif
    implicit none
#ifdef OMPPROFILE
    type            (hdf5Object                  )                          :: waitTimeGroup         , waitTimeDataset        , metaDataGroup
    character       (len={$hookNameLengthMaximum}), dimension({$hookCount}) :: eventHookNames
    double precision                              , dimension({$hookCount}) :: eventHookReadWaitTimes, eventHookWriteWaitTimes

"""
_WAIT_TIMES_TAIL = """    ! Open output group.
    !$ call hdf5Access%set()
    metaDataGroup=outputFile%openGroup('metaData','Galacticus meta data.'           )
    waitTimeGroup=metaDataGroup       %openGroup('openMP'  ,'Meta-data on OpenMP performance.')
    ! Write wait time data.
    call waitTimeGroup%writeDataset(eventHookNames         ,"eventHookNames"         ,"Names of event hooks"                                                              )
    call waitTimeGroup%writeDataset(eventHookReadWaitTimes ,"eventHookReadWaitTimes" ,"Total time spent waiting to read-lock event hooks" ,datasetReturned=waitTimeDataset)
    call waitTimeDataset%writeAttribute(unitType(1.0d0,"seconds","s"),"units")
    call waitTimeDataset%close()
    call waitTimeGroup%writeDataset(eventHookWriteWaitTimes,"eventHookWriteWaitTimes","Total time spent waiting to write-lock event hooks",datasetReturned=waitTimeDataset)
    call waitTimeDataset%writeAttribute(unitType(1.0d0,"seconds","s"),"units")
    call waitTimeDataset%close()
    ! Close output groups.
    call waitTimeGroup%close()
    call metaDataGroup%close()
    !$ call hdf5Access%unset()
#endif
   return
end subroutine eventsHooksWaitTimes
"""


def _emit_module_level_subs(hooks, parent_node, manager_parent):
    """Synthesize every module-level helper subroutine (copy-out, copy-in,
    restore, copy-done, filter, initialize, wait-times) from the hook list.
    """
    loc_expr = location(parent_node, parent_node.get('line', 0))

    # copy-out
    code = _COPY_OUT_HEAD
    for h in hooks:
        code += f"   {h['name']}Event_={h['name']}Event\n"
    code += _COPY_OUT_TAIL
    insert_post_contains(manager_parent, [_code_node(code)])

    # copy-in
    code = _COPY_IN_HEAD
    for h in hooks:
        code += _substitute(_COPY_IN_PER_HOOK, {'name': h['name']})
    code += _COPY_IN_TAIL
    insert_post_contains(manager_parent, [_code_node(code)])

    # restore
    code = _RESTORE_HEAD
    for h in hooks:
        code += _substitute(_RESTORE_PER_HOOK, {
            'name':          h['name'],
            'interfaceType': _interface_type_get(h),
            'location':      loc_expr,
        })
    code += _RESTORE_TAIL
    insert_post_contains(manager_parent, [_code_node(code)])

    # copy-done
    insert_post_contains(manager_parent, [_code_node(_COPY_DONE)])

    # filter
    code = _FILTER_HEAD
    for h in hooks:
        code += f"   call {h['name']}Event%filter()\n"
    code += _FILTER_TAIL
    insert_post_contains(manager_parent, [_code_node(code)])

    # initializer
    code = _INITIALIZE_HEAD
    for h in hooks:
        name = h['name']
        code += f"   {name}EventGlobal%isGlobal=.true.\n"
        code += f"   {name}EventGlobal%lock_   =ompReadWriteLock()\n"
        code += f"   {name}Event      %label   ='{name}'\n"
        code += f"   {name}EventGlobal%label   ='{name} (global)'\n"
    code += _INITIALIZE_TAIL
    insert_post_contains(manager_parent, [_code_node(code)])
    set_visibility(manager_parent, 'eventsHooksInitialize', 'public')

    # wait-times writer (OMPPROFILE-gated).
    max_name_len = max((len(h['name']) for h in hooks), default=0)
    code = _substitute(_WAIT_TIMES_HEAD, {
        'hookCount':             len(hooks),
        'hookNameLengthMaximum': max_name_len,
    })
    for i, h in enumerate(hooks, start=1):
        code += f"   eventHookNames         ({i})='{h['name']}'\n"
        code += f"   eventHookReadWaitTimes ({i})={h['name']}Event%waitTimeRead\n"
        code += f"   eventHookWriteWaitTimes({i})={h['name']}Event%waitTimeWrite\n"
    code += _WAIT_TIMES_TAIL
    insert_post_contains(manager_parent, [_code_node(code)])
    set_visibility(manager_parent, 'eventsHooksWaitTimes', 'public')


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

def _process_event_hook_manager(node):
    """Handle a single `<eventHookManager/>` directive.

    Mirrors EventHooks.pm:33-561.  For each non-duplicate `<eventHook>`
    discovered across all files in directiveLocations/eventHook/file:
      - If the hook declares an `<interface>`, emit the dedicated
        hook{InterfaceType} / eventHook{InterfaceType} types plus the
        attach/detach/isAttached subroutines (once per unique interface
        signature — though in practice every hook with its own interface
        gets its own md5-hex type name).
      - Always emit the module-level event variables + threadprivate
        declaration at the directive site.
    Then synthesize the module-level helper subroutines (copy-out/in/done,
    restore, filter, initializer, wait-times).
    """
    directive = node.setdefault('directive', {})
    directive['processed'] = True
    manager_parent = node['parent']

    directive_locations = _load_directive_locations()
    hooks = _collect_hooks(directive_locations)

    for hook in hooks:
        type_block = ''
        interface_type = _interface_type_get(hook)
        if interface_type != 'Unspecified':
            type_block = _emit_typed_hook_block(hook, node, manager_parent)

        # Module-level event variables for this hook.
        type_block += (
            f"type(eventHook{interface_type}), public  :: {hook['name']}Event"
            f"                 , {hook['name']}EventGlobal\n"
            f"type(eventHook{interface_type})          :: {hook['name']}Event_\n"
            f"type(eventHookList                    ), pointer :: "
            f"{hook['name']}EventBackups => null()\n"
            f"!$omp threadprivate ({hook['name']}Event,"
            f"{hook['name']}EventBackups)\n"
        )
        sub_tree = parse_code(type_block, name='EventHooks')
        kids = children(sub_tree)
        for k in kids:
            k['parent'] = None
        insert_after_node(node, kids)

    _emit_module_level_subs(hooks, node, manager_parent)


def process_event_hooks(tree, options):
    """Mirrors Process_EventHooks() from EventHooks.pm."""
    # Materialise the walk up front because we mutate the tree as we go.
    for node in list(walk_tree(tree)):
        ntype = node.get('type')
        # Read-only access to `directive` here; we only call setdefault inside
        # the per-type handlers that *will* write back to it, so we don't
        # accidentally add an empty `'directive': {}` key to every node in the
        # tree (which would later trip post_process_directives).
        directive = node.get('directive') or {}

        if ntype == 'eventHookManager' and not directive.get('processed'):
            _process_event_hook_manager(node)
            continue

        if ntype == 'eventHook' and not directive.get('processed'):
            _process_event_hook_call_site(node)
            continue

        if (ntype == 'openMP'
                and node.get('name') == 'parallel'
                and 'eventFilterInserted' not in node):
            if not node.get('isCloser'):
                _process_openmp_opener(node)
            else:
                _process_openmp_closer(node)


register_process('eventHooks', process_event_hooks)
