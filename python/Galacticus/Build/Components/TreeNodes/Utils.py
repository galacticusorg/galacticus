"""`treeNode` utility methods: copyNodeTo, moveComponentsTo, massDistribution,
massBaryonic.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/TreeNodes/Utils.pm.  Four
`functions`-phase hooks.  This file lands in three commit-sized
steps; this step ports `Copy` and `Move` only.  The remaining two
hooks are stubbed and wired up in follow-up commits.
"""



from Galacticus.Build.Components.Utils import register


_COPY_VALUES = (
    'uniqueIdValue', 'indexValue', 'timeStepValue', 'subsamplingWeightValue',
)
_COPY_POINTERS = (
    'parent', 'firstChild', 'sibling', 'firstSatellite', 'mergeTarget',
    'firstMergee', 'siblingMergee', 'event', 'formationNode', 'hostTree',
)


def Tree_Node_Copy(build):
    """Generate `treeNodeCopyNodeTo`.  Mirrors `Tree_Node_Copy`.

    Copies value-typed fields and pointer-typed fields from `self` onto
    `targetNode`, then re-binds each component's `hostNode` pointer back
    to `targetNode`.  Two optional flags suppress the formation-node
    and event copies.
    """
    function = {
        'type':        'void',
        'name':        'treeNodeCopyNodeTo',
        'description': r"Make a copy of \mono{self} in \mono{targetNode}.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['targetNode'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['skipFormationNode', 'skipEvent'],
            },
            {
                'intrinsic':  'logical',
                'variables':  ['skipFormationNodeActual', 'skipEventActual'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = (
        "skipFormationNodeActual=.false.\n"
        "skipEventActual        =.false.\n"
        "if (present(skipFormationNode)) skipFormationNodeActual=skipFormationNode\n"
        "if (present(skipEvent        )) skipEventActual        =skipEvent\n"
    )
    content += ''.join(
        f"targetNode%{v} =  self%{v}\n" for v in _COPY_VALUES
    )
    content += ''.join(
        f"targetNode%{p} => self%{p}\n" for p in _COPY_POINTERS
    )
    content += (
        "if (skipFormationNodeActual) targetNode%formationNode => null()\n"
        "if (skipEventActual        ) targetNode%event         => null()\n"
    )
    for class_name in build.get('componentClassListActive') or []:
        cap = _ucfirst(class_name)
        content += f"targetNode%component{cap} = self%component{cap}\n"

    content += (
        "select type (targetNode)\n"
        "type is (treeNode)\n"
    )
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"   do i=1,size(self%component{cap})\n"
            f"     targetNode%component{cap}(i)%hostNode => targetNode\n"
            f"   end do\n"
        )
    content += "end select\n"

    function['content'] = content

    # Note: the Perl original adds `returnType` and `arguments` keys to the
    # binding (TreeNodes/Utils.pm:103-104).  The base framework writes those
    # out as part of the documentation block — match it here.
    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':        'procedure',
        'descriptor':  function,
        'name':        'copyNodeTo',
        'returnType':  r"\void",
        'arguments':   (
            r"\textcolor{red}{\textless class(treeNode)\textgreater} "
            r"targetNode\arginout, \logicalzero\ [skipFormationNode]\argin"
        ),
    })


def Tree_Node_Move(build):
    """Generate `treeNodeComponentsMove`.  Mirrors `Tree_Node_Move`.

    For each active class, destroys any existing `targetNode%component<X>`
    array, then `move_alloc`s `self%component<X>` into it and re-binds
    each component's `hostNode` pointer to `targetNode`.
    """
    function = {
        'type':        'void',
        'name':        'treeNodeComponentsMove',
        'description': (
            r"Move components from \mono{self} to \mono{targetNode}."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['intent(inout)', 'target'],
                'variables':  ['targetNode'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = ''
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(targetNode%component{cap})) then\n"
            f"  do i=1,size(targetNode%component{cap})\n"
            f"    call targetNode%component{cap}(i)%destroy()\n"
            f"  end do\n"
            f"  deallocate(targetNode%component{cap})\n"
            f"end if\n"
            f"if (allocated(self      %component{cap})) then\n"
            f"   call Move_Alloc(self%component{cap},"
            f"targetNode%component{cap})\n"
            f"   do i=1,size(targetNode%component{cap})\n"
            f"     targetNode%component{cap}(i)%hostNode => targetNode\n"
            f"   end do\n"
            f"end if\n"
        )
    function['content'] = content

    _bind(build, 'moveComponentsTo', function)


def Tree_Node_Mass_Distribution(build):
    """Generate `treeNodeMassDistribution`.  Mirrors `Tree_Node_Mass_Distribution`.

    Bulk of the body is a static template that maintains a small
    `massDistributions__` cache keyed by `(uniqueID, componentType,
    massType, weightBy, weightIndex)`.  The only per-class section is
    the inner loop that gathers each component's
    `massDistribution(...)` into a `massDistributionList_` linked list.
    """
    function = {
        'type':        'class(massDistributionClass), pointer => massDistribution_',
        'name':        'treeNodeMassDistribution',
        'description': (
            r"Construct and return the mass distribution associated with "
            r"\mono{self}."
        ),
        'modules':     [
            ("Mass_Distributions        , only : massDistributionClass"
             "       , massDistributionComposite, massDistributionList"
             "   , massDistributionZero, kinematicsDistributionClass"
             ", kinematicsDistributionIsothermal"),
            ("Galactic_Structure_Options, only : enumerationComponentTypeType"
             ", enumerationMassTypeType  , enumerationWeightByType"
             ", componentTypeAll    , componentTypeDarkMatterOnly"
             ", massTypeAll                      , massTypeDark, weightByMass"),
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationComponentTypeType',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['componentType'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationMassTypeType',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['massType'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationWeightByType',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['weightBy'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['weightIndex'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationComponentTypeType',
                'variables':  ['componentType_', 'componentType__'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationMassTypeType',
                'variables':  ['massType_', 'massType__'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationWeightByType',
                'variables':  ['weightBy_'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['weightIndex_'],
            },
            {
                'intrinsic':  'class',
                'type':       'massDistributionClass',
                'attributes': ['pointer'],
                'variables':  ['massDistributionComponent'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'iMassDistribution', 'iMassDistributionAll',
                              'iMassDistributionConstruct', 'iEmpty'],
            },
            {
                'intrinsic':  'logical',
                'variables':  ['construct', 'isDarkMatterOnly'],
            },
            {
                'intrinsic':  'integer',
                'type':       'kind_int8',
                'variables':  ['uniqueID', 'uniqueIDParent'],
            },
            {
                'intrinsic':  'type',
                'type':       'massDistributionList',
                'attributes': ['pointer'],
                'variables':  ['massDistributionList_', 'massDistributionListCopy_',
                              'next_', 'nextCopy_'],
            },
            {
                'intrinsic':  'class',
                'type':       'kinematicsDistributionClass',
                'attributes': ['pointer'],
                'variables':  ['kinematicsDistribution_'],
            },
        ],
    }

    # Per-class inner section.  Wedged between the `head_template` and
    # `tail_template` strings.
    per_class = ''
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        per_class += (
            f"  if (allocated(self%component{cap})) then\n"
            f"     do i=1,size(self%component{cap})\n"
            f"       massDistributionComponent => self%component{cap}(i)"
            f"%massDistribution(componentType__,massType__,weightBy_,weightIndex_)\n"
            "       if (associated(massDistributionComponent)) then\n"
            "          if (associated(massDistributionList_)) then\n"
            "            allocate(next_    %next)\n"
            "            allocate(nextCopy_%next)\n"
            "            next_     => next_    %next\n"
            "            nextCopy_ => nextCopy_%next\n"
            "          else\n"
            "            allocate(massDistributionList_    )\n"
            "            allocate(massDistributionListCopy_)\n"
            "            next_     => massDistributionList_\n"
            "            nextCopy_ => massDistributionListCopy_\n"
            "          end if\n"
            "          next_    %massDistribution_ => massDistributionComponent\n"
            "          nextCopy_%massDistribution_ => massDistributionComponent\n"
            "          next_    %next              => null()\n"
            "          nextCopy_%next              => null()\n"
            "       end if\n"
            "     end do\n"
            "  end if\n"
        )

    function['content'] = (
        _MASS_DISTRIBUTION_HEAD + per_class + _MASS_DISTRIBUTION_TAIL
    )
    _bind(build, 'massDistribution', function)


# ---------------------------------------------------------------------------
# Static body fragments for `Tree_Node_Mass_Distribution`.  Verbatim from
# Utils.pm:264-368 (head) and Utils.pm:398-496 (tail).  The only dynamic
# block is the per-class loop above.
# ---------------------------------------------------------------------------

_MASS_DISTRIBUTION_HEAD = """! Set defaults.
if (present(componentType)) then
 componentType_=componentType
else
 componentType_=componentTypeAll
end if
if (present(massType)) then
 massType_=massType
else
 massType_=massTypeAll
end if
if (present(weightBy)) then
 weightBy_=weightBy
else
 weightBy_=weightByMass
end if
if (present(weightIndex)) then
 weightIndex_=weightIndex
else
 weightIndex_=-1
end if
! Search for a match to our ID.
iMassDistribution   =0
iMassDistributionAll=0
iEmpty              =0
uniqueID            =self       %uniqueID()
if (associated(self%parent)) then
 uniqueIDParent  =self%parent%uniqueID()
else
 uniqueIDParent  =-1_kind_int8
end if
do i=1,massDistributionsCount
 if   (                                                        &
    &   massDistributions__(i)%uniqueID      == uniqueID       &
    &  .and.                                                   &
    &   massDistributions__(i)%componentType == componentType_ &
    &  .and.                                                   &
    &   massDistributions__(i)%massType      == massType_      &
    &  .and.                                                   &
    &   massDistributions__(i)%weightBy      == weightBy_      &
    &  .and.                                                   &
    &   massDistributions__(i)%weightIndex   == weightIndex_   &
    & ) then
  iMassDistribution=i
  exit
 else if (massDistributions__(i)%uniqueID < 0_kind_int8 .and. iEmpty == 0) then
  iEmpty           =i
 end if
 if   (                                                          &
    &   massDistributions__(i)%uniqueID      == uniqueID         &
    &  .and.                                                     &
    &   massDistributions__(i)%componentType == componentTypeAll &
    &  .and.                                                     &
    &   massDistributions__(i)%massType      == massTypeAll      &
    &  .and.                                                     &
    &   massDistributions__(i)%weightBy      == weightBy_        &
    &  .and.                                                     &
    &   massDistributions__(i)%weightIndex   == weightIndex_     &
    & ) then
  iMassDistributionAll=i
 end if
end do
! If we found no match, we need to create the distribution.
construct=.false.
if (iMassDistribution == 0) then
 isDarkMatterOnly=.false.
 if      (componentType_       == componentTypeDarkMatterOnly) then
  isDarkMatterOnly=.true.
  construct       =.true.
  componentType__ =componentTypeDarkMatterOnly
  massType__      =massType_
 else if (iMassDistributionAll == 0                          ) then
  construct       =.true.
  componentType__ =componentTypeAll
  massType__      =massTypeAll
 end if
 ! If no existing all/all mass distribution matched.....
 if (construct) then
  if      (iEmpty     /= 0) then
   ! If we have an empty slot, use that.
   iMassDistributionConstruct=iEmpty
  else
   ! Simply use the next slot, unless it is occupied by a parent node massDistribution (unless we have no choice because we have run out of slots).
   do i=1,massDistributionsCount+1
    massDistributionsLast=mod(massDistributionsLast,massDistributionsCount)+1
    if (massDistributions__(massDistributionsLast)%uniqueID == uniqueIDParent) cycle
    exit
   end do
   iMassDistributionConstruct=    massDistributionsLast
   !![
   <objectDestructor name="massDistributions__(iMassDistributionConstruct)%massDistribution_"/>
   !!]
   massDistributions__(massDistributionsLast)%uniqueID=-huge(kind_int8)
  end if
  if (isDarkMatterOnly) then
   iMassDistribution   =iMassDistributionConstruct
  else
   iMassDistributionAll=iMassDistributionConstruct
  end if
  if (.not.associated(massDistributions__(iMassDistributionConstruct)%massDistribution_)) then
   massDistributionList_     => null()
   massDistributionListCopy_ => null()
   next_                     => null()
   nextCopy_                 => null()
"""


_MASS_DISTRIBUTION_TAIL = """   allocate(massDistributionComposite :: massDistributions__(iMassDistributionConstruct)%massDistribution_)
   select type (massDistribution__ => massDistributions__(iMassDistributionConstruct)%massDistribution_)
   type is (massDistributionComposite)
     !![
     <referenceConstruct isResult="yes" object="massDistribution__" constructor="massDistributionComposite(massDistributionList_)"/>
     !!]
   end select
   massDistributions__(iMassDistributionConstruct)%uniqueID     =self%uniqueID        ()
   massDistributions__(iMassDistributionConstruct)%componentType=     componentType__
   massDistributions__(iMassDistributionConstruct)%massType     =     massType__
   massDistributions__(iMassDistributionConstruct)%weightBy     =     weightBy_
   massDistributions__(iMassDistributionConstruct)%weightIndex  =     weightIndex_
   next_ => massDistributionListCopy_
   do while (associated(next_))
    !![
    <objectDestructor name="next_%massDistribution_" nullify="no"/>
    !!]
    nextCopy_ => next_%next
    deallocate(next_)
    next_ => nextCopy_
   end do
   nullify(massDistributionList_)
  end if
 end if
 if (isDarkMatterOnly) then
  ! We already have the relevant mass distribution constructed in the iMassDistribution slot - nothing more to do here.
 else if (componentType_ == componentTypeAll .and. massType_ == massTypeAll) then
  ! The all/all mass distribution was required - we have just created it, so return it.
  iMassDistribution=iMassDistributionAll
 else
  ! Some other mass distribution was required - get it as a subset of the all/all mass distribution.
  iEmpty=0
  do i=1,massDistributionsCount
   if (massDistributions__(i)%uniqueID < 0_kind_int8 .and. iEmpty == 0) iEmpty=i
  end do
  if (iEmpty /= 0) then
   ! If we have an empty slot, use that.
   iMassDistribution=iEmpty
  else
   ! Simply use the next slot, unless it is occupied by a parent node massDistribution (unless we have no choice because we have run out of slots).
   do i=1,massDistributionsCount+1
    massDistributionsLast=mod(massDistributionsLast,massDistributionsCount)+1
    massDistributionsLast=mod(massDistributionsLast,massDistributionsCount)+1
    ! But never replace the all/all distribution.
    if   (                                                                              &
       &   massDistributions__(massDistributionsLast)%uniqueID      == uniqueID         &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%componentType == componentTypeAll &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%massType      == massTypeAll      &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%weightBy      == weightBy_        &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%weightIndex   == weightIndex_     &
       & ) cycle
    if   (                                                                              &
       &   massDistributions__(massDistributionsLast)%uniqueID      == uniqueIDParent   &
       & ) cycle
    exit
   end do
   iMassDistribution=massDistributionsLast
   !![
   <objectDestructor name="massDistributions__(iMassDistribution)%massDistribution_"/>
   !!]
   massDistributions__(massDistributionsLast)%uniqueID=-huge(kind_int8)
  end if
  massDistributions__(iMassDistribution)%massDistribution_ => massDistributions__(iMassDistributionAll)%massDistribution_%subset        (componentType_,massType_)
  massDistributions__(iMassDistribution)%uniqueID          =  self                                                       %uniqueID      (                        )
  massDistributions__(iMassDistribution)%componentType     =                                                              componentType_
  massDistributions__(iMassDistribution)%massType          =                                                              massType_
  massDistributions__(iMassDistribution)%weightBy          =                                                              weightBy_
  massDistributions__(iMassDistribution)%weightIndex       =                                                              weightIndex_
  if (.not.associated(massDistributions__(iMassDistribution)%massDistribution_)) then
   allocate(massDistributionZero :: massDistributions__(iMassDistribution)%massDistribution_)
   select type (massDistributions___ => massDistributions__(iMassDistribution)%massDistribution_)
   type is (massDistributionZero)
    !![
    <referenceConstruct object="massDistributions___" constructor="massDistributionZero(dimensionless=.false.)"/>
    !!]
    allocate(kinematicsDistributionIsothermal :: kinematicsDistribution_)
    select type (kinematicsDistribution_)
    type is (kinematicsDistributionIsothermal)
     !![
     <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionIsothermal(velocityDispersion_=0.0d0)"/>
     !!]
    end select
    call massDistributions___%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
   end select
  end if
 end if
end if
!![
<referenceAcquire target="massDistribution_" source="massDistributions__(iMassDistribution)%massDistribution_"/>
!!]
"""


# ---------------------------------------------------------------------------
# Helpers (shared with sister modules)
# ---------------------------------------------------------------------------

def _active_classes(build):
    active = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        if class_dict['name'] in active:
            yield class_dict


def _bind(build, method_name, function):
    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


def Tree_Node_Mass_Baryonic(build):
    """Generate `treeNodeMassBaryonic`.  Mirrors `Tree_Node_Mass_Baryonic`.

    Returns the sum of `massBaryonic()` over every component instance on
    every active class.
    """
    function = {
        'type':        'double precision',
        'name':        'treeNodeMassBaryonic',
        'description': (
            r"Return the total baryonic mass associated with \mono{self}."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = "treeNodeMassBaryonic=0.0d0\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"  if (allocated(self%component{cap})) then\n"
            f"     do i=1,size(self%component{cap})\n"
            f"       treeNodeMassBaryonic=treeNodeMassBaryonic"
            f"+self%component{cap}(i)%massBaryonic()\n"
            f"     end do\n"
            f"  end if\n"
        )
    function['content'] = content

    _bind(build, 'massBaryonic', function)


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Utils.pm:21-25.
# ---------------------------------------------------------------------------

register('treeNodeUtils', 'functions', Tree_Node_Copy)
register('treeNodeUtils', 'functions', Tree_Node_Move)
register('treeNodeUtils', 'functions', Tree_Node_Mass_Distribution)
register('treeNodeUtils', 'functions', Tree_Node_Mass_Baryonic)
