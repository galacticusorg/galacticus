"""Linked-list code generation helpers for the functionClass pipeline.

Andrew Benson (ported to Python 2026)
"""



from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


def _substitute(template, subs):
    """Replace `{$key}` placeholders in `template` with values from `subs`.

    Same targeted substituter used by Process/EventHooks.py — deliberately
    avoids pulling in a templating engine.
    """
    out = template
    for key, value in subs.items():
        out = out.replace('{$' + key + '}', str(value))
    return out


def linked_list_register_variable(linked_list, linked_list_variables,
                                  *variable_names):
    """Add a `type({ll.type}), pointer :: <names>` declaration to
    `linked_list_variables`, unless a declaration of the same type already
    exists.
    """
    for existing in linked_list_variables:
        if (existing.get('type') is not None
                and existing['type'] == linked_list['type']):
            return
    linked_list_variables.append({
        'intrinsic':  'type',
        'type':       linked_list['type'],
        'attributes': ['pointer'],
        'variables':  list(variable_names),
    })


def linked_list_module(linked_list):
    """Return the `module` attribute of the linked-list spec, or None."""
    return linked_list.get('module')


# ---------------------------------------------------------------------------
# deepCopy
# ---------------------------------------------------------------------------

_DC_INIT = """{$type}item             => destination%{$variable}
do while (associated({$type}item))
   ! Undo the reference count increment that resulted from the initial intrinsic assignment.
   referenceCount___={$type}item%{$object}%referenceCountDecrement()
   nullify({$type}item%{$object})
   {$type}itemNew => {$type}item%next
   deallocate({$type}item)
   {$type}item => {$type}itemNew
end do
destination%{$variable} => null            ()
"""

_DC_WALK_HEAD = """{$type}destination      => null            ()
{$type}item             => self%{$variable}
do while (associated({$type}item))
"""

_DC_FIRST_OBJECT = """   allocate({$type}itemNew)
   if (associated({$type}destination)) then
      {$type}destination%{$next}     => {$type}itemNew
      {$type}destination             => {$type}itemNew
   else
      destination       %{$variable} => {$type}itemNew
      {$type}destination             => {$type}itemNew
   end if
"""

_DC_LATER_OBJECT = """   if (associated({$type}destination)) then
      {$type}itemNew     => {$type}destination%{$next}
      {$type}destination => {$type}destination%{$next}
   else
      {$type}itemNew     => destination       %{$variable}
      {$type}destination => destination       %{$variable}
   end if
"""

_DC_BODY = """      nullify({$type}itemNew%{$object})
      if (associated({$type}item%{$object})) then
       if (associated({$type}item%{$object}%copiedSelf)) then
        select type(s => {$type}item%{$object}%copiedSelf)
        class is ({$objectType})
         {$type}itemNew%{$object} => s
        class default
         call Error_Report('copiedSelf has incorrect type'//{$location})
        end select
        call {$type}item%{$object}%copiedSelf%referenceCountIncrement()
       else
        allocate({$type}itemNew%{$object},mold={$type}item%{$object})
        call {$type}item%{$object}%deepCopy({$type}itemNew%{$object})
        {$type}item%{$object}%copiedSelf => {$type}itemNew%{$object}
        call {$type}itemNew%{$object}%autoHook()
       end if
      end if
   {$type}item => {$type}item%{$next}
end do
"""

_DC_RESET = """{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%deepCopyReset()
   {$type}item => {$type}item%{$next}
end do
"""

_DC_FINALIZE = """{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%deepCopyFinalize()
   {$type}item => {$type}item%{$next}
end do
"""


def deep_copy_linked_list(class_record, non_abstract_class,
                          linked_list_variables,
                          linked_list_reset_variables,
                          linked_list_finalize_variables):
    """Generate deep-copy / reset / finalize Fortran code for the class's
    `<linkedList>` block.

    Returns `(deepCopy, deepCopyReset, deepCopyFinalize, module)`.  Every
    returned code string is empty (and module is None) when the class has
    no `<linkedList>` metadata.
    """
    if 'linkedList' not in class_record:
        return '', '', '', None

    linked_list = class_record['linkedList']
    objects      = (linked_list.get('object')     or '').split()
    object_types = (linked_list.get('objectType') or '').split()

    # Variables used by the deep-copy walk.
    linked_list_register_variable(
        linked_list, linked_list_variables,
        linked_list['type'] + 'item',
        linked_list['type'] + 'destination',
        linked_list['type'] + 'itemNew',
    )
    if not any(
            (v.get('variables') or [None])[0] == 'referenceCount___'
            for v in linked_list_variables):
        linked_list_variables.append({
            'intrinsic': 'integer',
            'variables': ['referenceCount___'],
        })
    linked_list_register_variable(
        linked_list, linked_list_reset_variables,
        linked_list['type'] + 'item')
    linked_list_register_variable(
        linked_list, linked_list_finalize_variables,
        linked_list['type'] + 'item')

    deep_copy_code          = ''
    deep_copy_reset_code    = ''
    deep_copy_finalize_code = ''
    node = class_record.get('node') or {}
    loc_expr = location(node, node.get('line', 0))

    for i, obj in enumerate(objects):
        object_type = object_types[i] if i < len(object_types) else ''
        subs = {
            'type':       linked_list['type'],
            'variable':   linked_list['variable'],
            'next':       linked_list['next'],
            'object':     obj,
            'objectType': object_type,
            'location':   loc_expr,
        }
        if i == 0:
            deep_copy_code += _substitute(_DC_INIT, subs)
        deep_copy_code += _substitute(_DC_WALK_HEAD, subs)
        deep_copy_code += _substitute(
            _DC_FIRST_OBJECT if i == 0 else _DC_LATER_OBJECT, subs)
        deep_copy_code += _substitute(_DC_BODY, subs)
        deep_copy_reset_code    += _substitute(_DC_RESET, subs)
        deep_copy_finalize_code += _substitute(_DC_FINALIZE, subs)

    return (deep_copy_code, deep_copy_reset_code, deep_copy_finalize_code,
            linked_list_module(linked_list))


# ---------------------------------------------------------------------------
# stateStore / stateRestore
# ---------------------------------------------------------------------------

_SS_RESTORE = """{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%stateRestore(stateFile,gslStateFile,stateOperationID)
   {$type}item => {$type}item%{$next}
end do
"""

_SS_STORE = """{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%stateStore(stateFile,gslStateFile,stateOperationID)
   {$type}item => {$type}item%{$next}
end do
"""


def state_store_linked_list(class_record, non_abstract_class,
                            linked_list_variables):
    """Generate state-store / state-restore Fortran code for the linked-list
    block.  Returns `(input_code, output_code, module)`; the no-op path is
    normalised to the same 3-tuple shape.
    """
    if 'linkedList' not in class_record:
        return '', '', None
    linked_list = class_record['linkedList']
    objects = (linked_list.get('object') or '').split()
    linked_list_register_variable(
        linked_list, linked_list_variables, linked_list['type'] + 'item')

    input_code  = ''
    output_code = ''
    for obj in objects:
        subs = {
            'type':     linked_list['type'],
            'variable': linked_list['variable'],
            'next':     linked_list['next'],
            'object':   obj,
        }
        input_code  += _substitute(_SS_RESTORE, subs)
        output_code += _substitute(_SS_STORE,  subs)
    return input_code, output_code, linked_list_module(linked_list)


# ---------------------------------------------------------------------------
# allowedParameters
# ---------------------------------------------------------------------------

_AP_ITERATOR = """{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%allowedParameters(allowedParameters,'{$source}',.true.)
   {$type}item => {$type}item%{$next}
end do
"""


def allowed_parameters_linked_list(non_abstract_class, linked_list_variables,
                                   source):
    """Generate the allowed-parameters walk for a non-abstract class's
    linked list.  Returns `(iterator_code, module)` or `('', None)` if the
    class has no linked list.
    """
    if 'linkedList' not in non_abstract_class:
        return '', None
    linked_list = non_abstract_class['linkedList']
    objects = (linked_list.get('object') or '').split()
    linked_list_register_variable(
        linked_list, linked_list_variables, linked_list['type'] + 'item')
    iterator = ''
    for obj in objects:
        iterator += _substitute(_AP_ITERATOR, {
            'type':     linked_list['type'],
            'variable': linked_list['variable'],
            'next':     linked_list['next'],
            'object':   obj,
            'source':   source,
        })
    return iterator, linked_list_module(linked_list)


# ---------------------------------------------------------------------------
# auto-descriptor
# ---------------------------------------------------------------------------

_AD_ITERATOR = """{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%descriptor(parameters)
   {$type}item => {$type}item%{$next}
end do
"""


def auto_descriptor_linked_list(linked_list, linked_list_variables):
    """Generate the auto-descriptor walk.  Returns `(iterator_code, module)`."""
    objects = (linked_list.get('object') or '').split()
    linked_list_register_variable(
        linked_list, linked_list_variables, linked_list['type'] + 'item')
    iterator = ''
    for obj in objects:
        iterator += _substitute(_AD_ITERATOR, {
            'type':     linked_list['type'],
            'variable': linked_list['variable'],
            'next':     linked_list['next'],
            'object':   obj,
        })
    return iterator, linked_list_module(linked_list)


# ---------------------------------------------------------------------------
# assigner
# ---------------------------------------------------------------------------

_ASN_HEAD = """nullify(self%{$variable})
{$type}itemFrom => from%{$variable}
if (associated({$type}itemFrom)) then
   allocate(self%{$variable})
   {$type}itemSelf => self%{$variable}
   do while (associated({$type}itemFrom))
"""

_ASN_PER_OBJECT = """      {$type}itemSelf%{$object} => {$type}itemFrom%{$object}
      call {$type}itemSelf%{$object}%referenceCountIncrement()
"""

_ASN_TAIL = """      {$type}itemFrom => {$type}itemFrom%{$next}
      if (associated({$type}itemFrom)) allocate({$type}itemSelf%{$next})
      {$type}itemSelf => {$type}itemSelf%{$next}
   end do
end if
"""


def assigner_linked_list(linked_list, linked_list_variables):
    """Generate the assignment walk for the linked list.  Returns
    `(iterator_code, module)`.
    """
    objects = (linked_list.get('object') or '').split()
    linked_list_register_variable(
        linked_list, linked_list_variables,
        linked_list['type'] + 'itemSelf',
        linked_list['type'] + 'itemFrom',
    )
    subs_base = {
        'type':     linked_list['type'],
        'variable': linked_list['variable'],
        'next':     linked_list['next'],
    }
    iterator = _substitute(_ASN_HEAD, subs_base)
    for obj in objects:
        iterator += _substitute(_ASN_PER_OBJECT, {**subs_base, 'object': obj})
    iterator += _substitute(_ASN_TAIL, subs_base)
    return iterator, linked_list_module(linked_list)
