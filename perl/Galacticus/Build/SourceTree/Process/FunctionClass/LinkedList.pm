# Contains a Perl module which provides linked list code generation utilities for functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass::LinkedList;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use Exporter 'import';

our @EXPORT_OK = qw(
    deepCopyLinkedList
    stateStoreLinkedList
    allowedParametersLinkedList
    autoDescriptorLinkedList
    assignerLinkedList
    linkedListRegisterVariable
    linkedListModule
);

sub linkedListRegisterVariable {
    # Register pointer variables of the linked list type in the given variable list, unless already present.
    my $linkedList          = shift();
    my $linkedListVariables = shift();
    my @variableNames       = @_;
    push(
	@{$linkedListVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => \@variableNames
	}
	)
	unless ( grep {exists($_->{'type'}) && $_->{'type'} eq $linkedList->{'type'}} @{$linkedListVariables} );
}

sub linkedListModule {
    # Return the module associated with a linked list, or undef if none.
    my $linkedList = shift();
    return exists($linkedList->{'module'}) ? $linkedList->{'module'} : undef();
}

sub deepCopyLinkedList {
    # Create deep-copy instructions for linked list objects.
    my $class                       = shift();
    my $nonAbstractClass            = shift();
    my $linkedListVariables         = shift();
    my $linkedListResetVariables    = shift();
    my $linkedListFinalizeVariables = shift();
    return ("","","",undef())
	unless ( exists($class->{'linkedList'}) );
    my $linkedList = $class->{'linkedList'};
    # Get object names and types.
    my @objects     = split(" ",$linkedList->{'object'    });
    my @objectTypes = split(" ",$linkedList->{'objectType'});
    # Add variables needed for linked list processing.
    &linkedListRegisterVariable($linkedList,$linkedListVariables,
	$linkedList->{'type'}.'item',$linkedList->{'type'}.'destination',$linkedList->{'type'}.'itemNew');
    push(
	@{$linkedListVariables},
	{
	    intrinsic  => 'integer',
	    variables  => [ 'referenceCount___' ]
	}
	)
	unless ( grep {$_->{'variables'}->[0] eq 'referenceCount___'} @{$linkedListVariables} );
    &linkedListRegisterVariable($linkedList,$linkedListResetVariables,   $linkedList->{'type'}.'item');
    &linkedListRegisterVariable($linkedList,$linkedListFinalizeVariables,$linkedList->{'type'}.'item');
    # Generate code for the walk through the linked list.
    my $deepCopyCode;
    my $deepCopyResetCode;
    my $deepCopyFinalizeCode;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type            =                                            $linkedList->{'type'    };
	$code::variable        =                                            $linkedList->{'variable'};
	$code::next            =                                            $linkedList->{'next'    };
	$code::object          =                                            $objects    [$i]         ;
	$code::objectType      =                                            $objectTypes[$i]         ;
	$code::location        = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($class->{'node'},$class->{'node'}->{'line'});
	if ( $i == 0 ) {
	    $deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item             => destination%{$variable}
do while (associated({$type}item))
   ! Undo the reference count increment that resulted from the initial intrinsic assignment.
   referenceCount___={$type}item%{$object}%referenceCountDecrement()
   nullify({$type}item%{$object})
   {$type}itemNew => {$type}item%next
   deallocate({$type}item)
   {$type}item => {$type}itemNew
end do
destination%{$variable} => null            ()
CODE
	}
	$deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}destination      => null            ()
{$type}item             => self%{$variable}
do while (associated({$type}item))
CODE
	if ( $i == 0 ) {
	    $deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
   allocate({$type}itemNew)
   if (associated({$type}destination)) then
      {$type}destination%{$next}     => {$type}itemNew
      {$type}destination             => {$type}itemNew
   else
      destination       %{$variable} => {$type}itemNew
      {$type}destination             => {$type}itemNew
   end if
CODE
	} else {
	    $deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
   if (associated({$type}destination)) then
      {$type}itemNew     => {$type}destination%{$next}
      {$type}destination => {$type}destination%{$next}
   else
      {$type}itemNew     => destination       %{$variable}
      {$type}destination => destination       %{$variable}
   end if
CODE
	}
	$deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
      nullify({$type}itemNew%{$object})
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
CODE
	$deepCopyResetCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%deepCopyReset()
   {$type}item => {$type}item%{$next}
end do
CODE
	$deepCopyFinalizeCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%deepCopyFinalize()
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = &linkedListModule($linkedList);
    return ($deepCopyCode,$deepCopyResetCode,$deepCopyFinalizeCode,$deepCopyModule);
}

sub stateStoreLinkedList {
    # Create state store/restore instructions for linked list objects.
    my $class               = shift();
    my $nonAbstractClass    = shift();
    my $linkedListVariables = shift();
    return ("","","",undef())
	unless ( exists($class->{'linkedList'}) );
    my $linkedList = $class->{'linkedList'};
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    &linkedListRegisterVariable($linkedList,$linkedListVariables,$linkedList->{'type'}.'item');
    # Generate code for the walk through the linked list.
    my $inputCode;
    my $outputCode;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type     = $linkedList->{'type'    };
	$code::variable = $linkedList->{'variable'};
	$code::next     = $linkedList->{'next'    };
	$code::object   = $objects[$i];
	$inputCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%stateRestore(stateFile,gslStateFile,stateOperationID)
   {$type}item => {$type}item%{$next}
end do
CODE
	$outputCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%stateStore(stateFile,gslStateFile,stateOperationID)
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = &linkedListModule($linkedList);
    return ($inputCode,$outputCode,$deepCopyModule);
}

sub allowedParametersLinkedList {
    # Create allowed parameter instructions for linked list objects.
    my $nonAbstractClass    = shift();
    my $linkedListVariables = shift();
    my $source              = shift();
    return ("",undef())
	unless ( exists($nonAbstractClass->{'linkedList'}) );
    my $linkedList = $nonAbstractClass->{'linkedList'};
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    &linkedListRegisterVariable($linkedList,$linkedListVariables,$linkedList->{'type'}.'item');
    # Generate code for the walk through the linked list.
    my $iterator;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type     = $linkedList->{'type'    };
	$code::variable = $linkedList->{'variable'};
	$code::next     = $linkedList->{'next'    };
	$code::object   = $objects[$i];
	$code::source   = $source;
	$iterator .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%allowedParameters(allowedParameters,'{$source}',.true.)
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = &linkedListModule($linkedList);
    return ($iterator,$deepCopyModule);
}

sub autoDescriptorLinkedList {
    # Create auto-descriptor instructions for linked list objects.
    my $linkedList          = shift();
    my $linkedListVariables = shift();
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    &linkedListRegisterVariable($linkedList,$linkedListVariables,$linkedList->{'type'}.'item');
    # Generate code for the walk through the linked list.
    my $iterator;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type     = $linkedList->{'type'    };
	$code::variable = $linkedList->{'variable'};
	$code::next     = $linkedList->{'next'    };
	$code::object   = $objects[$i];
	$iterator .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%descriptor(parameters)
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = &linkedListModule($linkedList);
    return ($iterator,$deepCopyModule);
}

sub assignerLinkedList {
    # Create assignment instructions for linked list objects.
    my $linkedList          = shift();
    my $linkedListVariables = shift();
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    &linkedListRegisterVariable($linkedList,$linkedListVariables,
	$linkedList->{'type'}.'itemSelf',$linkedList->{'type'}.'itemFrom');
    # Generate code for the walk through the linked list.
    my $iterator;
    $code::type      = $linkedList->{'type'    };
    $code::variable  = $linkedList->{'variable'};
    $code::next      = $linkedList->{'next'    };
    $iterator       .= fill_in_string(<<'CODE', PACKAGE => 'code');
nullify(self%{$variable})
{$type}itemFrom => from%{$variable}
if (associated({$type}itemFrom)) then
   allocate(self%{$variable})
   {$type}itemSelf => self%{$variable}
   do while (associated({$type}itemFrom))
CODE
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::object  = $objects[$i];
	$iterator     .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$type}itemSelf%{$object} => {$type}itemFrom%{$object}
      call {$type}itemSelf%{$object}%referenceCountIncrement()
CODE
    }
    $iterator .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$type}itemFrom => {$type}itemFrom%{$next}
      if (associated({$type}itemFrom)) allocate({$type}itemSelf%{$next})
      {$type}itemSelf => {$type}itemSelf%{$next}
   end do
end if
CODE
    my $assignerModule = &linkedListModule($linkedList);
    return ($iterator,$assignerModule);
}

1;
