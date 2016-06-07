module key_mod
  type, abstract :: key
    integer :: id
    contains
      procedure :: setID
      procedure :: getID
      procedure :: deletekey
      procedure(eq), deferred :: equals
      generic :: operator(==) => equals
  end type key

  abstract interface
      function eq(this, other)
        import
        class(key), intent(in) :: this
        class(key), intent(in) :: other
        logical :: eq
      end function eq
  end interface

  contains

    subroutine deletekey(this)
        class(key) :: this
    end subroutine deletekey

    subroutine setID(this, i)
        class(key) :: this
        this%id = i
    end subroutine setID

    function getID(this)
        class(key) :: this
        integer :: getID
        getID = this%id
    end function getID

end module key_mod


module stringkey_mod
  use key_mod
  type, extends (key) :: string_key 
    character(len=100) :: value
    contains
      procedure :: setValue => setStrKeyVal
      procedure :: getValue => getStrKeyVal
      procedure :: setStrKeyEqual
      procedure :: equals => strKeyEquals
      generic :: assignment(=) => setStrKeyEqual
  end type string_key

  interface string_key
    module procedure strconstructor
  end interface

  contains


    function strconstructor(i, val)
        class(string_key), pointer :: strconstructor
        character(len=*) :: val
        integer :: i
        allocate(strconstructor)
        strconstructor%value = trim(val)
        strconstructor%id = i
    end function strconstructor

    subroutine setStrKeyVal(this, val)
      class(string_key) :: this
      character(len=*) :: val
      this%value = val
    end subroutine setStrKeyVal

    function getStrKeyVal(this)
      class(string_key) :: this
      character(len=100) getStrKeyVal
      getStrKeyVal = this%value
    end function getStrKeyVal

    subroutine setStrKeyEqual(this, other)
      class(string_key), intent(out) :: this
      class(string_key), intent(in) :: other
      this%value = other%value
    end subroutine setStrKeyEqual

    function strKeyEquals(this, other)
      class(string_key), intent(in) :: this
      class(key), intent(in) :: other
      logical :: strKeyEquals
      select type(other)
      type is (string_key)
      if (this%value == other%value) then
        strKeyEquals = .true.
      else 
        strKeyEquals = .false.
      endif
      class default
        strKeyEquals = .false.
      end select
    end function strKeyEquals

end module stringkey_mod

module intkey_mod
  use key_mod
  type, extends (key) :: int_key 
    integer :: value
    contains
      procedure :: setValue => setIntKeyVal
      procedure :: getValue => getIntKeyVal
      procedure :: setIntKeyEqual
      procedure :: equals => intKeyEquals
      generic :: assignment(=) => setIntKeyEqual
  end type int_key

  interface int_key
    module procedure intconstructor
  end interface

  contains


    function intconstructor(i, val)
        class(int_key), pointer :: intconstructor
        integer :: val
        integer :: i
        allocate(intconstructor)
        intconstructor%value = val
        intconstructor%id = i
    end function intconstructor

    subroutine setIntKeyVal(this, val)
      class(int_key) :: this
      integer :: val
      this%value = val
    end subroutine setIntKeyVal

    function getIntKeyVal(this)
      class(int_key) :: this
      integer :: getIntKeyVal
      getIntKeyVal = this%value
    end function getIntKeyVal

    subroutine setIntKeyEqual(this, other)
      class(int_key), intent(out) :: this
      class(int_key), intent(in) :: other
      this%value = other%value
    end subroutine setIntKeyEqual

    function intKeyEquals(this, other)
      class(int_key), intent(in) :: this
      class(key), intent(in) :: other
      logical :: intKeyEquals
      select type (other)
      type is (int_key)
      if (this%value == other%value) then
        intKeyEquals = .true.
      else 
        intKeyEquals = .false.
      endif
      class default
      intKeyEquals = .false.
      end select
    end function intKeyEquals

end module intkey_mod


module node_mod
    use key_mod
implicit none
    type node
        class(key), pointer :: mykey => null()
        class(node), pointer :: next => null()
        class(node), pointer :: prev => null()
        contains
        procedure :: getNext
        procedure :: getPrev
        procedure :: setKey
        procedure :: getkey
        procedure :: setNext
        procedure :: setPrev
        procedure :: delete
    end type node

    interface node
        module procedure constructor
    end interface

    contains

    subroutine delete(this)
        class(node) :: this
        !call this%mykey%deletekey()
        if(associated(this%mykey)) deallocate(this%mykey)
        this%next => null()
        this%prev => null()
    end subroutine delete

    function getNext(this)
        class(node) :: this
        class(node), pointer :: getNext
        getNext => this%next
    end function getNext    
    
    function getPrev(this)
        class(node) :: this
        class(node), pointer :: getPrev
        getPrev => this%prev
    end function getPrev    

    subroutine setKey(this, val)
        class(node) :: this
        class(key) :: val
        if (associated(this%mykey)) then
            deallocate(this%mykey)
        end if
        allocate(this%mykey, source=val)
    end subroutine setKey

    function getKey(this)
        class(node) :: this
        class(key), pointer :: getKey
        allocate(getKey, source=this%mykey)
    end function getKey

    subroutine setNext(this, n)
        class(node), pointer :: n
        class(node) :: this
        this%next => n
    end subroutine setNext

    subroutine setPrev(this, n)
        class(node), pointer :: n
        class(node) :: this
        this%prev => n
    end subroutine setPrev

    function constructor(value, n, p)
        class(node), pointer :: constructor
        class(key), pointer :: value
        class(node), pointer :: n, p
        allocate(constructor)
        allocate(constructor%mykey, source=value)
        constructor%next => n
        constructor%prev => p
    end function constructor
end module node_mod



module list_mod
    use node_mod
    use key_mod
    implicit none
    type list
        class(node), pointer :: head => null()
        class(node), pointer :: tail => null()
        integer, public :: s = 0
        contains
        procedure :: insert
        procedure :: isEmpty
        procedure :: contains
        procedure :: size
        procedure :: remove
        procedure :: removeFirst
        procedure :: removeLast
    end type list


    contains


    subroutine insert(this, value)
        class(list) :: this
        class(key), pointer :: value
        class(node), pointer :: new
        if (this%s == 0) then 
            this%head => node(value, this%tail, this%head)
            this%head%prev => null()
            this%head%next => null()
            this%tail => this%head
            this%s = this%s + 1
        else 
            new => null()
            new => node(value, new, this%tail)
            call this%tail%setNext(new)
            this%tail => new
            this%s = this%s + 1
        end if
    end subroutine insert

    function contains(this, t)
        class(list) :: this
        class(key), pointer :: t
        class(node), pointer :: curr
        class(key), pointer :: temp
        logical :: contains
        integer :: i
        contains = .false.
        i = this%s
        curr => this%head
        do while(i > 0)
            temp => curr%getKey()
            if (t == temp) then
                contains = .true.
                exit
            else
                curr => curr%getNext()
                i = i-1
            end if
            deallocate(temp)
        end do
        if (associated(temp)) then
            deallocate(temp)
        end if
    end function contains

    function remove(this, t)
        class(list) :: this
        class(key), pointer :: t
        class(key), pointer :: remove
        class(node), pointer :: curr 
        class(node), pointer :: temp, temp2
        integer :: i
        i = this%s
        curr => this%head 
        do while (i > 0)
            remove => curr%getKey()
            if (i==this%s .and. remove == t) then
                deallocate(remove)
                remove => this%removeFirst()
                exit
            end if
            if (i==1) then
                if(remove == t) then
                deallocate(remove)
                remove => this%removeLast()
                end if
                exit
            end if
            if (remove == t) then
                temp => curr%getPrev()
                temp2 => curr%getNext()
                if(associated(temp)) call temp%setNext(temp2)
                temp => curr%getNext()
                temp2 => curr%getPrev()
                if(associated(temp)) call temp%setPrev(temp2)
                call curr%delete()
                deallocate(curr)
                this%s = this%s - 1
                exit
            else
                deallocate(remove)
                curr => curr%getNext()
                i = i-1
            end if
        end do
    end function remove
        
        

    function isEmpty(this)
        class(list) :: this
        logical :: isEmpty
        if (this%s == 0) then 
            isEmpty = .true.
        else 
            isEmpty = .false.
        end if
    end function isEmpty
    
    function size(this)
        class(list) :: this
        integer :: size
        size = this%s
    end function size

    function removeFirst(this)
        class(list) :: this
        class(key), pointer :: removeFirst
        class(node), pointer :: temp
        removeFirst => this%head%getKey()

        temp => this%head
        this%head => this%head%getNext()
        call temp%delete()
        deallocate(temp)
        this%s = this%s - 1
    end function removeFirst

    function removeLast(this)
        class(list) :: this
        class(key), pointer :: removeLast
        class(node), pointer :: temp
        allocate(removeLast, source=this%tail%getKey())
        temp => this%tail
        this%tail => this%tail%getPrev()
        deallocate(temp)
        this%s = this%s - 1
    end function removeLast

end module list_mod
