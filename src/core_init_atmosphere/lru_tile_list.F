module tilekey_mod
    use key_mod
    integer, parameter :: KKIND=4
    type, extends (key) :: tile_key
        integer :: i, j
        integer :: nx, ny
        real(kind=KKIND), dimension(:,:,:), pointer :: tile_data
        contains
            procedure :: setTile
            procedure :: getData
            procedure :: equals => tileKeyEquals
    end type tile_key

    interface tile_key
        module procedure tileconstructor
    end interface

    contains

        function tileconstructor(i, j, nx, ny, t)
            class(tile_key), pointer :: tileconstructor
            integer :: i, j, nx, ny
            real(kind=KKIND), dimension(:,:,:), pointer, optional :: t

            allocate(tileconstructor)
            tileconstructor%i = i
            tileconstructor%j = j
            tileconstructor%nx = nx
            tileconstructor%ny = ny
            if(present(t)) tileconstructor%tile_data => t

        end function tileconstructor

        subroutine setTile(this, i, j, nx, ny, t)
            class(tile_key) :: this
            integer :: i, j, nx, ny
            real(kind=KKIND), dimension(:,:,:), pointer :: t

            this%i = i
            this%j = j
            this%nx = nx
            this%ny = ny
            this%tile_data => t
        end subroutine setTile

        function getData(this)
            class(tile_key) :: this
            real(kind=KKIND), dimension(:,:,:), pointer :: getData

            getData => this%tile_data
        end function

        function tileKeyEquals(this, other)
            class(tile_key), intent(in) :: this
            class(key), intent(in) :: other
            logical :: tileKeyEquals

            select type(other)
            type is (tile_key)
                tileKeyEquals = .true.
                if (this%i .ne. other%i) tileKeyEquals = .false.
                if (this%j .ne. other%j) tileKeyEquals = .false.
                if (this%nx .ne. other%nx) tileKeyEquals = .false.
                if (this%ny .ne. other%ny) tileKeyEquals = .false.
            class default
                tileKeyEquals = .false.
            end select
        end function tileKeyEquals

end module tilekey_mod


module lrulist_mod
    use list_mod
    use node_mod
    use key_mod
    use tilekey_mod
    implicit none

    type, extends(list) :: lru_list
        integer, public :: capacity
        integer :: label
        contains
        procedure :: lru_insert
        procedure :: hit
        procedure :: lru_remove
        procedure :: get
    end type lru_list

    contains

    subroutine lru_insert(this, value)
        class(lru_list), intent(inout) :: this
        class(key), pointer, intent(in) :: value
        class(node), pointer :: new
        if (this%s == this%capacity) then
            call this%lru_remove()
        end if
        if (this%s > this%capacity) print *, "Size error"
        if (this%s == 0) then
            call value%setID(1)
            call this%insert(value)
            this%label = 1
        else
            call value%setID(this%label+1)
            this%label = this%label+1
            call this%insert(value)
        end if
    end subroutine lru_insert
function hit(this, value)
        class(lru_list) :: this
        class(key) :: value
        class(node), pointer :: curr
        class(key), pointer :: temp
        logical :: hit
        integer :: i, x, y

        i=this%s
        hit = .false.
        curr => this%head
        do while(i > 0)
            temp => curr%getKey()
            if (value == temp) then
                hit = .true.
                call temp%setID(this%label+1)
                this%label = this%label+1
                call curr%setKey(temp)
                deallocate(temp)
                temp => curr%getKey()
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
    end function hit
                       
    subroutine lru_remove(this)
        class(lru_list) :: this
        integer :: x, y, i
        class(key), pointer :: del, temp
        class(node), pointer :: curr
        real(kind=8), dimension(:,:,:), pointer :: arr

        curr => this%head
        del => curr%getKey()
        temp => del
        i = this%s - 1
        do while(i > 0)
            curr => curr%getNext()
            temp => curr%getKey()
            if (temp%getID() < del%getID()) then
                deallocate(del)
                del => temp
            else
                if(associated(temp)) deallocate(temp)
            end if
            i = i-1
        end do
        temp => this%remove(del)
        deallocate(temp)
        deallocate(del)
    end subroutine lru_remove

    function get(this, value)
        class(lru_list) :: this
        class(key), pointer :: get, value
        integer :: i
        class(node), pointer :: curr

        curr => this%head
        i=this%s
        do while (i > 0)
            get => curr%getKey()
            if (get == value) then
                exit
            else
                curr => curr%getNext()
                deallocate(get)
            end if
            i = i-1
        end do
    end function get

end module lrulist_mod

