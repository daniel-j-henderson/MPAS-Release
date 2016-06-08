module kd_tree_mod

    type :: tree_node
        real(kind=RKIND), dimension(:), allocatable :: coords
        type(tree_node), pointer :: left_child => null()
        type(tree_node), pointer :: right_child => null()
        type(tree_node), pointer :: parent => null()
        
        contains
        
    end type

    interface tree_node
        module procedure treeconstructor
    end interface
    
    type :: kd_tree
        type(tree_node), pointer :: root   
        
        contains
        procedure :: build_tree 
    end type kd_tree

    contains

    subroutine build_tree(this, points)
        class(kd_tree) :: this
        real(kind=RKIND), dimension(:,:) :: points

            
    
        this%root => tree_node(points)
    end subroutine build_tree

    function treeconstructor(points)
        class(tree_node), intent(out) :: treeconstructor
        real(kind=RKIND), dimension(:,:) :: points
        integer :: median, n, k
    
        allocate(treeconstructor)
        
        (/k, n/) = shape(points)
        if(n == 1) then
            allocate(treeconstructor%coords(k), source=points(:,1))
            return
        end if    
        median = (n+1) / 2
        allocate(treeconstructor%coords(k), source=points(:,median))

        if(n > 2) treeconstructor%left_child => tree_node(points(:,1:median-1))
        treeconstructor%right_child => tree_node(points(:,median+1,n))
