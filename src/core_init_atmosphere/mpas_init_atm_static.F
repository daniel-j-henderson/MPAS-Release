! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_init_atm_static
!==================================================================================================
 use atm_advection
! use mpas_configure
 use mpas_dmpar
 use mpas_pool_routines
 use init_atm_hinterp
 use init_atm_llxy

 use mpas_atmphys_utilities

 implicit none
 private :: foo, foo2, get_tile, con_lx, isClosest, div_and_conq
 public:: init_atm_static,           &
          init_atm_static_orogwd,    &
          init_atm_check_read_error, &
          nearest_cell,              &
          sphere_distance               

 contains



 function get_tile(i, j, cache, nx, ny, nz, nh, data_p, sub_p, isigned, endian, scalefactor, wordsize)
    use tilekey_mod
    use lrulist_mod
    implicit none
    real(kind=4), dimension(:,:,:), pointer :: get_tile
    integer, intent(in) :: i, j, nx, ny, nh, nz
    type(lru_list), intent(inout) :: cache
    character(len=StrKIND+1), intent(in) :: data_p, sub_p
    real(kind=4), intent(in) :: scalefactor
    integer, intent(in) :: endian,isigned,wordsize

    character(len=StrKIND) :: fname
    class(key), pointer :: keyptr, tptr
    integer :: istatus, u, l
    
    keyptr => tile_key(i, j, nx, ny)
    if (cache%hit(keyptr)) then
        tptr => cache%get(keyptr)
        select type (tptr)
        type is (tile_key)
            get_tile => tptr%getData()
        class default
            write(0,*), "error, wrong key type"
        end select
    else
        write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(data_p)// &
                     trim(sub_p),mod(i,43200),'-',mod(i+nx-1-1,43200)+1,'.',j,'-',j+ny-1
        write(0,*) trim(fname)
        allocate(get_tile(nx+2*nh, ny+2*nh, 1))
        call read_geogrid(fname, len_trim(fname), get_tile, nx+2*nh, ny+2*nh, nz, isigned, endian, scalefactor, wordsize, istatus)
        call init_atm_check_read_error(istatus, fname)
        write (0,*) "Read success"
        tptr => tile_key(i, j, nx, ny, get_tile)
        call cache%lru_insert(tptr)
    end if

    deallocate(tptr)
    deallocate(keyptr)

 end function get_tile


 function isClosest(lat_pt, lon_pt, latCell, lonCell, nEoC, CoC, iCell)
    logical :: isClosest
    real(kind=RKIND), intent(in) :: lat_pt, lon_pt
    real(kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell
    integer, intent(in) :: nEoC, iCell
    integer, dimension(:), intent(in) :: CoC

    real(kind=RKIND) :: dist
    integer :: i


    isClosest = .true.

    dist = sphere_distance(lat_pt, lon_pt, latCell(iCell), lonCell(iCell), 100.0)
    do i=1,nEoC
        if (sphere_distance(lat_pt, lon_pt, latCell(CoC(i)), lonCell(CoC(i)), 100.0) < dist) then
            isClosest = .false.
            exit
        end if
    end do
 end function isClosest

 recursive subroutine div_and_conq(mask, istart, iend, jstart, jend, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
    implicit none
    logical, dimension(:,:,:), intent(out) :: mask
    integer, intent(in) :: istart, iend, jstart, jend, ioffset, joffset, nEoC, iCell
    integer, dimension(:), intent(in) :: CoC
    real(kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell

    real(kind=RKIND) :: lat_pt, lon_pt, lat_cell, lon_cell, x_cell, y_cell, z_cell, dist, tdim=0.0083333333, r=100.0
    real(kind=RKIND), dimension(:,:,:), allocatable :: lat_arr, lon_arr, dist2_icell_arr, dist2_ocell_arr
    logical :: closest = .true.
    integer :: isplit, jsplit, i, j, l, nz=1

    closest = .true.
    lat_pt = -89.99583 + jstart * tdim
    lat_pt = lat_pt * PI / 180.0
    lon_pt = -179.99583 + istart * tdim
    lon_pt = lon_pt * PI / 180.0
    if (.not. isClosest(lat_pt, lon_pt, latCell, lonCell, nEoC, CoC, iCell)) then
        closest = .false. 
    end if

    lat_pt = -89.99583 + jend * tdim
    lat_pt = lat_pt * PI / 180.0
    lon_pt = -179.99583 + istart * tdim
    lon_pt = lon_pt * PI / 180.0
    if (closest .and. .not. isClosest(lat_pt, lon_pt, latCell, lonCell, nEoC, CoC, iCell)) then
        closest = .false. 
    end if

    lat_pt = -89.99583 + jstart * tdim
    lat_pt = lat_pt * PI / 180.0
    lon_pt = -179.99583 + iend * tdim
    lon_pt = lon_pt * PI / 180.0
    if (closest .and. .not. isClosest(lat_pt, lon_pt, latCell, lonCell, nEoC, CoC, iCell)) then
        closest = .false. 
    end if

    lat_pt = -89.99583 + jend * tdim
    lat_pt = lat_pt * PI / 180.0
    lon_pt = -179.99583 + iend * tdim
    lon_pt = lon_pt * PI / 180.0
    if (closest .and. .not. isClosest(lat_pt, lon_pt, latCell, lonCell, nEoC, CoC, iCell)) then
        closest = .false. 
    end if

    if (closest) then
!        mask((istart-ioffset):(iend-ioffset), (jstart-joffset):(jend-joffset),:) = .true.
        return
    else
        if (iend - istart < 5 .or. jend - jstart < 5) then
            if(iend - istart > 5) then
                isplit = (iend - istart) / 2 + istart
                call div_and_conq(mask, istart, isplit, jstart, jend, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
                call div_and_conq(mask, isplit+1, iend, jstart, jend, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
            else if(jend - jstart > 5) then
                jsplit = (jend - jstart) / 2 + jstart
                call div_and_conq(mask, istart, iend, jsplit+1, jend, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
                call div_and_conq(mask, istart, iend, jstart, jsplit, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
            else
                do j=jstart, jend
                lat_pt = -89.99583 + (j-1) * tdim
                lat_pt = lat_pt * PI / 180.0
                do i=istart, iend
                    lon_pt = -179.99583 + (i-1) * tdim
                    lon_pt = lon_pt * PI / 180.0
                    if(.not. isClosest(lat_pt, lon_pt, latCell, lonCell, nEoC, CoC, iCell)) mask(i-ioffset,j-joffset,:) = .false.
                end do
                end do
                return
            end if
        else
            isplit = (iend - istart) / 2 + istart
            jsplit = (jend - jstart) / 2 + jstart
            call div_and_conq(mask, istart, isplit, jstart, jsplit, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
            call div_and_conq(mask, istart, isplit, jsplit+1, jend, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
            call div_and_conq(mask, isplit+1, iend, jstart, jsplit, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
            call div_and_conq(mask, isplit+1, iend, jsplit+1, jend, ioffset, joffset, nEoC, CoC, latCell, lonCell, iCell)
        end if
    end if

 end subroutine div_and_conq

 subroutine foo2(mesh, dims, configs, geog_data_path, geog_sub_path, nx, ny, nz, nh, &
                 field, ncats, isigned, endian, scalefactor, wordsize)
    
 use intkey_mod
 use list_mod
 use kd_tree_mod
 type (mpas_pool_type), intent(inout) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs
 integer, pointer :: nCells, nEdges, nVertices, nCellsSolve, vertexDegree
 
 character(len=StrKIND+1), intent(in) :: geog_data_path, geog_sub_path
 character(len=StrKIND) :: fname
 integer, intent(in) :: endian,isigned,wordsize
 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell, cellsOnVertex, verticesOnCell, edgesOnVertex
 real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
 real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
 real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
 real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
 real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
 real (kind=RKIND), pointer :: radius
 real(kind=4), intent(in) :: scalefactor
 integer, dimension(:), allocatable :: borderVertices
 type(list) ::bvertslist
 class(key), pointer :: keyptr
 type(kd_tree) :: cell_tree
 real(kind=RKIND), dimension(:,:), allocatable :: points
 real (kind=RKIND), dimension(:), pointer, intent(out) :: field
 real(kind=4), dimension(:,:,:), allocatable :: rarray
 integer, dimension(:,:), allocatable :: num
 integer :: nBorderVertices, temp
 real(kind=RKIND), dimension(2) :: c_lats, c_lons
 real(kind=RKIND), dimension(3) :: pt 

 integer :: i, j, k, v, l, it, jt, nx, ny, nz, nh, nti, ntj, jts, ncats, istatus, iCheckStart, iCheckEnd, jCheckStart, jCheckEnd, &
            iTileStart, iTileEnd, jTileStart, jTileEnd, iRegionStart, iRegionEnd, jRegionStart, jRegionEnd 
 real(kind=RKIND) :: lat_pt, lon_pt, lat_min, lon_min, lat_max, lon_max, lat_cell, lon_cell, tdim=.00833333
 

 call mpas_pool_get_array(mesh, 'xCell', xCell)
 call mpas_pool_get_array(mesh, 'yCell', yCell)
 call mpas_pool_get_array(mesh, 'zCell', zCell)
 call mpas_pool_get_array(mesh, 'latCell', latCell)
 call mpas_pool_get_array(mesh, 'lonCell', lonCell)
 call mpas_pool_get_array(mesh, 'latEdge', latEdge)
 call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
 call mpas_pool_get_array(mesh, 'latVertex', latVertex)
 call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
 
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
 call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
 call mpas_pool_get_config(mesh, 'sphere_radius', radius)

 call mpas_pool_get_dimension(dims, 'nCells', nCells)
 call mpas_pool_get_dimension(dims, 'nEdges', nEdges)
 call mpas_pool_get_dimension(dims, 'nVertices', nVertices)
 call mpas_pool_get_dimension(dims, 'vertexDegere', vertexDegree)
 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve) 

 write(0,*) "Begin foo2(). nCells=", nCells, ", nCellsSolve=", nCellsSolve

 allocate(points(3,nCells))
 allocate(field(nCells))
 allocate(num(0:ncats, nCells))
 allocate(rarray(nx+2*nh, ny+2*nh, nz))
 do i=1, nCells-1
    points(:,i) = (/xCell(i), yCell(i), zCell(i)/)
 end do
 call cell_tree%create_tree(points)
 write (0,*) "Created Tree of Cells"

 do i=nCellsSolve+1, nCells
 do j=1,nEdgesOnCell(i)
 do k=1,vertexDegree
    temp = cellsOnVertex(k, verticesOnCell(j,i))
    if(temp <= nCellsSolve) then
        keyptr => int_key(0,temp)
        if (.not. bvertslist%contains(keyptr)) then
            call bvertslist%insert(keyptr)
            lon_cell = lonCell(i)
            lat_cell = latCell(i)
            if (lon_cell > lon_max) then
                lon_max = lon_cell
            else if(lon_cell < lon_min) then
                lon_min = lon_cell
            end if

            if (lat_cell > lat_max) then
                lat_max = lat_cell
            else if (lat_cell < lat_min) then   
                lat_min = lat_cell
            end if

            exit
        end if
        deallocate(keyptr)
    end if
 end do
 end do
 end do

 nBorderVertices = bvertslist%size()
 allocate(borderVertices(nBorderVertices))
 do i=1,nBorderVertices
     keyptr => bvertslist%removeFirst()
    select type(keyptr)
    type is (int_key)
    borderVertices(i) = keyptr%getValue()
    class default
    end select 
    deallocate (keyptr)
 end do
 write (0,*) "Gathered ", nBorderVertices, " border vertices"

 lat_min = lat_min * 180.0 / PI
 lat_max = lat_max * 180.0 / PI
 lon_min = lon_min * 180.0 / PI
 lon_max = lon_max * 180.0 / PI
 iRegionStart = floor((lon_min + 180.0) / tdim)
 iRegionEnd = ceiling((lon_max + 180.0) / tdim)
 jRegionStart = floor((lat_min + 90.0) / tdim)
 jRegionEnd = ceiling((lat_max + 90.0) / tdim)
 ! These are the i, j bounding values in the space of the geographical
 ! dataset. I expand the bounding box by an arbitrary amount on all sides
 ! to avoid cutting off any necessary data. This amount may need
 ! calibrated.

! Check to see if the pole is in the region. If so, expand region start/end to
! the pole.
! polecap=.false.
! if (isClosest(PI/2, 1.0, latCell, lonCell, nEoC, CoC, iCell)) then
!     lat_max = 180.0
!     jCheckEnd = 21600
!     polecap=.true.
!     write (0,*) "North Pole"
! end if
! if (isClosest(PI/(-2), 1.0, latCell, lonCell, nEoC, CoC, iCell)) then
!     jCheckStart = 1
!     lat_min = -180.0
!     polecap = .true.
!     write (0,*) "South Pole"
! end if
! if (polecap) then
!     iCheckstart=1
!     iCheckEnd=43200
! end if

 ! This just makes finding nti, ntj simpler
 if (mod(iCheckStart, nx) == 0) iRegionStart = iRegionStart - 1 
 if (mod(jCheckStart, ny) == 0) jRegionStart = jRegionStart - 1

 nti = (iRegionEnd - 1) / nx - iRegionStart / nx + 1
 ntj = (jRegionEnd - 1) / ny - jRegionStart / ny + 1
 iTileStart = iRegionStart - mod(iRegionStart, nx) + 1 
 jTileStart = jRegionStart - mod(jRegionStart, ny) + 1 

 lat_min = lat_min * PI / 180.0
 lat_max = lat_max * PI / 180.0

 lon_min = lon_min * PI / 180.0
 lon_max = lon_max * PI / 180.0

 write (0,*) "lat min, max:", lat_min, lat_max
 write (0,*) "lon min, max:", lon_min, lon_max
 write (0,*) "number of tiles to check:", nti*ntj
 write (0,*) "iRegionStart, End:", iRegionStart, iRegionEnd
 write (0,*) "jRegionStart, End:", jRegionStart, jRegionEnd

 jts = jTileStart

 do it=1,nti
 jTileStart = jts
 do jt=1,ntj
    c_lats(1) = -89.99583 + (jTileStart + (jt-1)*ny) * tdim
    c_lats(2) = -89.99583 + (jTileStart + jt*ny - 1) * tdim
    c_lons(1) = -179.99583 + (iTileStart + (it-1)*nx) * tdim
    c_lons(2) = -179.99583 + (iTileStart + it*nx-1) * tdim
    c_lats = c_lats * PI / 180.0
    c_lons = c_lons * PI / 180.0
    temp = nInsideCorners(c_lats, c_lons, radius, nCellsSolve, cell_tree) !other stuff, needs to be able to
                                          !determine nearest cell to each corner
    if(temp == 4) then
        iCheckStart = 1
        iCheckEnd = nx
        jCheckStart = 1
        jCheckEnd = ny
    else if(temp == 0) then
        lat_min = 100.0
        lat_max = -100.0
        lon_max = 100.0
        lon_min = -100.0
        do v=1,nBorderVertices
            if(latVertex(borderVertices(v)) <= c_lats(2) .and. latVertex(borderVertices(v)) >= c_lats(1) .and.&
               lonVertex(borderVertices(v)) <= c_lons(2) .and. lonVertex(borderVertices(v)) >= c_lons(1)) then
                
                iCheckStart = 1
                iCheckEnd = nx
                jCheckStart = 1
                jCheckEnd = ny
                exit
            end if
        end do
    else
        cycle
    end if

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
                 trim(geog_sub_path),mod(iTileStart,43200),'-',mod(iTileStart+nx-1-1,43200)+1,'.',jTileStart,'-',jTileStart+ny-1
    write(0,*) trim(fname)
    call read_geogrid(fname, len_trim(fname), rarray, nx+2*nh, ny+2*nh, nz, isigned, endian, scalefactor, wordsize, istatus)
    call init_atm_check_read_error(istatus, fname)
    
    do j=jCheckStart, jCheckEnd
    do i=iCheckStart, iCheckEnd
        lat_pt = -89.99583 + (jTileStart + j - 1) * tdim
        lon_pt = -179.99583 + (iTileStart + i - 1) * tdim
        lat_pt = lat_pt * PI / 180.0
        lon_pt = lon_pt * PI / 180.0
        call con_lx(lat_pt, lon_pt, radius, pt(1), pt(2), pt(3)) 
        l = cell_tree%nearest_cell(pt)
        if (l <= nCellsSolve) then
            if (ncats == 1) then
                field(l) = field(l) + rarray(i, j, 1)
                num(1, l) = num(1, l) + 1
            else
                num(rarray(i, j, 1), l) = num(rarray(i, j, 1), l) + 1
            end if
        end if
    end do
    end do

 jTileStart = jTileStart + ny
 end do
 iTileStart = iTileStart + nx
 end do

 if (ncats == 1) then
    field(:) = field(:) / num(1,:)
 else
    do i=1,nCellsSolve
        field(i) = 0
        do j=1,ncats
            if (num(j, i) > num(int(field(i)), i)) then
                field(i) = j
            end if
        end do
    end do
 end if

 write (0,*) "End subroutine foo2"
 !TODO//Deallocate stuff
 end subroutine foo2

 function nInsideCorners(c_lats, c_lons, radius, nCellsSolve, tree)
    use kd_tree_mod
    real(kind=RKIND), dimension(2), intent(in) :: c_lats, c_lons
    class(kd_tree) :: tree
    real(kind=RKIND), intent(in) :: radius
    integer, intent(in) :: nCellsSolve
    integer :: nInsideCorners

    real(kind=RKIND), dimension(3) :: pt
   
    nInsideCorners = 0 
    call con_lx(c_lats(1), c_lons(1), radius, pt(1), pt(2), pt(3))
    if (tree%nearest_cell(pt) <= nCellsSolve) nInsideCorners = nInsideCorners+1
    call con_lx(c_lats(1), c_lons(2), radius, pt(1), pt(2), pt(3))
    if (tree%nearest_cell(pt) <= nCellsSolve) nInsideCorners = nInsideCorners+1
    call con_lx(c_lats(2), c_lons(1), radius, pt(1), pt(2), pt(3))
    if (tree%nearest_cell(pt) <= nCellsSolve) nInsideCorners = nInsideCorners+1
    call con_lx(c_lats(2), c_lons(2), radius, pt(1), pt(2), pt(3))
    if (tree%nearest_cell(pt) <= nCellsSolve) nInsideCorners = nInsideCorners+1

 end function nInsideCorners

 

 subroutine foo(geog_data_path, geog_sub_path, nx, ny, nz, nh, nCells, nCellsSolve, nVertices, latCell, lonCell, &
                latVertex, lonVertex, verticesOnCell, cellsOnCell, nEdgesOnCell, field, nCats, isigned, endian, scalefactor, wordsize)
    use IEEE_ARITHMETIC
    use tilekey_mod
    use lrulist_mod

    implicit none
 
    character(len=StrKIND+1), intent(in) :: geog_data_path, geog_sub_path
    character(len=StrKIND) :: fname
    integer, intent(in) :: nx, ny, nz, nh, nCells, nVertices, nCats, nCellsSolve
    integer, intent(in) :: endian,isigned,wordsize
    integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
    integer, dimension(:,:), pointer, intent(in) :: cellsOnCell, verticesOnCell
    real (kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell
    real (kind=RKIND), dimension(:), pointer, intent(in) :: latVertex, lonVertex
    real (kind=RKIND), dimension(:), pointer, intent(out) :: field

    real(kind=4), intent(in) :: scalefactor
    real(kind=4), dimension(:,:,:), allocatable :: rarray
    real(kind=RKIND), dimension(:,:,:), allocatable :: lat_arr, lon_arr, dist2_icell_arr, dist2_ocell_arr
    real(kind=4), dimension(:,:,:), pointer :: temparray
    integer, dimension(:,:,:), allocatable :: int_arr
    integer :: iTileStart, jTileStart, iTileEnd, jTileEnd, iCheckStart, iCheckEnd, jCheckStart, jCheckEnd,&
                i, j, l, u, iCell, nti, ntj, istatus, stride, it, jt, x, x2, y, y2, tempint, nEoC
    real (kind=RKIND) :: lat_pt, lon_pt, lat_min, lat_max, lon_min, lon_max, tdim, tile_w, tile_h, &
                         dist, lat_cell, lon_cell, x_cell, y_cell, z_cell, r 
    integer, dimension(:), allocatable :: num, CoC
    logical :: closest, polecap
    logical, dimension(:,:,:), allocatable :: mask_arr

    real(kind=8) :: extime
    integer(kind=8) :: t1, t2, t3, t4, t5, t6, cntrate 
    integer, parameter :: cacheSize = 5
    
    type(lru_list) :: tilecache
    character :: mask_method = 'A'
    integer :: solve_method = 1

    namelist /solve/ solve_method, mask_method
    open(41,file='namelist.daniel')
    read(41,solve)
    close(41)

    r = 100.0
    extime = 0.0
    tilecache%tinsert = 0.0
    tilecache%thit = 0.0
    tilecache%tremove = 0.0
    tilecache%tget = 0.0

    call system_clock(t1, cntrate)

    allocate(field(nCells))
    allocate(num(0:nCats))
    field(:) = 0.0
    num(:) = 0
    tdim = 0.0083333333
    tile_w = nx * tdim
    tile_h = ny * tdim
    tilecache%capacity = 60
    

    MAIN: do iCell=1, nCellsSolve
        nEoC = nEdgesOnCell(iCell)
        allocate(CoC(nEoC))
        CoC(:) = cellsOnCell(1:nEoC, iCell)
        call system_clock(t2)
        lat_min = latVertex(verticesOnCell(1, iCell))
        lat_max = latVertex(verticesOnCell(1, iCell))
        lon_max = lonVertex(verticesOnCell(1, iCell))
        lon_min = lonVertex(verticesOnCell(1, iCell))
        ! Find the bounding box lat/lon
        EXPAND: do j=2, nEdgesOnCell(iCell)
            l = verticesOnCell(j, iCell)
            if (latVertex(l) > lat_max) then
                lat_max = latVertex(l)
            else if (latVertex(l) < lat_min) then
                lat_min = latVertex(l)
            end if

            if (lonVertex(l) > lon_max) then
                lon_max = lonVertex(l)
            else if (lonVertex(l) < lon_min) then
                lon_min = lonVertex(l)
            end if
        end do EXPAND
        ! Check to see if the cell crosses the date line. If so, shift the
        ! western half to the east and re-make the bounding box lat lon
        if ((lon_max - lon_min .gt. PI) .and. ((abs(lonCell(iCell) - lon_max) .lt. PI/6.0) .or. (abs(lon_min - lonCell(iCell)) .lt. PI/6.0))) then
            lon_max = lonVertex(verticesOnCell(1, iCell))
            lon_min = lonVertex(verticesOnCell(1, iCell))
            if (lon_min < PI) then
                lon_min = lon_min + 2*PI
                lon_max = lon_max + 2*PI
            end if
            EXPANDAGAIN: do j=2, nEdgesOnCell(iCell)
                l = verticesOnCell(j, iCell)
                lon_pt = lonVertex(l)
                if (lon_pt < PI) then 
                    lon_pt = lon_pt + 2*PI
                end if
                if (lon_pt > lon_max) then
                    lon_max = lon_pt
                else if (lon_pt < lon_min) then
                    lon_min = lon_pt
                end if
            end do EXPANDAGAIN
        end if
        lat_min = lat_min * 180.0 / PI
        lat_max = lat_max * 180.0 / PI
        lon_min = lon_min * 180.0 / PI
        lon_max = lon_max * 180.0 / PI
        iCheckStart = floor((lon_min + 180.0) / tdim)
        iCheckEnd = ceiling((lon_max + 180.0) / tdim)
        jCheckStart = floor((lat_min + 90.0) / tdim)
        jCheckEnd = ceiling((lat_max + 90.0) / tdim)
        ! These are the i, j bounding values in the space of the geographical
        ! dataset. I expand the bounding box by an arbitrary amount on all sides
        ! to avoid cutting off any necessary data. This amount may need
        ! calibrated.
        iCheckStart = max(iCheckStart - 30, 1)
        jCheckStart = max(jCheckStart - 30, 1)
        iCheckEnd = iCheckEnd + 30
        jCheckEnd = min(jCheckEnd + 30, 21600)

        ! At the poles, artificially extend the j part of the bounding box as
        ! high/low as the pole.
        polecap=.false.
        if (isClosest(PI/2, 1.0, latCell, lonCell, nEoC, CoC, iCell)) then
            lat_max = 180.0
            jCheckEnd = 21600
            polecap=.true.
            write (0,*) "North Pole"
        end if
        if (isClosest(PI/(-2), 1.0, latCell, lonCell, nEoC, CoC, iCell)) then
            jCheckStart = 1
            lat_min = -180.0
            polecap = .true.
            write (0,*) "South Pole"
        end if
        if (polecap) then
            iCheckstart=1
            iCheckEnd=43200
        end if

        ! This just makes finding nti, ntj simpler
        if (mod(iCheckStart, nx) == 0) iCheckStart = iCheckStart - 1 
        if (mod(jCheckStart, ny) == 0) jCheckStart = jCheckStart - 1

        nti = (iCheckEnd - 1) / nx - iCheckStart / nx + 1
        ntj = (jCheckEnd - 1) / ny - jCheckStart / ny + 1
        iTileStart = iCheckStart - mod(iCheckStart, nx) + 1 
        jTileStart = jCheckStart - mod(jCheckStart, ny) + 1 
        lat_cell = latCell(iCell)
        lon_cell = lonCell(iCell)

        lat_min = lat_min * PI / 180.0
        lat_max = lat_max * PI / 180.0
        lon_min = lon_min * PI / 180.0
        lon_max = lon_max * PI / 180.0
        num(:) = 0

        select case (solve_method)
    
        case (1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       SOLVE METHOD 1      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! This method first creates a mask over rarray (the bounding box) then
        ! uses that to put the data into the field(iCell). There are 3 ways to
        ! make this mask. 

        allocate(rarray(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))
        allocate(mask_arr(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))
        allocate(int_arr(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz), source=0)
        mask_arr(:,:,:) = .true.
        rarray(:,:,:) = 0.0

        ! Gather rarray, the geographical data within the bounding box.

        call system_clock(t4)
        ICHECK: do it=1,nti
        jTileStart = jCheckStart - mod(jCheckStart, ny) + 1
        JCHECK: do jt=1,ntj
            !iTileStart = mod(iTileStart, 43200)
            !jTileStart = mod(jTileStart, 21600)
            temparray => get_tile(mod(iTileStart, 43200), mod(jTileStart, 21600), tilecache, nx, ny, nz, nh, geog_data_path, geog_sub_path, isigned, endian, scalefactor, wordsize)
            iTileEnd = iTileStart + nx - 1
            jTileEnd = jTileStart + ny - 1
            x = max(iCheckStart, iTileStart)
            x2 = min(iCheckEnd, iTileEnd)
            y = max(jCheckStart, jTileStart)
            y2 = min(jCheckEnd, jTileEnd)

            rarray(x:x2, y:y2, :) = temparray(x - iTileStart + 1 + nh:x2 - iTileStart + 1 + nh,y - jTileStart + 1 + nh:y2 - jTileStart + 1 + nh,:) 
            jTileStart = jTileStart+ny
        end do JCHECK
        iTileStart = iTileStart + nx
        end do ICHECK

        

        call system_clock(t3)
!        write(0,*) "Time to gather rarray", real(t3-t4) / real(cntrate)

        call system_clock(t5)

        select case (mask_method)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       MASK METHOD A       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! This is a recursive method. It divides the mask into quadrants. If all
        ! 4 corners of the quadrant are within iCell, the whole quadrant is kept
        ! as .true. At the teminating quadrant, each point is checked
        ! individually. Fastest but least correct method (up to 25m max error on
        ! very coarse meshes).
        case ('A')

        call system_clock(t3)
        call div_and_conq(mask_arr, iCheckStart, iCheckEnd, jCheckStart, jCheckEnd, iCheckStart-1, jCheckStart-1, nEoC, CoC, latCell, lonCell, iCell)
        call system_clock(t4)
!        write(0,*) "Time D+C: ", real(t4-t3) / real(cntrate)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       MASK METHOD B       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! This is an array-based method. Instead of checking each point
        ! individually, perform array operations and compare the distance
        ! from iCell to each neighbor to create the mask. Slower than A.

        case ('B')

        allocate(lat_arr(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))
        allocate(lon_arr(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))
        allocate(dist2_icell_arr(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))
        allocate(dist2_ocell_arr(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))
        call con_lx(lat_cell, lon_cell, r, x_cell, y_cell, z_cell)

        do i=iCheckStart,iCheckEnd
            lon_pt = -179.99583 + (i-1) * tdim
            lon_pt = lon_pt * PI / 180.0

            lon_arr(i,:,:) = lon_pt
        end do
        
        do j=jCheckStart,jCheckEnd
            lat_pt = -89.99583 + (j-1) * tdim
            lat_pt = lat_pt * PI / 180.0         
            
            lat_arr(:,j,:) = lat_pt
        end do

        call arrDist2(lat_arr, lon_arr, r, x_cell, y_cell, z_cell, dist2_icell_arr)

        do l=1,nEdgesOnCell(iCell)
            call con_lx(latCell(cellsOnCell(l,iCell)), lonCell(cellsOnCell(l,iCell)), r, x_cell, y_cell, z_cell)
            call arrDist2(lat_arr, lon_arr, r, x_cell, y_cell, z_cell, dist2_ocell_arr)
            where(dist2_ocell_arr < dist2_icell_arr) mask_arr = .false.
        end do

        deallocate(lat_arr, lon_arr, dist2_icell_arr, dist2_ocell_arr) 


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       MASK METHOD C       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Point by point, determine if it is owned by iCell and if not, set to
        ! false. Slower than B.        
        case ('C')

        do j=jCheckStart, jCheckEnd
        lat_pt = -89.99583 + (j-1) * tdim
        lat_pt = lat_pt * PI / 180.0
        do i=iCheckStart, iCheckEnd
            lon_pt = -179.99583 + (i-1) * tdim
            lon_pt = lon_pt * PI / 180.0

            dist = sphere_distance(lat_pt, lon_pt, lat_cell, lon_cell, 1.0)
            closest = .true.
            DISTC: do l=1, nEdgesOnCell(iCell)
                if (sphere_distance(lat_pt, lon_pt, latCell(cellsOnCell(l, iCell)), lonCell(cellsOnCell(l, iCell)), 1.0) &
                    < dist) then
                    closest = .false.
                    exit
                end if
            end do DISTC
            
            if (.not. closest) mask_arr(i,j,:) = .false.
        end do
        end do

        case default
        end select

        !MASK MADE


        call system_clock(t3)

        do j=jCheckStart, jCheckEnd
        do i=iCheckStart, iCheckEnd
            if (mask_arr(i,j,1) .and. ncats == 1) then
                field(iCell) = field(iCell) + rarray(i,j,1)
                num(1) = num(1) + 1
            else if (mask_arr(i,j,1) .and. ncats > 1) then
                num(int(rarray(i,j,1))) = num(int(rarray(i,j,1))) + 1
            end if
        end do
        end do

! Alternate version of above, produces error at the south pole for some reason.

!        if (ncats == 1) then
!            field(iCell) = sum(rarray(:,:,1), MASK=mask_arr(:,:,1))
!            where(mask_arr(:,:,1)) 
!                int_arr(:,:,1) = 1
!            elsewhere
!                int_arr(:,:,1) = 0
!            end where
!            num(1) = sum(int_arr(:,:,1))        
!        else
!            do l=0,ncats
!                where(mask_arr(:,:,1) .and. int(rarray(:,:,1)) == l) 
!                    int_arr(:,:,1) = 1
!                elsewhere
!                    int_arr(:,:,1) = 0
!                end where
!                num(l) = sum(int_arr(:,:,1))
!            end do
!        end if
    

        call system_clock(t6)
        write(0,*) "Time to make mask iCell ", iCell, ": ", real(t3-t5) / real(cntrate)

        
        deallocate(rarray, mask_arr, int_arr, CoC)





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       SOLVE METHOD 2      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! For each tile, get it using get_tile and process the relevant points

        case (2)

        call system_clock(t5)
        lat_cell = latCell(iCell)
        lon_cell = lonCell(iCell)
        ITCHECK: do it=1,nti
        j = jCheckStart
        jTileStart = jCheckStart - mod(jCheckStart, ny) + 1 
        JTCHECK: do jt=1,ntj
            temparray => get_tile(mod(iTileStart, 43200), mod(jTileStart, 21600), tilecache, nx, ny, nz, nh, geog_data_path, geog_sub_path, isigned, endian, scalefactor, wordsize)
            iTileEnd = iTileStart + nx - 1
            jTileEnd = jTileStart + ny - 1
            i = max(iCheckStart, iTileStart)
            tempint = i
            do while (j <= jCheckEnd .and. j <= jTileEnd)
            lat_pt = -89.99583 + (j-1) * tdim
            lat_pt = lat_pt * PI / 180.0
            do while (i <= iCheckEnd .and. i <= iTileEnd)
                u = u+1
                lon_pt = -179.99583 + (i-1) * tdim
                lon_pt = lon_pt * PI / 180.0
                closest = isClosest(lat_pt, lon_pt, latCell, lonCell, nEdgesOnCell(iCell), cellsOnCell(:,iCell), iCell)
                x = i - iTileStart + nh + 1
                y = j - jTileStart + nh + 1
                if (closest) then
                    if (ncats == 1) then
                        field(iCell) = field(iCell) + temparray(x,y,1)
                        num(1) = num(1) + 1
                    else 
                        num(int(temparray(x,y,1))) = num(int(temparray(x,y,1))) + 1
                    end if
                end if
                i = i+1
            end do
            j = j+1 
            i = tempint
            end do 
            jTileStart = jTileStart + ny
        end do JTCHECK
            iTileStart = iTileStart + nx
        end do ITCHECK

        deallocate(CoC)




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       SOLVE METHOD 3      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Same as Method 2, only copy all relevant points into rarray first then
        ! go over all of rarray and process the points. Not sure yet if faster
        ! or slower.
        case (3)
       
        allocate(rarray(iCheckStart:iCheckEnd, jCheckStart:jCheckEnd, nz))

        I3CHECK: do it=1,nti
        jTileStart = jCheckStart - mod(jCheckStart, ny) + 1
        J3CHECK: do jt=1,ntj
            temparray => get_tile(mod(iTileStart, 43200), mod(jTileStart, 21600), tilecache, nx, ny, nz, nh, geog_data_path, geog_sub_path, isigned, endian, scalefactor, wordsize)
            iTileEnd = iTileStart + nx - 1
            jTileEnd = jTileStart + ny - 1
            x = max(iCheckStart, iTileStart)
            x2 = min(iCheckEnd, iTileEnd)
            y = max(jCheckStart, jTileStart)
            y2 = min(jCheckEnd, jTileEnd)

            rarray(x:x2, y:y2, :) = temparray(x - iTileStart + 1 + nh:x2 - iTileStart + 1 + nh,y - jTileStart + 1 + nh:y2 - jTileStart + 1 + nh,:) 
            jTileStart = jTileStart+ny
        end do J3CHECK
        iTileStart = iTileStart + nx
        end do I3CHECK

        do j=jCheckStart, jCheckEnd
        lat_pt = -89.99583 + (j-1) * tdim
        lat_pt = lat_pt * PI / 180.0
        do i=iCheckStart, iCheckEnd
            lon_pt = -179.99583 + (i-1) * tdim
            lon_pt = lon_pt * PI / 180.0
            if (isClosest(lat_pt, lon_pt, latCell, lonCell, nEdgesOnCell(iCell), cellsOnCell(:,iCell), iCell)) then
                if (ncats == 1) then
                    field(iCell) = field(iCell) + rarray(i,j,1)
                    num(1) = num(1) + 1
                else 
                    num(int(rarray(i,j,1))) = num(int(rarray(i,j,1))) + 1
                end if
            end if
        end do 
        end do 

        deallocate(CoC, rarray)
        case default
        end select


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            END SOLVE METHODS               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Put the finishing touches on the value for field(iCell)
        
        if (ncats == 1) then
            if (num(1) > 0) then
                field(iCell) = field(iCell) / real(num(1))
            else
                write(0,*) "End MAIN, num = 0 for iCell = ", iCell
            end if
        else
            field(iCell) = 0
            do i=1,ncats
                if (num(i) > num(int(field(iCell)))) field(iCell) = i
            end do
        end if

        call system_clock(t3)
        write (0,*) "Time for iCell ", iCell, ": ", real(t3-t2) / real(cntrate)
    end do MAIN

    write (0,*) "Total time spent in tilecache methods = ", tilecache%tget + tilecache%tremove + tilecache%tinsert + tilecache%thit
    extime = real(t3-t1) / real(cntrate)
    write (0,*) "Total execution time for foo:", extime
    write (0,*) "Per cell avg: ", extime / nCellsSolve
    deallocate(num)
    do while (tilecache%size() > 0)
        call tilecache%lru_remove()
    end do
    write (0,*) "Exit foo"

 end subroutine foo




              
!==================================================================================================
 subroutine init_atm_static(mesh, dims, configs)
!==================================================================================================

 use mpi, only: MPI_BARRIER

!inout arguments:
 type (mpas_pool_type), intent(inout) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!local variables:
 type(proj_info):: proj

 character(len=StrKIND) :: fname
 character(len=StrKIND), pointer :: config_geog_data_path
 character(len=StrKIND), pointer :: config_landuse_data
 character(len=StrKIND+1) :: geog_data_path      ! same as config_geog_data_path, but guaranteed to have a trailing slash
 character(len=StrKIND+1) :: geog_sub_path       ! subdirectory names in config_geog_data_path, with trailing slash

 integer:: isice_lu,iswater_lu,ismax_lu

 integer:: nx,ny,nz
 integer:: endian,isigned,istatus,wordsize
 integer:: i,j,k
 integer:: iCell,iEdge,iVtx,iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(5) :: interp_list
 integer,dimension(:),allocatable  :: nhs
 integer,dimension(:,:),allocatable:: ncat
      
 real(kind=4):: scalefactor
 real(kind=4),dimension(:,:,:),allocatable:: rarray

 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND):: lat_pt,lon_pt
 real(kind=RKIND),dimension(:,:),allocatable  :: soiltemp_1deg
 real(kind=RKIND),dimension(:,:),allocatable  :: maxsnowalb
 real(kind=RKIND),dimension(:,:,:),allocatable:: vegfra

 integer, pointer :: nCells, nEdges, nVertices, maxEdges, nCellsSolve
 logical, pointer :: on_a_sphere
 real (kind=RKIND), pointer :: sphere_radius
 
 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell, verticesOnCell
 real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
 real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
 real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
 real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge
 real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
 real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
 real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
 real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
 real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
 real (kind=RKIND), dimension(:), pointer :: fEdge, fVertex

 real (kind=RKIND), dimension(:), pointer :: ter, field
 real (kind=RKIND), dimension(:), pointer :: soiltemp
 real (kind=RKIND), dimension(:), pointer :: snoalb
 real (kind=RKIND), dimension(:), pointer :: shdmin, shdmax
 real (kind=RKIND), dimension(:,:), pointer :: greenfrac
 real (kind=RKIND), dimension(:,:), pointer :: albedo12m
 integer, dimension(:), pointer :: lu_index
 integer, dimension(:), pointer :: soilcat_top
 integer, dimension(:), pointer :: soilcat_bot
 integer, dimension(:), pointer :: landmask
 character(len=StrKIND), pointer :: mminlu

!--------------------------------------------------------------------------------------------------


 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static:'

 call mpas_pool_get_config(configs, 'config_geog_data_path', config_geog_data_path)
 call mpas_pool_get_config(configs, 'config_landuse_data', config_landuse_data)

 write(geog_data_path, '(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if (geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 end if

!
! Scale all distances and areas from a unit sphere to one with radius sphere_radius
!


 call mpas_pool_get_array(mesh, 'xCell', xCell)
 call mpas_pool_get_array(mesh, 'yCell', yCell)
 call mpas_pool_get_array(mesh, 'zCell', zCell)
 call mpas_pool_get_array(mesh, 'xVertex', xVertex)
 call mpas_pool_get_array(mesh, 'yVertex', yVertex)
 call mpas_pool_get_array(mesh, 'zVertex', zVertex)
 call mpas_pool_get_array(mesh, 'xEdge', xEdge)
 call mpas_pool_get_array(mesh, 'yEdge', yEdge)
 call mpas_pool_get_array(mesh, 'zEdge', zEdge)
 call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
 call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
 call mpas_pool_get_array(mesh, 'areaCell', areaCell)
 call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
 call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
 call mpas_pool_get_array(mesh, 'latCell', latCell)
 call mpas_pool_get_array(mesh, 'lonCell', lonCell)
 call mpas_pool_get_array(mesh, 'latEdge', latEdge)
 call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
 call mpas_pool_get_array(mesh, 'latVertex', latVertex)
 call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
 call mpas_pool_get_array(mesh, 'fEdge', fEdge)
 call mpas_pool_get_array(mesh, 'fVertex', fVertex)
 
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
 call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)

 call mpas_pool_get_array(mesh, 'ter', ter)
 call mpas_pool_get_array(mesh, 'lu_index', lu_index)
 call mpas_pool_get_array(mesh, 'mminlu', mminlu)
 call mpas_pool_get_array(mesh, 'soilcat_top', soilcat_top)
 call mpas_pool_get_array(mesh, 'soilcat_bot', soilcat_bot)
 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'soiltemp', soiltemp)
 call mpas_pool_get_array(mesh, 'snoalb', snoalb)
 call mpas_pool_get_array(mesh, 'greenfrac', greenfrac)
 call mpas_pool_get_array(mesh, 'albedo12m', albedo12m)
 call mpas_pool_get_array(mesh, 'shdmin', shdmin)
 call mpas_pool_get_array(mesh, 'shdmax', shdmax)

 call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
 call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

 call mpas_pool_get_dimension(dims, 'nCells', nCells)
 call mpas_pool_get_dimension(dims, 'nEdges', nEdges)
 call mpas_pool_get_dimension(dims, 'nVertices', nVertices)
 call mpas_pool_get_dimension(dims, 'maxEdges', maxEdges)
 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve) 

 xCell = xCell * sphere_radius
 yCell = yCell * sphere_radius
 zCell = zCell * sphere_radius
 xVertex = xVertex * sphere_radius
 yVertex = yVertex * sphere_radius
 zVertex = zVertex * sphere_radius
 xEdge = xEdge * sphere_radius
 yEdge = yEdge * sphere_radius
 zEdge = zEdge * sphere_radius
 dvEdge = dvEdge * sphere_radius
 dcEdge = dcEdge * sphere_radius
 areaCell = areaCell * sphere_radius**2.0
 areaTriangle = areaTriangle * sphere_radius**2.0
 kiteAreasOnVertex = kiteAreasOnVertex * sphere_radius**2.0

 where (lonVertex < 0) lonVertex = lonVertex+2*PI
 where (lonCell < 0) lonCell = lonCell + 2*PI

!
! Initialize Coriolis parameter field on edges and vertices
!
 do iEdge=1,nEdges
    fEdge(iEdge)  = 2.0 * omega * sin(latEdge(iEdge))
 end do
 do iVtx=1,nVertices
    fVertex(iVtx) = 2.0 * omega * sin(latVertex(iVtx))
 end do

!
! Compute weights used in advection and deformation calculation
!
 call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius) 
 call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius) 


!
! Set land use and soil category parameters for water and ice
!
 surface_input_select0: select case(trim(config_landuse_data))
    case('USGS')
       isice_lu = 24
       iswater_lu = 16
       ismax_lu = 24
       write(mminlu,'(a)') 'USGS'
    case('MODIFIED_IGBP_MODIS_NOAH')
       isice_lu = 15
       iswater_lu = 17
       ismax_lu = 20
       write(mminlu,'(a)') 'MODIFIED_IGBP_MODIS_NOAH'
    case default
         write(0,*) '*****************************************************************'
         write(0,*) 'Invalid land use dataset '''//trim(config_landuse_data)//''' selected for config_landuse_data'
         write(0,*) '   Possible options are: ''USGS'', ''MODIFIED_IGBP_MODIS_NOAH'''
         write(0,*) '*****************************************************************'
         call mpas_dmpar_global_abort('Please correct the namelist.')
 end select surface_input_select0


!
! Interpolate HGT
!
!nx = 126
!ny = 126
 nx = 1200
 ny = 1200
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 2
 scalefactor = 1.0
 geog_sub_path = 'topo_30s/'

 call foo2(mesh, dims, configs, geog_data_path, geog_sub_path, nx, ny, nz, 3, &
                 field, 1, isigned, endian, scalefactor, wordsize)

! call foo(geog_data_path, geog_sub_path, nx, ny, nz, 3, nCells, nCellsSolve, nVertices, latCell, lonCell, latVertex, lonVertex, &
!                verticesOnCell, cellsOnCell, nEdgesOnCell, field, 1, isigned, endian, scalefactor, wordsize)
 ter(1:nCells) = field(:)
 if (associated(field)) deallocate(field)
 call MPI_BARRIER(MPI_COMM_WORLD, istatus)
 write(0,*) '--- end interpolate TER'
! isigned  = 1
! endian   = 0
! wordsize = 2
! scalefactor = 1.0
! allocate(rarray(nx,ny,nz))
! allocate(nhs(nCells))
! nhs(:) = 0
! ter(:) = 0.0
!
! do jTileStart = 1,20401,ny-6
!    jTileEnd = jTileStart + ny - 1 - 6
!
!    do iTileStart=1,42001,nx-6
!       iTileEnd = iTileStart + nx - 1 - 6
!       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
!             'topo_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
!       write(0,*) trim(fname)
!
!       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
!                         scalefactor,wordsize,istatus)
!       call init_atm_check_read_error(istatus, fname)
!
!       iPoint = 1
!       do j=4,ny-3
!       do i=4,nx-3
!          lat_pt = -89.99583  + (jTileStart + j - 5) * 0.0083333333
!          lon_pt = -179.99583 + (iTileStart + i - 5) * 0.0083333333
!          lat_pt = lat_pt * PI / 180.0
!          lon_pt = lon_pt * PI / 180.0
!
!          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
!                                nEdgesOnCell,cellsOnCell, &
!                                latCell,lonCell)
!          ter(iPoint) = ter(iPoint) + rarray(i,j,1)
!          nhs(iPoint) = nhs(iPoint) + 1
!       end do
!       end do
!
!    end do
! end do
!
! do iCell = 1,nCells
!    ter(iCell) = ter(iCell) / real(nhs(iCell))
! end do
! deallocate(rarray)
! deallocate(nhs)
! write(0,*) '--- end interpolate TER'


!
! Interpolate LU_INDEX
!
 surface_input_select1: select case(trim(config_landuse_data))
    case('USGS')
       geog_sub_path = 'landuse_30s/'
    case('MODIFIED_IGBP_MODIS_NOAH')
       geog_sub_path = 'modis_landuse_20class_30s/'
    case default
       write(0,*) '*****************************************************************'
       write(0,*) 'Invalid land use dataset '''//trim(config_landuse_data)//''' selected for config_landuse_data'
       write(0,*) '   Possible options are: ''USGS'', ''MODIFIED_IGBP_MODIS_NOAH'''
       write(0,*) '*****************************************************************'
       call mpas_dmpar_global_abort('Please correct the namelist.')
 end select surface_input_select1
 nx = 1200
 ny = 1200
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 1
 scalefactor = 1.0

 call foo2(mesh, dims, configs, geog_data_path, geog_sub_path, nx, ny, nz, 0, &
                 field, ismax_lu, isigned, endian, scalefactor, wordsize)
! call foo(geog_data_path, geog_sub_path, nx, ny, nz, 0, nCells, nCellsSolve, nVertices, latCell, lonCell, latVertex, lonVertex, &
!                verticesOnCell, cellsOnCell, nEdgesOnCell, field, ismax_lu, isigned, endian, scalefactor, wordsize)
 lu_index(1:nCells) = int(field(:))
 if (associated(field)) deallocate(field)

 where(lu_index == 0) lu_index = 1
 call MPI_BARRIER(MPI_COMM_WORLD, istatus)
 !allocate(rarray(nx,ny,nz))
 !allocate(ncat(ismax_lu,nCells))
 !ncat(:,:) = 0
 !lu_index(:) = 0.0
!
! do jTileStart = 1,20401,ny
!    jTileEnd = jTileStart + ny - 1
!
!    do iTileStart = 1,42001,nx
!       iTileEnd = iTileStart + nx - 1
!       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
!             trim(geog_sub_path),iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
!       write(0,*) trim(fname)
!
!       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
!                         scalefactor,wordsize,istatus)
!       call init_atm_check_read_error(istatus, fname)
!
!       iPoint = 1
!       do j=1,ny
!       do i=1,nx
!!
!! The MODIS dataset appears to have zeros at the South Pole, possibly other places, too
!!
!if (rarray(i,j,1) == 0) cycle
!
!          lat_pt = -89.99583  + (jTileStart + j - 2) * 0.0083333333
!          lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
!          lat_pt = lat_pt * PI / 180.0
!          lon_pt = lon_pt * PI / 180.0
!
!          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
!                                nEdgesOnCell,cellsOnCell, &
!                                latCell,lonCell)
!          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
!       end do
!       end do
!
!    end do
! end do
!
! do iCell = 1,nCells
!    lu_index(iCell) = 1
!    do i = 2,ismax_lu
!       if(ncat(i,iCell) > ncat(lu_index(iCell),iCell)) then
!          lu_index(iCell) = i
!       end if
!    end do
! end do
! deallocate(rarray)
! deallocate(ncat)
 write(0,*) '--- end interpolate LU_INDEX'


!
! Interpolate SOILCAT_TOP
!
 nx = 1200
 ny = 1200
 nz = 1
 isigned     = 1
 endian      = 0
 wordsize    = 1
 scalefactor = 1.0
 soilcat_top(:) = 0.0

! geog_sub_path = 'soiltype_top_30s/'
! call foo(geog_data_path, geog_sub_path, nx, ny, nz, 0, nCells, nCellsSolve, nVertices, latCell, lonCell, latVertex, lonVertex, &
!                verticesOnCell, cellsOnCell, nEdgesOnCell, field, ismax_lu, isigned, endian, scalefactor, wordsize)
! soilcat_top(1:nCells) = int(field(:))
! if (associated(field)) deallocate(field)

 write(0,*) '--- end interpolate SOILCAT_TOP'


!
! Interpolate SOILCAT_BOT
!
 nx = 1200
 ny = 1200
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 1
 scalefactor = 1.0
 soilcat_bot(:) = 0.0

! geog_sub_path = 'soiltype_bot_30s/'
! call foo(geog_data_path, geog_sub_path, nx, ny, nz, 0, nCells, nCellsSolve, nVertices, latCell, lonCell, latVertex, lonVertex, &
!                verticesOnCell, cellsOnCell, nEdgesOnCell, field, ismax_lu, isigned, endian, scalefactor, wordsize)
! soilcat_bot(1:nCells) = int(field(:))
! if (associated(field)) deallocate(field)
!
 write(0,*) '--- end interpolate SOILCAT_BOT'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KLUDGE TO FIX SOIL TYPE OVER ANTARCTICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 where (lu_index == isice_lu) soilcat_top = 16
 where (lu_index == isice_lu) soilcat_bot = 16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CORRECT INCONSISTENT SOIL AND LAND USE DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iCell = 1,nCells
    if (lu_index(iCell) == iswater_lu .or. &
        soilcat_top(iCell) == 14 .or. &
        soilcat_bot(iCell) == 14) then
        if (lu_index(iCell) /= iswater_lu) then
            write(0,*) 'Turning lu_index into water at ', iCell
            lu_index(iCell) = iswater_lu
        end if
        if (soilcat_top(iCell) /= 14) then
            write(0,*) 'Turning soilcat_top into water at ', iCell
            soilcat_top(iCell) = 14
        end if
        if (soilcat_bot(iCell) /= 14) then
            write(0,*) 'Turning soilcat_bot into water at ', iCell
            soilcat_bot(iCell) = 14
        end if
    end if
 end do


!
! Derive LANDMASK
!
landmask(:) = 0
 do iCell=1, nCells
    if (lu_index(iCell) /= iswater_lu) landmask(iCell) = 1
 end do
 write(0,*) '--- end interpolate LANDMASK'


!
! Interpolate SOILTEMP:
!
 nx = 186
 ny = 186
 nz = 1
 isigned  = 0
 endian   = 0
 wordsize = 2
 scalefactor = 0.01
 allocate(rarray(nx,ny,nz))
 allocate(soiltemp_1deg(-2:363,-2:183))
soiltemp(:) = 0.0

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'soiltemp_1deg/',1,'-',180,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned, endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus, fname)
 soiltemp_1deg(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
            'soiltemp_1deg/',181,'-',360,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname, len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                        scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 soiltemp_1deg(181:363,-2:183) = rarray(4:186,1:186,1)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
  
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if (x >= 360.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 179.0) y = 179.0
       soiltemp(iCell) = interp_sequence(x,y,1,soiltemp_1deg,-2,363,-2,183, &
                                           1,1,0.0_RKIND,interp_list,1)
    else
       soiltemp(iCell) = 0.0
    end if

 end do
 deallocate(rarray)
 deallocate(soiltemp_1deg)
 write(0,*) '--- end interpolate SOILTEMP'

!soiltemp(:) = 1.0
!
! Interpolate SNOALB
!
 nx = 186
 ny = 186
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(maxsnowalb(-2:363,-2:183))
 snoalb(:) = 0.0

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'maxsnowalb/',1,'-',180,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 maxsnowalb(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'maxsnowalb/',181,'-',360,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus, fname)
 maxsnowalb(181:363,-2:183) = rarray(4:186,1:186,1)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
  
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if (x >= 360.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 179.0) y = 179.0
       snoalb(iCell) = interp_sequence(x,y,1,maxsnowalb,-2,363,-2,183, &
                                         1,1,0.0_RKIND,interp_list,1)
    else
       snoalb(iCell) = 0.0
    end if

 end do
 snoalb(:) = snoalb(:) / 100.0
 deallocate(rarray)
 deallocate(maxsnowalb)
 write(0,*) '--- end interpolate SNOALB'
snoalb(:) = 1.0

!
! Interpolate GREENFRAC
!
 nx = 1256
 ny = 1256
 nz = 12
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(vegfra(-2:2503,-2:1253,12))
 greenfrac(:,:) = 0.0

 call map_set(PROJ_LATLON, proj,    &
              latinc = 0.144_RKIND, &
              loninc = 0.144_RKIND, &
              knowni = 1.0_RKIND,   &
              knownj = 1.0_RKIND,   &
              lat1 = -89.928_RKIND, &
              lon1 = -179.928_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'greenfrac/',1,'-',1250,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'greenfrac/',1251,'-',2500,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 vegfra(1251:2503,-2:1253,1:12) = rarray(4:1256,1:1256,1:12)

 do iCell = 1,nCells

    if (landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if(x >= 2500.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 1249.0) y = 1249.0
       do k = 1,12
          greenfrac(k,iCell) = interp_sequence(x,y,k,vegfra,-2,2503,-2,1253, &
                                                 1,12,-1.e30_RKIND,interp_list,1)
       end do
    else
       greenfrac(:,iCell) = 0.0
    end if
    shdmin(iCell) = minval(greenfrac(:,iCell))
    shdmax(iCell) = maxval(greenfrac(:,iCell))
      
 end do
 deallocate(rarray)
 deallocate(vegfra)
 write(0,*) '--- end interpolate GREENFRAC'


!
! Interpolate ALBEDO12M
!
 nx = 1256
 ny = 1256
 nz = 12
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(vegfra(-2:2503,-2:1253,12))
 albedo12m(:,:) = 0.0

 call map_set(PROJ_LATLON, proj,    &
              latinc = 0.144_RKIND, &
              loninc = 0.144_RKIND, &
              knowni = 1.0_RKIND,   &
              knownj = 1.0_RKIND,   &
              lat1 = -89.928_RKIND, &
              lon1 = -179.928_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'albedo_ncep/',1,'-',1250,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor, wordsize, istatus)
 call init_atm_check_read_error(istatus,fname)
 vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'albedo_ncep/',1251,'-',2500,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 vegfra(1251:2503,-2:1253,1:12) = rarray(4:1256,1:1256,1:12)

 do iCell = 1,nCells

    if (landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if(x >= 2500.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 1249.0) y = 1249.0
       do k = 1,12
          albedo12m(k,iCell) = interp_sequence(x,y,k,vegfra,-2,2503,-2,1253, &
                                                 1,12,0.0_RKIND,interp_list,1)
       end do
    else
       albedo12m(:,iCell) = 8.0
    end if
 end do
 deallocate(rarray)
 deallocate(vegfra)
 write(0,*) '--- end interpolate ALBEDO12M'

 call MPI_BARRIER(MPI_COMM_WORLD, istatus)
 write (0,*) "Barrier exited with status ", istatus

 end subroutine init_atm_static

!==================================================================================================
 subroutine init_atm_static_orogwd(mesh, dims, configs)
!==================================================================================================

!inout arguments:
 type (mpas_pool_type), intent(inout) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!local variables:
 type(proj_info):: proj

 character(len=StrKIND) :: mess
 character(len=StrKIND) :: fname
 character(len=StrKIND) :: dir_gwdo
 character(len=StrKIND), pointer :: config_geog_data_path
 character(len=StrKIND+1) :: geog_data_path, geog_sub_path      ! same as config_geog_data_path, but guaranteed to have a trailing slash

 integer, pointer :: nCells, nEdges, maxEdges

 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell
 integer, dimension(:), pointer :: landmask

 integer:: nx,ny,nz
 integer:: endian,isigned,istatus,wordsize
 integer:: i,j
 integer:: iCell,iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(5) :: interp_list
 integer,dimension(:),allocatable:: nhs

 real(kind=4):: scalefactor
 real(kind=4),dimension(:,:,:),allocatable:: rarray

 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND):: lat_pt,lon_pt
 real(kind=RKIND):: dx,dy,known_lat,known_lon,known_x,known_y
 real(kind=RKIND):: minMeshD,maxMeshD
 real(kind=RKIND):: mindcEdge,maxdcEdge
 real(kind=RKIND),dimension(:,:),allocatable:: xarray

 real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
 real(kind=RKIND), dimension(:), pointer :: meshDensity
 real(kind=RKIND), dimension(:), pointer :: dcEdge
 real(kind=RKIND), dimension(:), pointer :: varsso
 real(kind=RKIND), dimension(:), pointer :: con, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, var2d


 call mpas_pool_get_dimension(dims, 'nCells', nCells)
 call mpas_pool_get_dimension(dims, 'nEdges', nEdges)
 call mpas_pool_get_dimension(dims, 'maxEdges', maxEdges)

 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'latCell', latCell)
 call mpas_pool_get_array(mesh, 'lonCell', lonCell)
 call mpas_pool_get_array(mesh, 'varsso', varsso)
 call mpas_pool_get_array(mesh, 'meshDensity', meshDensity)
 call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
 call mpas_pool_get_array(mesh, 'con', con)
 call mpas_pool_get_array(mesh, 'oa1', oa1)
 call mpas_pool_get_array(mesh, 'oa2', oa2)
 call mpas_pool_get_array(mesh, 'oa3', oa3)
 call mpas_pool_get_array(mesh, 'oa4', oa4)
 call mpas_pool_get_array(mesh, 'ol1', ol1)
 call mpas_pool_get_array(mesh, 'ol2', ol2)
 call mpas_pool_get_array(mesh, 'ol3', ol3)
 call mpas_pool_get_array(mesh, 'ol4', ol4)
 call mpas_pool_get_array(mesh, 'var2d', var2d)


 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static_orogwd:'

 call mpas_pool_get_config(configs, 'config_geog_data_path', config_geog_data_path)

 write(geog_data_path, '(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if (geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 end if

!
! Interpolate VARSSO:
 geog_sub_path = 'varsso'
 varsso(:) = 0.0_RKIND
 nx = 600
 ny = 600
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 4
 scalefactor = 1.0

! dx = 0.00833333
! dy = 0.00833333
! known_x = 1.0
! known_y = 1.0
! known_lat = -59.99583
! known_lon = -179.99583
!  
! allocate(rarray(nx,ny,nz))
! allocate(nhs(nCells))
! nhs(:) = 0
! rarray(:,:,:) = 0._RKIND
! do jTileStart = 1,13801,ny
!    jTileEnd = jTileStart + ny - 1
!
!    do iTileStart = 1,42601,nx
!       iTileEnd = iTileStart + nx -1
!       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)//'varsso/', &
!             iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
!       write(0,*) trim(fname)
!
!       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
!                         scalefactor,wordsize,istatus)
!       call init_atm_check_read_error(istatus,fname)
!
!       iPoint = 1
!       do j = 1,ny
!       do i = 1,nx
!          lat_pt = known_lat + (jTileStart + j - 2) * dy
!          lon_pt = known_lon + (iTileStart + i - 2) * dx
!          lat_pt = lat_pt * PI / 180.0
!          lon_pt = lon_pt * PI / 180.0
!
!          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
!                                nEdgesOnCell,cellsOnCell, &
!                                latCell,lonCell)
!          varsso(iPoint) = varsso(iPoint) + rarray(i,j,1)
!          nhs(iPoint) = nhs(iPoint) + 1
!       enddo
!       enddo
!
!    enddo
! enddo
!
! do iCell = 1,nCells
!    if(nhs(iCell) .gt. 0) &
!       varsso(iCell) = varsso(iCell) / real(nhs(iCell))
! enddo
! deallocate(rarray)
! deallocate(nhs)
 write(0,*) '--- end interpolate VARSSO'

!... statistic fields needed for the parameterization of gravity wavwe drag over orography. The
!input directory depends on the mesh resolution, and the mesh must be a uniform mesh.
 minMeshD  = minval(meshDensity(1:nCells))
 maxMeshD  = maxval(meshDensity(1:nCells))
 mindcEdge = minval(dcEdge(1:nEdges))
 maxdcEdge = maxval(dcEdge(1:nEdges))

 write(0,*)
 write(0,*) 'BEGIN INTERPOLATION OF STATISTICAL FIELDS FOR GRAVITY WAVE DRAG OVER OROGRAPHY'
 write(0,*) 'min MeshD  =', minMeshD
 write(0,*) 'max MeshD  =', maxMeshD
 write(0,*) 'min dcEdge =', mindcEdge
 write(0,*) 'max dcEdge =', maxdcEdge

 dir_gwdo = '   '
 if(minMeshD == 1.0_RKIND .and. maxMeshD == 1.0_RKIND) then
    !... uniform 10242 mesh:
    if(mindcEdge .ge. 200000._RKIND .and. maxdcEdge .lt. 260000._RKIND) then
       dir_gwdo = 'orogwd_2deg'
    elseif(mindcEdge .ge. 90000._RKIND .and. maxdcEdge .lt. 150000_RKIND) then
       dir_gwdo = 'orogwd_1deg'
    elseif(mindcEdge .ge. 40000._RKIND .and. maxdcEdge .lt. 70000._RKIND) then
       dir_gwdo = 'orogwd_30m'
    else
       write(0,*)
!      write(mess,*) 'GWDO: Interpolation not available. The initialization will abort'
!      call physics_error_fatal(mess)
       write(mess,*) 'GWDO: Interpolation not available. Set config_gwdo_scheme = .false.'
       return
    endif
 else
    write(0,*)
!   write(mess,*) 'GWDO: The input mesh must be a uniform mesh. The initialization will abort'
!   call physics_error_fatal(mess)
    write(mess,*) 'GWDO: The input mesh must be a uniform mesh. Set config_gwdo_scheme = .false.'
    return
 endif
 write(0,*) 'dir_gwdo   =    ', trim(dir_gwdo)
 write(0,*)

!
! Interpolate CON:
!
 con(:) = 0.0_RKIND

 con_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select con_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/con/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       con(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate CON'

!
! Interpolate OA1:
!
 oa1(:) = 0.0_RKIND

 oa1_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa1_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa1/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa1(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA1'

!
! Interpolate OA2:
 oa2(:) = 0.0_RKIND

 oa2_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa2_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa2/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa2(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
     endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA2'

!
! Interpolate OA3:
!
 oa3(:) = 0.0_RKIND

 oa3_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa3_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa3(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA3'

!
! Interpolate OA4:
!
 oa4(:) = 0.0_RKIND

 oa4_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa4_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa4/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa4(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA4'

!
! Interpolate OL1:
!
 ol1(:) = 0.0_RKIND

 ol1_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol1_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol1/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol1(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL1'

!
! Interpolate OL2:
!
 ol2(:) = 0.0_RKIND

 ol2_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol2_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol2/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol2(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL2'

!
! Interpolate OL3:
!
 ol3(:) = 0.0_RKIND

 ol3_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol3_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol3(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL3'

!
! Interpolate OL4:
!
 ol4(:) = 0.0_RKIND

 ol4_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol4_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol4/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol4(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL4'

!
! Interpolate VAR2D:
!
 var2d(:) = 0.0_RKIND

 var2d_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny = 90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 4
       scalefactor = 0.02
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 4
       scalefactor = 0.02
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 4
       scalefactor = 0.02
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.02
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select var2d_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/var/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)


 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       var2d(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                      0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate VAR2D'

 end subroutine init_atm_static_orogwd

!==================================================================================================
 subroutine init_atm_check_read_error(istatus, fname)
!==================================================================================================
 implicit none

 integer, intent(in) :: istatus
 character (len=*), intent(in) :: fname

 if (istatus /= 0) then
     call mpas_dmpar_global_abort('ERROR: Could not read file '//trim(fname))
 end if

 end subroutine init_atm_check_read_error

!==================================================================================================
 integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                               nEdgesOnCell, cellsOnCell, latCell, lonCell)
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: target_lat, target_lon
 integer, intent(in) :: start_cell
 integer, intent(in) :: nCells, maxEdges
 integer, dimension(nCells), intent(in) :: nEdgesOnCell
 integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
 real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell

 integer :: i
 integer :: iCell
 integer :: current_cell
 real (kind=RKIND) :: current_distance, d
 real (kind=RKIND) :: nearest_distance

 nearest_cell = start_cell
 current_cell = -1

 do while (nearest_cell /= current_cell)
    current_cell = nearest_cell
    current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                       target_lon, 1.0_RKIND)
    nearest_cell = current_cell
    nearest_distance = current_distance
    do i = 1, nEdgesOnCell(current_cell)
       iCell = cellsOnCell(i,current_cell)
       if (iCell <= nCells) then
          d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_RKIND)
          if (d < nearest_distance) then
             nearest_cell = iCell
             nearest_distance = d
          end if
       end if
    end do
 end do

 end function nearest_cell

!==================================================================================================
 real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

!Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
!sphere with given radius.
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
 real (kind=RKIND) :: arg1

 arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
              cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
 sphere_distance = 2.*radius*asin(arg1)

 end function sphere_distance


 subroutine arrDist2(lat, lon, radius, x, y, z, dist2)

    real(kind=RKIND), dimension(:,:,:), intent(in) :: lat, lon
    real(kind=RKIND), dimension(:,:,:), intent(inout) :: dist2
    real(kind=RKIND), intent(in) :: radius, x, y, z

    dist2 = (radius * sin(lat)-z)**2 
    dist2 = dist2 + (radius * cos(lon) * cos(lat)-x)**2
    dist2 = dist2 + (radius * sin(lon) * cos(lat)-y)**2
 end subroutine arrDist2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CONVERT_LX
   !
   ! Convert (lat,lon) to an (x, y, z) location on a sphere with specified
   ! radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine con_lx(lat, lon, radius, x, y, z)

      implicit none

      real (kind=RKIND), intent(in) :: radius, lat, lon
      real (kind=RKIND), intent(out) :: x, y, z

      z = radius * sin(lat)
      x = radius * cos(lon) * cos(lat)
      y = radius * sin(lon) * cos(lat)
   end subroutine con_lx

    


!==================================================================================================
 end module mpas_init_atm_static
!==================================================================================================
