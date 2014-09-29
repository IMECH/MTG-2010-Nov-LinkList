program LinkedList
!
! Purpose:
!       To give An example about half neighbors search base on Linked-list cells
!
! Record of revisions:
!       Date            Programmer              Description of change
!       ====            ==========              =====================
!    Nov 1, 2010        Zhou Lvwen                  Original Code
!


implicit none
integer, parameter :: natom = 100                       ! Number of atoms or particles
integer, dimension(3), parameter :: Nc = [3, 3, 3]      ! Number of cells in x, y and z dimension
real   , dimension(3), parameter :: boxl = [6, 6, 6]    ! Box longth in x, y, and z dimension
real   , dimension(natom,3) :: r                        ! Position of atoms
integer, dimension(natom,3) :: CellSub                  ! Multiple subscripts of cells: atoms 
integer, dimension(3) :: CellSub_i, CellSub_j           ! Multiple subscripts of cells: atom i,j
integer, dimension(Nc(1)*Nc(2)*Nc(3)) :: head  = 0      ! 
integer, dimension(natom) ::  link = 0                  !
integer :: i, j, k, l                                   ! Loop index
integer :: CellInd_i, CellInd_j                         ! Index of cell i,j
integer :: atom_i, atom_j                               ! Index of atom i,j
integer, dimension(14,3) :: Dc                          ! 14 Neighbors of cell(0,0,0) 
integer :: sub2ind                                      ! Function: Linear index from multiple subscripts


! Build neighbor cells list
l = 1
do i = -1, 1
   do j = -1, 1
     do k = -1, 1
        if ((i + j + k)>=0 .and. k >= 0) then
           Dc(l,:) = [i, j, k]
           l = l+1
        end if
     end do
   end do
end do

! Generate the random position of all atoms or particles 
do i = 1, 3
   do j = 1, natom
      call random_number(r(j,i))
   end do
end do

! Calculate cell subscript values of all atoms or particles
do i = 1, 3
   r(:, i) = r(:, i) * boxl(i)
   CellSub(:,i) = ceiling( r(:,i)/boxl(i) * Nc(i))
end do

! Generate the link array and head array by linked-list cells
do i = 1, natom
   CellInd_i = sub2ind(Nc, CellSub(i,:))
   write(*,*) CellSub(i,:), CellInd_i
   ! Link to the previous occupant (or EMPTY if you're the 1st)
   link(i) = head(CellInd_i)

   ! The last one goes to the header
   head(CellInd_i) = i
end do



do i = 1, Nc(1)*Nc(2)*Nc(3)
   CellInd_i = i
   atom_i = head(CellInd_i)
   do while (atom_i /= 0)

      do j = 1, 14
         CellSub_j =  CellSub(atom_i,:) + Dc(j,:)
        
         ! Periodic boundary condition
         do k = 1, 3
            if (CellSub_j(k) <  1    ) CellSub_j(k) = CellSub_j(k) + Nc(k)
            if (CellSub_j(k) >  Nc(k)) CellSub_j(k) = CellSub_j(k) - Nc(k)
         end do

         CellInd_j = sub2ind(Nc, CellSub_j)
         
         if (CellInd_i == CellInd_j) then
            atom_j = link(atom_i)
         else
            atom_j = head(CellInd_j)
         end if

         do while(atom_j /= 0)
            !*************************************************
            !compute the force between atom_i and atom_j in here
            ! 
            ! r_ij = r(atom_i,:) - r(atom_j,:)
            !
            !*************************************************
           atom_j = link(atom_j)
         end do
      end do
      atom_i = link(atom_i)
   end do
end do


end program LinkedList


integer function sub2ind(siz,sub)
!
! Purpose:
! 	To determine  the equivalent single index corresponding  to a 
!       given set of subscript values.
!       
!       siz =[siz_x, siz_y, siz_z]   sub = [sub_x, sub_y, sub_z]
!                                       
!                                         /  sub_x  \
!       ind = [siz_x, siz_y-1, siz_z-1] * |  sub_y  |
!                                         \  sub_z  /    
! Record of revisions:
! 	Date		Programmer		Description of change
!       ====            ==========              =====================
!    Nov 3, 2010        Zhou Lvwen                  Original Code
!

implicit none

integer, intent(in), dimension(3) :: siz	! arrary size
integer, intent(in), dimension(3) :: sub	! subscript value


sub2ind = dot_product(sub-[0,1,1],[1, siz(1), siz(2)*siz(3)])

end function sub2ind
