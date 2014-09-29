Program LocalData
!
! Purpose:
!   To generate the NewList and CellEnd array by linked-list cells
!
! Record of revisions:
! 	Date		Programmer		Description of change
!       ====            ==========              =====================
!    Nov 3, 2010        Zhou Lvwen                  Original Code
!

implicit none
integer, parameter :: Nc = 4		! Number of cells
! Cell Index of all atoms or particles
integer, dimension(8) :: CellInd = [4, 1, 4, 1, 3, 3, 4, 1]
integer, dimension(Nc+1) :: CellEnd = 0 ! CellEnd array
integer, dimension(Nc+1) :: Temp = 0	! A copy of CellEnd for generate NewList
integer, dimension(8) :: NewList = 0    ! NewList array
integer :: i				! Loop index
integer :: Icell			! Cell index of particle i

do i = 1, 8
   Icell = CellInd(i)
   CellEnd(Icell+1:) = CellEnd(Icell+1:) + 1
end do

Temp = CellEnd

do i = 8, 1, -1
   Icell = CellInd(i)
   NewList( CellEnd(Icell+1) ) = i
   Temp(Icell+1) = Temp(Icell+1) - 1
end do

write(*,'(A, 8(1X,I2), 1X, A)') ' NewList = [', NewList, ']'
write(*,'(A, 5(1X,I2), 1X, A)') ' CellEnd = [', CellEnd, ']'

End program LocalData


