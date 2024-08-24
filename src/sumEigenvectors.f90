! This file automatically generated from sumEigenvectors.bf90 with bpp.
!
! =======================================================================
! ============ Optimized routine for eveSolver ====================
!              SUM EIGENVECTOR EXPANSION
! 
! =======================================================================
!




! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : beginLoops
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------




! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================




! Argument list


subroutine sumEigenvectors( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,u,w,amp,numEigs,ierr )
!======================================================================
!             u = SUM_{ieg} amp(ieg)*w(.,ieg)
!======================================================================
    implicit none
    integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

    real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:*)
    real w(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
    real amp(0:*)

  ! integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
    integer numEigs,ierr

  !     ---- local variables -----
    integer i1,i2,i3,ieg
    real ampi

    write(*,'(" sumEigenvectors: numEigs=",i6," n1a,n1b,n2a,n2b,n3a,n3b=",6i4)') numEigs,n1a,n1b,n2a,n2b,n3a,n3b

    if( .true. )then
    ! --- This way is faster than the version below ---
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
            u(i1,i2,i3,0)=0.
          end do
          end do
          end do

        do ieg=0,numEigs-1
            ampi = amp(ieg)
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                u(i1,i2,i3,0) = u(i1,i2,i3,0) + ampi*w(i1,i2,i3,ieg)
              end do
              end do
              end do
        end do

    else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b

            u(i1,i2,i3,0)=0.
            do ieg=0,numEigs-1
                u(i1,i2,i3,0) = u(i1,i2,i3,0) + amp(ieg)*w(i1,i2,i3,ieg)
            end do

          end do
          end do
          end do
    end if 
        
    return
    end