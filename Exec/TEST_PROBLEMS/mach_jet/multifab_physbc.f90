module multifab_physbc_module

  use multifab_module
  use define_bc_module
  use bl_types
  use bl_error_module
  use bc_module
  use bl_constants_module
  use inlet_bc_module

  implicit none

  private

  public :: multifab_physbc

contains

    subroutine multifab_physbc(s,start_scomp,start_bccomp,ncomp,the_bc_level, &
                               time_in,dx_in,prob_lo_in,prob_hi_in)

      integer        , intent(in   )           :: start_scomp,start_bccomp,ncomp
      type(multifab) , intent(inout)           :: s
      type(bc_level) , intent(in   )           :: the_bc_level
      real(kind=dp_t), intent(in   ), optional :: time_in,dx_in(:),prob_lo_in(:),prob_hi_in(:)

      ! Local
      integer                  :: lo(get_dim(s)),hi(get_dim(s))
      integer                  :: i,ng,dm,scomp,bccomp
      real(kind=dp_t)          :: time,dx(get_dim(s))
      real(kind=dp_t)          :: prob_lo(get_dim(s)),prob_hi(get_dim(s))
      real(kind=dp_t), pointer :: sp(:,:,:,:)

      write(*,*) "HELLO????"

      if (nghost(s) == 0) return

      ! set optional arguments
      time    = ZERO
      dx      = ZERO
      prob_lo = ZERO
      prob_hi = ZERO
      if (present(time_in))       time = time_in
      if (present(dx_in))           dx = dx_in
      if (present(prob_lo_in)) prob_lo = prob_lo_in
      if (present(prob_hi_in)) prob_hi = prob_hi_in

      ng = nghost(s)
      dm = get_dim(s)

      !$OMP PARALLEL DO PRIVATE(i,sp,lo,hi,scomp,bccomp)
      do i=1,nfabs(s)
         sp => dataptr(s,i)
         lo = lwb(get_box(s,i))
         hi = upb(get_box(s,i))
         select case (dm)
         case (1)
            do scomp=start_scomp,start_scomp+ncomp-1
               bccomp = start_bccomp + scomp - start_scomp
               call physbc_1d(sp(:,1,1,scomp), lo, hi, ng, &
                              the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                              time, dx, prob_lo, prob_hi, bccomp)
            end do
         case (2)
            do scomp=start_scomp,start_scomp+ncomp-1
               bccomp = start_bccomp + scomp - start_scomp
               call physbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                              the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                              time, dx, prob_lo, prob_hi, bccomp)
            end do
         case (3)
            do scomp=start_scomp,start_scomp+ncomp-1
               bccomp = start_bccomp + scomp - start_scomp
               call physbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                              the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                              time, dx, prob_lo, prob_hi, bccomp)
            end do
         end select
      end do
      !$OMP END PARALLEL DO

    end subroutine multifab_physbc

  subroutine physbc_1d(s,lo,hi,ng,bc,time,dx,prob_lo,prob_hi,icomp)

      integer        , intent(in   ) :: lo(:),hi(:),ng
      real(kind=dp_t), intent(inout) :: s(lo(1)-ng:)
      integer        , intent(in   ) :: bc(:,:)
      real(kind=dp_t), intent(in   ) :: time,dx(:),prob_lo(:),prob_hi(:)
      integer        , intent(in   ) :: icomp

    call bl_error("physbc_1d not written")

  end subroutine physbc_1d

  subroutine physbc_2d(s,lo,hi,ng,bc,time,dx,prob_lo,prob_hi,icomp)

    use probdata_module, only: inlet_mach
    use base_state_geometry_module, only: nr_fine, dr

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: time,dx(:),prob_lo(:),prob_hi(:)
    integer        , intent(in   ) :: icomp

    !     Local variables
    integer :: i

    real(kind=dp_t) :: A,B,x

    write(*,*) "I'm here!!"

    A = 4.5d-2
    B = 1.d2

    if (ng == 0) return

    !--------------------------------------------------------------------------
    ! lower X
    !--------------------------------------------------------------------------
    if (bc(1,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,1) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call bl_error("physbc_2d: bc(1,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper X
    !--------------------------------------------------------------------------
    if (bc(1,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else if (bc(1,2) .eq. PERIODIC) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call bl_error("physbc_2d: bc(1,2) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! lower Y
    !--------------------------------------------------------------------------
    if (bc(2,1) .eq. EXT_DIR) then
       ! xvel
       if (icomp .eq. 1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = 0.d0
       ! yvel
       if (icomp .eq. 2) then
          do i=lo(1)-ng,hi(1)+ng
             x = (dble(i)+0.5d0)*dr(nr_fine)
             ! inflow is Mach number 0.01 front with a Mach number 0.1 bump in the middle
             s(i,lo(2)-ng:lo(2)-1) = (inlet_mach/1.d-1)* &
                  INLET_CS*(1.d-2 + A*(tanh(B*(x-0.40d0)) + tanh(B*(0.6d0-x))))
          end do
       end if
       ! rho
       if (icomp .eq. 3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHO
       ! rhoh
       if (icomp .eq. 4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHOH
       ! species
       if (icomp .eq. 5) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_RHO
       ! temperature
       if (icomp .eq. 6) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_TEMP
       ! tracer
       if (icomp .eq. 7) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = 0.d0
    else if (bc(2,1) .eq. FOEXTRAP) then
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call bl_error("physbc_2d: bc(2,1) not yet supported")
    end if

    !--------------------------------------------------------------------------
    ! upper Y
    !--------------------------------------------------------------------------
    if (bc(2,2) .eq. FOEXTRAP) then
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = &
               ( 15.d0 * s(i,hi(2)  ) &
                -10.d0 * s(i,hi(2)-1) &
                + 3.d0 * s(i,hi(2)-2) ) * EIGHTH
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! nothing to do - these ghost cells are filled with either
       ! multifab_fill_boundary or multifab_fill_ghost_cells
    else
       call bl_error("physbc_2d: bc(2,2) not yet supported")
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,hi,ng,bc,time,dx,prob_lo,prob_hi,icomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    real(kind=dp_t), intent(in   ) :: time,dx(:),prob_lo(:),prob_hi(:)
    integer        , intent(in   ) :: icomp

    call bl_error("physbc_3d not written")

  end subroutine physbc_3d

end module multifab_physbc_module
