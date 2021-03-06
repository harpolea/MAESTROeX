module make_scal_force_module

  use base_state_geometry_module, only:  max_radial_level, nr_fine, dr, nr, base_cutoff_density_coord
  use meth_params_module, only: enthalpy_pred_type, use_exact_base_state
  use fill_3d_data_module, only: put_1d_array_on_cart_sphr

  implicit none

  private

contains

  subroutine mkrhohforce(lev, lo, hi, &
                         rhoh_force, f_lo, f_hi, &
#if (AMREX_SPACEDIM == 1)
                         umac, u_lo, u_hi, &
#elif (AMREX_SPACEDIM == 2)
                         vmac, v_lo, v_hi, &
#elif (AMREX_SPACEDIM == 3)
                         wmac, w_lo, w_hi, &
#endif
                         thermal, t_lo, t_hi, &
                         p0, rho0, grav, psi, &
                         is_prediction, add_thermal) &
                         bind(C,name="mkrhohforce")

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer         , intent(in   ) :: lev,lo(3),hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: rhoh_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
#if (AMREX_SPACEDIM == 1)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
#elif (AMREX_SPACEDIM == 2)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#elif (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(in   ) :: thermal(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    double precision, intent(in   ) ::   p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: grav(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::  psi(0:max_radial_level,0:nr_fine-1)
    integer         ,  intent(in   ) :: is_prediction, add_thermal

    ! Local variable
    double precision :: gradp0,veladv
    integer :: i,j,k,r

    integer, parameter :: predict_rhoh             = 0;
    integer, parameter :: predict_rhohprime        = 1;
    integer, parameter :: predict_h                = 2;
    integer, parameter :: predict_T_then_rhohprime = 3;
    integer, parameter :: predict_T_then_h         = 4;
    integer, parameter :: predict_hprime           = 5;
    integer, parameter :: predict_Tprime_then_h    = 6;

    !
    ! Add wtilde d(p0)/dr
    !
#if (AMREX_SPACEDIM == 1)
    j = 0
    k = 0
    do i = lo(1),hi(1)

       if (i .lt. base_cutoff_density_coord(lev)) then
          gradp0 = rho0(lev,i) * grav(lev,i)
       else if (i.eq.nr(lev)-1) then
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = ( p0(lev,i) - p0(lev,i-1) ) / dr(lev)
       else
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = ( p0(lev,i+1) - p0(lev,i) ) / dr(lev)
       end if

       veladv = 0.5d0*(umac(i)+umac(i+1))
       rhoh_force(i,j,k) = veladv * gradp0           
    enddo
#elif (AMREX_SPACEDIM == 2)
    k = 0
    do j = lo(2),hi(2)

       if (j .lt. base_cutoff_density_coord(lev)) then
          gradp0 = rho0(lev,j) * grav(lev,j)
       else if (j.eq.nr(lev)-1) then
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = ( p0(lev,j) - p0(lev,j-1) ) / dr(lev)
       else
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = ( p0(lev,j+1) - p0(lev,j) ) / dr(lev)
       end if

       do i = lo(1),hi(1)
          veladv = 0.5d0*(vmac(i,j,k)+vmac(i,j+1,k))
          rhoh_force(i,j,k) =  veladv * gradp0           
       end do
    enddo

#elif (AMREX_SPACEDIM == 3)
    do k = lo(3),hi(3)

       if (k .lt. base_cutoff_density_coord(lev)) then
          gradp0 = rho0(lev,k) * grav(lev,k)
       else if (k.eq.nr(lev)-1) then
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = ( p0(lev,k) - p0(lev,k-1) ) / dr(lev)
       else
          ! NOTE: this should be zero since p0 is constant up here
          gradp0 = ( p0(lev,k+1) - p0(lev,k) ) / dr(lev)
       end if

       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          veladv = 0.5d0*(wmac(i,j,k)+wmac(i,j,k+1))
          rhoh_force(i,j,k) = veladv * gradp0 
       end do
       end do
    enddo
#endif

    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h or rhoh
    if ((is_prediction .eq. 1 .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .eq. 1 .AND. enthalpy_pred_type == predict_rhoh) .OR. &
        (is_prediction .eq. 0)) then
       do k = lo(3),hi(3)
#if (AMREX_SPACEDIM == 3)
       r = k
#endif
       do j = lo(2),hi(2)
#if (AMREX_SPACEDIM == 2)
       r = j
#endif
       do i = lo(1),hi(1)
#if (AMREX_SPACEDIM == 1)
       r = i
#endif
          rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi(lev,r)
       end do
       end do
       enddo
    endif

    if (add_thermal .eq. 1) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
       end do
       end do
       end do
    end if

  end subroutine mkrhohforce

  subroutine mkrhohforce_sphr(lo, hi, &
                                rhoh_force, f_lo, f_hi, &
                                umac, u_lo, u_hi, &
                                vmac, v_lo, v_hi, &
                                wmac, w_lo, w_hi, &
                                thermal, t_lo, t_hi, &
                                p0_cart, p_lo, p_hi, &
                                p0macx, x_lo, x_hi, &
                                p0macy, y_lo, y_hi, &
                                p0macz, z_lo, z_hi, &
                                dx, psi, &
                                is_prediction, add_thermal, &
                                r_cc_loc, r_edge_loc, &
                                cc_to_r, ccr_lo, ccr_hi) &
             bind(C,name="mkrhohforce_sphr")

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer         , intent(in   ) :: lo(3),hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: rhoh_force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    double precision, intent(in   ) :: umac(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    double precision, intent(in   ) :: vmac(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: wmac(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    integer         , intent(in   ) :: t_lo(3), t_hi(3)
    double precision, intent(in   ) :: thermal(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer         , intent(in   ) :: p_lo(3), p_hi(3)
    double precision, intent(in   ) :: p0_cart(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) ::  p0macx(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) ::  p0macy(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    integer         , intent(in   ) :: z_lo(3), z_hi(3)
    double precision, intent(in   ) ::  p0macz(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: psi(0:max_radial_level,0:nr_fine-1)
    integer         , intent(in   ) :: is_prediction, add_thermal
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
                                               ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

    ! Local variable
    double precision, allocatable :: psi_cart(:,:,:,:)

    double precision :: divup, p0divu
    integer          :: i,j,k

    integer, parameter :: predict_rhoh             = 0;
    integer, parameter :: predict_rhohprime        = 1;
    integer, parameter :: predict_h                = 2;
    integer, parameter :: predict_T_then_rhohprime = 3;
    integer, parameter :: predict_T_then_h         = 4;
    integer, parameter :: predict_hprime           = 5;
    integer, parameter :: predict_Tprime_then_h    = 6;

    !
    ! Here we make u grad p = div (u p) - p div (u)
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k,divup,p0divu)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             divup = (umac(i+1,j,k) * p0macx(i+1,j,k) - umac(i,j,k) * p0macx(i,j,k)) / dx(1) + &
                     (vmac(i,j+1,k) * p0macy(i,j+1,k) - vmac(i,j,k) * p0macy(i,j,k)) / dx(2) + &
                     (wmac(i,j,k+1) * p0macz(i,j,k+1) - wmac(i,j,k) * p0macz(i,j,k)) / dx(3)

             p0divu = ( (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                        (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                        (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3) ) * p0_cart(i,j,k)

             rhoh_force(i,j,k) = divup - p0divu

          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! psi should always be in the force if we are doing the final update
    ! For prediction, it should not be in the force if we are predicting
    ! (rho h)', but should be there if we are predicting h or rhoh
    !
    if ((is_prediction .eq. 1 .AND. enthalpy_pred_type == predict_h) .OR. &
        (is_prediction .eq. 1 .AND. enthalpy_pred_type == predict_rhoh) .OR. &
        (is_prediction .eq. 0)) then

       allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

       call put_1d_array_on_cart_sphr(lo,hi,psi_cart,lo,hi,1,psi,dx,0,0, &
               r_cc_loc,r_edge_loc, cc_to_r,ccr_lo,ccr_hi)
          
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi_cart(i,j,k,1)
             enddo
          enddo
       enddo       
       !$OMP END PARALLEL DO
       deallocate(psi_cart)
    endif

    if (add_thermal .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine mkrhohforce_sphr

  subroutine modify_scal_force(lev, lo, hi, &
                               force, f_lo, f_hi, &
                               scal,  s_lo, s_hi, &
                               umac,  u_lo, u_hi, &
#if (AMREX_SPACEDIM >= 2)
                               vmac,  v_lo, v_hi, &
#if (AMREX_SPACEDIM == 3)
                               wmac,  w_lo, w_hi, &
#endif
#endif
                               s0, s0_edge, w0, dx, do_fullform) &
                               bind(C,name="modify_scal_force")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
#if (AMREX_SPACEDIM >= 2)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
#endif
#endif
    double precision, intent(inout) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
#if (AMREX_SPACEDIM >= 2)
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#endif
#endif
    double precision, intent(in   ) :: s0     (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: s0_edge(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: w0     (0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: do_fullform
    
    ! local
    integer :: i,j,k,r
    double precision :: divu,divs0u
    
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

       ! umac does not contain w0
#if (AMREX_SPACEDIM == 1)
       r = i
       divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1)
#elif (AMREX_SPACEDIM == 2)
       r = j
       divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
             +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2)
#elif (AMREX_SPACEDIM == 3)
       r = k
       divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
             +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
             +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
#endif

       ! add w0 contribution
       divu = divu + (w0(lev,r+1)-w0(lev,r))/dx(AMREX_SPACEDIM)

       if (do_fullform .eq. 1) then
          force(i,j,k) = force(i,j,k) - scal(i,j,k)*divu
       else

#if (AMREX_SPACEDIM == 1)
          divs0u = (umac(i+1,j,k) * s0_edge(lev,i+1) - &
                    umac(i  ,j,k) * s0_edge(lev,i  ) )/ dx(1)
#elif (AMREX_SPACEDIM == 2)
          divs0u = s0(lev,j)*(  umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                          +(vmac(i,j+1,k) * s0_edge(lev,j+1) - &
                            vmac(i,j  ,k) * s0_edge(lev,j  ) )/ dx(2)
#elif (AMREX_SPACEDIM == 3)
          divs0u = s0(lev,k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                          +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                          +(wmac(i,j,k+1) * s0_edge(lev,k+1) &
                          - wmac(i,j,k  ) * s0_edge(lev,k  ))/ dx(3)
#endif

          force(i,j,k) = force(i,j,k) - (scal(i,j,k)-s0(lev,r))*divu - divs0u 
       endif

    end do
    end do
    end do

  end subroutine modify_scal_force

  subroutine modify_scal_force_sphr(domlo, domhi, lo, hi, &
                                      force, f_lo, f_hi, &
                                      scal,  s_lo, s_hi, &
                                      umac,  u_lo, u_hi, &
                                      vmac,  v_lo, v_hi, &
                                      wmac,  w_lo, w_hi, &
                                      s0_cart,s0_lo, s0_hi, & 
                                      w0, dx, do_fullform, &
                                      r_cc_loc, r_edge_loc, &
                                      cc_to_r, ccr_lo, ccr_hi) &
             bind(C,name="modify_scal_force_sphr")

    integer         , intent(in   ) :: domlo(3), domhi(3), lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: u_lo(3), u_hi(3)
    integer         , intent(in   ) :: v_lo(3), v_hi(3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    integer         , intent(in   ) :: s0_lo(3), s0_hi(3) 
    double precision, intent(inout) :: force(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: scal (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision, intent(in   ) :: umac (u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    double precision, intent(in   ) :: vmac (v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    double precision, intent(in   ) :: wmac (w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    double precision, intent(in   ) :: s0_cart(s0_lo(1):s0_hi(1),s0_lo(2):s0_hi(2),s0_lo(3):s0_hi(3))
    double precision, intent(in   ) :: w0(0:max_radial_level,0:nr_fine)
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: do_fullform
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
                                               ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))
    
    ! Local variables
    integer :: i,j,k,r
    double precision :: divumac,divs0u
    double precision :: s0_xlo,s0_xhi
    double precision :: s0_ylo,s0_yhi
    double precision :: s0_zlo,s0_zhi
    
    double precision :: divu(0:max_radial_level,0:nr_fine-1)
    double precision, allocatable :: divu_cart(:,:,:,:)

    allocate(divu_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    !$OMP PARALLEL DO PRIVATE(r)
    if (use_exact_base_state) then
       divu(0,:) = 0.0d0
    else
       do r=0,nr_fine-1
          divu(0,r) = (r_edge_loc(0,r+1)**2 * w0(0,r+1) - &
                       r_edge_loc(0,r  )**2 * w0(0,r  ) ) / &
                      (dr(0)*r_cc_loc(0,r)**2)
       end do
    end if
    !$OMP END PARALLEL DO

    ! compute w0 contribution to divu
    call put_1d_array_on_cart_sphr(lo,hi,divu_cart,lo,hi,1,divu,dx,0,0,r_cc_loc,r_edge_loc, &
                                      cc_to_r,ccr_lo,ccr_hi)

    !$OMP PARALLEL DO PRIVATE(i,j,k,divumac,s0_xhi,s0_xlo,s0_yhi,s0_ylo) &
    !$OMP PRIVATE(s0_zhi,s0_zlo,divs0u)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             ! umac does not contain w0
             divumac = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                      +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                      +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)


             if (do_fullform .eq. 1) then

                force(i,j,k) = force(i,j,k) - scal(i,j,k)*(divumac+divu_cart(i,j,k,1))                 
             else
             
                if (i.lt.domhi(1)) then
                   s0_xhi = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i+1,j,k))
                else
                   s0_xhi = s0_cart(i,j,k)
                end if
                if (i.gt.domlo(1)) then
                   s0_xlo = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i-1,j,k))
                else
                   s0_xlo = s0_cart(i,j,k)
                end if
                
                if (j.lt.domhi(2)) then
                   s0_yhi = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j+1,k))
                else
                   s0_yhi = s0_cart(i,j,k)
                end if
                if (j.gt.domlo(2)) then
                   s0_ylo = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j-1,k))
                else
                   s0_ylo = s0_cart(i,j,k)
                end if
                
                if (k.lt.domhi(3)) then
                   s0_zhi = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j,k+1))
                else
                   s0_zhi = s0_cart(i,j,k)
                end if
                if (k.gt.domlo(3)) then
                   s0_zlo = 0.5d0 * (s0_cart(i,j,k) + s0_cart(i,j,k-1))
                else
                   s0_zlo = s0_cart(i,j,k)
                end if
                
                divs0u = (umac(i+1,j,k)*s0_xhi - umac(i,j,k)*s0_xlo)/dx(1) + &
                           (vmac(i,j+1,k)*s0_yhi - vmac(i,j,k)*s0_ylo)/dx(2) + &
                           (wmac(i,j,k+1)*s0_zhi - wmac(i,j,k)*s0_zlo)/dx(3)
             
                force(i,j,k) = force(i,j,k) - divs0u &
                     -(scal(i,j,k)-s0_cart(i,j,k))*(divumac+divu_cart(i,j,k,1)) 

             endif

          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    deallocate(divu_cart)

  end subroutine modify_scal_force_sphr

end module make_scal_force_module
