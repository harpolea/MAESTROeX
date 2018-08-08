! init_base_state is used to initialize the base state arrays from the
! model file.  The actual reading of the model file is handled by the
! model_parser_module in Util/
!
! Note: The initial base state quantities returned from this routine
! are only a temporary base state.  These quantities are mapped onto
! the full 2- or 3-d state in initscaldata.f90 and a new base state is
! created after initialization by averaging the density and calling
! enforce_HSE in initialize.f90.

module base_state_module

  use model_parser_module
  use eos_type_module
  use eos_module
  use extern_probin_module, only: eos_gamma
  use amrex_constants_module
  use simple_log_module
  use inlet_bc_module
  use fundamental_constants_module, only: Gconst
  use amrex_fort_module, only: amrex_spacedim
  use network, only: nspec
  use meth_params_module, only: nscal, model_file, spherical, base_cutoff_density, &
                                do_2d_planar_octant, do_planar_invsq_grav, rho_comp, &
                                rhoh_comp, spec_comp, temp_comp, grav_const, &
                                planar_invsq_mass, print_init_hse_diag, prob_lo, &
                                prob_hi, small_dens, small_temp, &
                                anelastic_cutoff, buoyancy_cutoff_factor
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use probdata_module, only: do_stratified, do_isentropic

  implicit none

  private

contains

  subroutine init_base_state(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init) &
       bind(C, name="init_base_state")

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)

    ! local variables
    integer         :: n,r,j
    real(kind=dp_t) :: H,z
    real(kind=dp_t) :: dens_zone, temp_zone, pres_zone
    real(kind=dp_t) :: xn_zone(nspec)

    type (eos_t) :: eos_state

    if (spherical .eq. 1) then
       call bl_error("ERROR: rt base_state is not valid for spherical")
    endif

    if (do_stratified) then

       ! use the EOS to make the state consistent
       dens_zone = 1.d-3
       pres_zone = 1.d6

       ! only initialize the first species
       xn_zone(:) = ZERO
       xn_zone(1) = 1.d0

       p0_init(:,0) = pres_zone

       ! H = pres_base / dens_base / abs(grav_const)
       H = 1.d6 / 1.d-3 / abs(grav_const)

       ! set an initial guess for the temperature -- this will be reset
       ! by the EOS
       temp_zone = 10.d0

        do n=0,max_radial_level

           do r=0,nr(n)-1

               z = (dble(r)+HALF) * dr(1)

               if (do_isentropic) then
                  dens_zone = 1.d-3*(grav_const*1.d-3*(eos_gamma - 1.0)*z/ &
                    (eos_gamma*1.d6) + 1.d0)**(1.d0/(eos_gamma - 1.d0))
               else
                  dens_zone = 1.d-3*exp(-z/H)
               end if

               s0_init(n,r, rho_comp) = dens_zone

               if (r.eq.0) then
                  p0_init(n,r) = p0_init(n,r) - &
                       dr(1) * HALF * s0_init(n,r,rho_comp) * &
                       abs(grav_const)
               else if (r.gt.0) then
                  p0_init(n,r) = p0_init(n,r-1) - &
                       dr(1) * HALF * (s0_init(n,r,rho_comp) + s0_init(n,r-1,rho_comp)) * &
                       abs(grav_const)
               end if

               pres_zone = p0_init(n,r)

               ! use the EOS to make the state consistent
               eos_state%rho   = dens_zone
               eos_state%T     = temp_zone
               eos_state%p     = pres_zone
               eos_state%xn(:) = xn_zone(:)

               ! (rho,p) --> T, h
               call eos(eos_input_rp, eos_state)

               s0_init(n,r, rho_comp) = dens_zone
               s0_init(n,r,rhoh_comp) = dens_zone*eos_state%h

               s0_init(n,r,spec_comp:spec_comp-1+nspec) = ZERO
               s0_init(n,r,spec_comp) = dens_zone
               s0_init(n,r,temp_comp) = eos_state%T

           end do

           ! initialize any inlet BC parameters
           call set_inlet_bcs()

        end do ! end loop over levels

    else

        do n=0,max_radial_level

            ! use the EOS to make the state consistent
            eos_state%T    = 10.d0
            eos_state%rho   = 1.d-3
            eos_state%p     = 1.d6
            eos_state%xn(:) = 1.d0

            ! (rho,p) --> T, h
            call eos(eos_input_rp, eos_state)

            s0_init(n,0:nr(n)-1, rho_comp) = eos_state%rho
            s0_init(n,0:nr(n)-1,rhoh_comp) = eos_state%rho * eos_state%h
            s0_init(n,0:nr(n)-1,spec_comp) = eos_state%rho
            s0_init(n,0:nr(n)-1,temp_comp) = eos_state%T

            p0_init(n,0:nr(n)-1) = eos_state%p

            ! initialize any inlet BC parameters
            call set_inlet_bcs()

        enddo

    endif

    ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    rho0 = s0_init(:,:,rho_comp)
    rhoh0 = s0_init(:,:,rhoh_comp)
    tempbar = s0_init(:,:,temp_comp)
    tempbar_init = s0_init(:,:,temp_comp)
    p0 = p0_init

  end subroutine init_base_state

  subroutine init_base_state_irreg(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init, &
                                     r_cc_loc, r_edge_loc) &
       bind(C, name="init_base_state_irreg")

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

  end subroutine  init_base_state_irreg

end module base_state_module
