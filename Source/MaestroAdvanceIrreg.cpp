
#include <Maestro.H>

using namespace amrex;

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStepIrreg (bool is_initIter) {

    // timer for profiling
    BL_PROFILE_VAR("Maestro::AdvanceTimeStepIrreg()",AdvanceTimeStepIrreg);

    // cell-centered MultiFabs needed within the AdvanceTimeStep routine
    Vector<MultiFab>      rhohalf(finest_level+1);
    Vector<MultiFab>       macrhs(finest_level+1);
    Vector<MultiFab>       macphi(finest_level+1);
    Vector<MultiFab>     S_cc_nph(finest_level+1);
    Vector<MultiFab> rho_omegadot(finest_level+1);
    Vector<MultiFab>     thermal1(finest_level+1);
    Vector<MultiFab>     thermal2(finest_level+1);
    Vector<MultiFab>     rho_Hnuc(finest_level+1);
    Vector<MultiFab>     rho_Hext(finest_level+1);
    Vector<MultiFab>           s1(finest_level+1);
    Vector<MultiFab>           s2(finest_level+1);
    Vector<MultiFab>       s2star(finest_level+1);
    Vector<MultiFab>      p0_cart(finest_level+1);
    Vector<MultiFab> delta_p_term(finest_level+1);
    Vector<MultiFab>       Tcoeff(finest_level+1);
    Vector<MultiFab>      hcoeff1(finest_level+1);
    Vector<MultiFab>     Xkcoeff1(finest_level+1);
    Vector<MultiFab>      pcoeff1(finest_level+1);
    Vector<MultiFab>      hcoeff2(finest_level+1);
    Vector<MultiFab>     Xkcoeff2(finest_level+1);
    Vector<MultiFab>      pcoeff2(finest_level+1);
    Vector<MultiFab>   scal_force(finest_level+1);
    Vector<MultiFab>    delta_chi(finest_level+1);
    Vector<MultiFab>       sponge(finest_level+1);

    // face-centered in the dm-direction (planar only)
    Vector<MultiFab> etarhoflux_dummy(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  umac(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > sedge(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > sflux(finest_level+1);

    ////////////////////////
    // needed for spherical routines only

    // cell-centered
    Vector<MultiFab> w0_force_cart_dummy(finest_level+1);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac_dummy(finest_level+1);


    // end spherical-only MultiFabs
    ////////////////////////

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    Vector<Real> grav_cell_nph   ( (max_radial_level+1)*nr_fine );
    Vector<Real> rho0_nph        ( (max_radial_level+1)*nr_fine );
    Vector<Real> p0_nph          ( (max_radial_level+1)*nr_fine );
    Vector<Real> peosbar         ( (max_radial_level+1)*nr_fine );
    Vector<Real> w0_force_dummy  ( (max_radial_level+1)*nr_fine );
    Vector<Real> Sbar            ( (max_radial_level+1)*nr_fine );
    Vector<Real> beta0_nph       ( (max_radial_level+1)*nr_fine );
    Vector<Real> gamma1bar_nph   ( (max_radial_level+1)*nr_fine );

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    Vector<Real> rho0_pred_edge_dummy( (max_radial_level+1)*(nr_fine+1) );

    // make sure C++ is as efficient as possible with memory usage
    grav_cell_nph       .shrink_to_fit();
    rho0_nph            .shrink_to_fit();
    p0_nph              .shrink_to_fit();
    peosbar             .shrink_to_fit();
    w0_force_dummy      .shrink_to_fit();
    Sbar                .shrink_to_fit();
    beta0_nph           .shrink_to_fit();
    gamma1bar_nph       .shrink_to_fit();
    rho0_pred_edge_dummy.shrink_to_fit();

    int is_predictor;

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    Print() << "\nTimestep " << istep << " starts with TIME = " << t_old
            << " DT = " << dt << std::endl << std::endl;

    if (maestro_verbose > 0) {
        Print() << "Cell Count:" << std::endl;
        for (int lev=0; lev<=finest_level; ++lev) {
            Print() << "Level " << lev << ", " << CountCells(lev) << " cells" << std::endl;
        }
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        // cell-centered MultiFabs
        rhohalf     [lev].define(grids[lev], dmap[lev],       1,    1);
        macrhs      [lev].define(grids[lev], dmap[lev],       1,    0);
        macphi      [lev].define(grids[lev], dmap[lev],       1,    1);
        S_cc_nph    [lev].define(grids[lev], dmap[lev],       1,    0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec,    0);
        thermal1    [lev].define(grids[lev], dmap[lev],       1,    0);
        thermal2    [lev].define(grids[lev], dmap[lev],       1,    0);
        rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1,    0);
        rho_Hext    [lev].define(grids[lev], dmap[lev],       1,    0);
        s1          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
        s2          [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
        s2star      [lev].define(grids[lev], dmap[lev],   Nscal, ng_s);
        p0_cart     [lev].define(grids[lev], dmap[lev],       1,    0);
        delta_p_term[lev].define(grids[lev], dmap[lev],       1,    0);
        Tcoeff      [lev].define(grids[lev], dmap[lev],       1,    1);
        hcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
        Xkcoeff1    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
        pcoeff1     [lev].define(grids[lev], dmap[lev],       1,    1);
        hcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
        Xkcoeff2    [lev].define(grids[lev], dmap[lev], NumSpec,    1);
        pcoeff2     [lev].define(grids[lev], dmap[lev],       1,    1);
        scal_force  [lev].define(grids[lev], dmap[lev],   Nscal,    1);
        delta_chi   [lev].define(grids[lev], dmap[lev],       1,    0);
        sponge      [lev].define(grids[lev], dmap[lev],       1,    0);

        // face-centered in the dm-direction (planar only)
        AMREX_D_TERM(etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);,
                     etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);,
                     etarhoflux_dummy[lev].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1););

        // face-centered arrays of MultiFabs
        AMREX_D_TERM(umac [lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1,     1);,
                     umac [lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1,     1);,
                     umac [lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1,     1););
        AMREX_D_TERM(sedge[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0);,
                     sedge[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0);,
                     sedge[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0););
        AMREX_D_TERM(sflux[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], Nscal, 0);,
                     sflux[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], Nscal, 0);,
                     sflux[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], Nscal, 0););
    }

#if (AMREX_SPACEDIM == 3)
    for (int lev=0; lev<=finest_level; ++lev) {
        w0mac_dummy[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
        w0mac_dummy[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
        w0mac_dummy[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
    }
    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            w0_force_cart_dummy[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        }
    }
#endif

    
    // set etarhoflux_dummy to zero
    for (int lev=0; lev<=finest_level; ++lev) {
        etarhoflux_dummy[lev].setVal(0.);
    }

#if (AMREX_SPACEDIM == 3)
    // initialize MultiFabs and Vectors to ZERO
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            w0mac_dummy[lev][d].setVal(0.);
        }
    }
    if (spherical == 1) {
        for (int lev=0; lev<=finest_level; ++lev) {
            w0_force_cart_dummy[lev].setVal(0.);
        }
    }
#endif

    std::fill(Sbar.begin(), Sbar.end(), 0.);
    
    // set dummy variables to zero
    std::fill(w0_force_dummy.begin(), w0_force_dummy.end(), 0.);
    std::fill(rho0_pred_edge_dummy.begin(), rho0_pred_edge_dummy.end(), 0.);
    std::fill(w0.begin(), w0.end(), 0.);

    // make the sponge for all levels
    if (do_sponge) {
        init_sponge(rho0_old.dataPtr());
        MakeSponge(sponge);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 1 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 1 : react state >>>" << std::endl;
    }

    React(sold,s1,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,0.5*dt);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 2 -- define average expansion at time n+1/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 2 : compute provisional S >>>" << std::endl;
    }

    if (t_old == 0.) {
        // this is either a pressure iteration or the first time step
        // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::LinComb(S_cc_nph[lev],0.5,S_cc_old[lev],0,0.5,S_cc_new[lev],0,0,1,0);
        }
    }
    else {
        // set S_cc_nph = S_cc_old + (dt/2) * dSdt
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::LinComb(S_cc_nph[lev],1.0,S_cc_old[lev],0,0.5*dt,dSdt[lev],0,0,1,0);
        }
    }
    // no ghost cells for S_cc_nph
    AverageDown(S_cc_nph,0,1);

    // compute delta_p_term = peos_old - p0_old (for RHS of projections)
    if (dpdt_factor > 0.0) {
	// peos_old (delta_p_term) now holds the thermodynamic p computed from sold(rho,h,X)
	PfromRhoH(sold,sold,delta_p_term);

	// no need to compute peosbar, p0_minus_peosbar since make_w0 is not called
	
	// compute p0_cart from p0
	Put1dArrayOnCart(p0_old, p0_cart, 0, 0, bcs_f, 0);

	// compute delta_p_term = peos_old - p0_old
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Subtract(delta_p_term[lev],p0_cart[lev],0,0,1,0);
	}
    }
    else {
        // these should have no effect if dpdt_factor <= 0
        for (int lev=0; lev<=finest_level; ++lev) {
            delta_p_term[lev].setVal(0.);
        }
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 3 -- construct the advective velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 3 : create MAC velocities >>>" << std::endl;
    }

    // compute unprojected MAC velocities
    AdvancePremac(umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy);
    
    for (int lev=0; lev<=finest_level; ++lev) {
        delta_chi[lev].setVal(0.);
        macphi   [lev].setVal(0.);
    }

    // Sbar = (1 / gamma1bar * p0) * dp/dt
    for (int i=0; i<Sbar.size(); ++i) {
	Sbar[i] = (p0_old[i] - p0_nm1[i])/(dtold*gamma1bar_old[i]*p0_old[i]);
    }
    
    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    is_predictor = 1;
    MakeRHCCforMacProj(macrhs,rho0_old,S_cc_nph,Sbar,beta0_old,gamma1bar_old,p0_old,
		       delta_p_term,delta_chi,is_predictor);

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac,macphi,macrhs,beta0_old,is_predictor);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 4 -- advect the full state through dt
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 4 : advect base >>>" << std::endl;
    }

    // no need to advect the base state density
    rho0_new = rho0_old;

    // thermal is the forcing for rhoh or temperature
    if (use_thermal_diffusion) {
	MakeThermalCoeffs(s1,Tcoeff,hcoeff1,Xkcoeff1,pcoeff1);

	MakeExplicitThermal(thermal1,s1,Tcoeff,hcoeff1,Xkcoeff1,pcoeff1,p0_old,
	                    temp_diffusion_formulation);
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            thermal1[lev].setVal(0.);
        }
    }

    // copy temperature from s1 into s2 for seeding eos calls
    // temperature will be overwritten later after enthalpy advance
    for (int lev=0; lev<=finest_level; ++lev) {
	s2[lev].setVal(0.);
        MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
    }

    if (maestro_verbose >= 1) {
        Print() << "            :  density_advance >>>" << std::endl;
        Print() << "            :   tracer_advance >>>" << std::endl;
    }

    // set sedge and sflux to zero
    for (int lev=0; lev<=finest_level; ++lev) {
	for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
	    sedge[lev][idim].setVal(0.);
	    sflux[lev][idim].setVal(0.);
	}
    }

    // advect rhoX, rho, and tracers
    DensityAdvance(1,s1,s2,sedge,sflux,scal_force,etarhoflux_dummy,umac,w0mac_dummy,rho0_pred_edge_dummy);

    // no need to compute etarho
    if (evolve_base_state) {
        // correct the base state density by "averaging"
	Average(s2, rho0_new, Rho);
	compute_cutoff_coords(rho0_new.dataPtr());
    }

    // update grav_cell_new
    if (evolve_base_state) {
	make_grav_cell(grav_cell_new.dataPtr(),
		       rho0_new.dataPtr(),
		       r_cc_loc.dataPtr(),
		       r_edge_loc.dataPtr());
    }
    else {
        grav_cell_new = grav_cell_old;
    }

    // base state pressure update
    if (evolve_base_state) {

	// set new p0 through HSE
	p0_new = p0_old;
	
	enforce_HSE(rho0_new.dataPtr(),
		    p0_new.dataPtr(),
		    grav_cell_new.dataPtr(),
		    r_cc_loc.dataPtr(),
		    r_edge_loc.dataPtr());
	
	// compute p0_nph
	for (int i=0; i<p0_nph.size(); ++i) {
	    p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
	}

	// no need for psi
    }
    else {
        p0_new = p0_old;
    }

    // base state enthalpy update
    if (evolve_base_state) {
	// compute rhoh0_old and rhoh0_new by "averaging"
	Average(s1, rhoh0_old, RhoH);
	Average(s2, rhoh0_new, RhoH);
    }
    else {
	rhoh0_new = rhoh0_old;
    }
    
    if (maestro_verbose >= 1) {
        Print() << "            : enthalpy_advance >>>" << std::endl;
    }
    
    EnthalpyAdvance(1,s1,s2,sedge,sflux,scal_force,umac,w0mac_dummy,thermal1);
    
    //////////////////////////////////////////////////////////////////////////////
    // STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 4a: thermal conduct >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
	ThermalConduct(s1,s2,hcoeff1,Xkcoeff1,pcoeff1,hcoeff1,Xkcoeff1,pcoeff1);
    }

    // pass temperature through for seeding the temperature update eos call
    // pi goes along for the ride
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
        MultiFab::Copy(s2[lev],s1[lev],  Pi,  Pi,1,ng_s);
    }

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s2,p0_new,0);
    }
    else {
        TfromRhoH(s2,p0_new);
    }

    if (use_thermal_diffusion) {
        // make a copy of s2star since these are needed to compute
        // coefficients in the call to thermal_conduct_full_alg
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Copy(s2star[lev],s2[lev],0,0,Nscal,ng_s);
	}
    }
    
    //////////////////////////////////////////////////////////////////////////////
    // STEP 5 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 5 : react state >>>" << std::endl;
    }

    React(s2,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_new,0.5*dt);

    if (evolve_base_state) {
        // compute beta0 and gamma1bar
        MakeGamma1bar(snew,gamma1bar_new,p0_new);
        make_beta0_irreg(beta0_new.dataPtr(), rho0_new.dataPtr(), p0_new.dataPtr(),
			 gamma1bar_new.dataPtr(), grav_cell_new.dataPtr(),
			 r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    }
    else {
        // Just pass beta0 and gamma1bar through if not evolving base state
        beta0_new = beta0_old;
        gamma1bar_new = gamma1bar_old;
    }

    for(int i=0; i<beta0_nph.size(); ++i) {
        beta0_nph[i] = 0.5*(beta0_old[i]+beta0_new[i]);
	gamma1bar_nph[i] = 0.5*(gamma1bar_old[i]+gamma1bar_new[i]);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 6 -- define a new average expansion rate at n+1/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 6 : make new S >>>" << std::endl;
    }

    if (evolve_base_state) {
        // reset cutoff coordinates to old time value
        compute_cutoff_coords(rho0_old.dataPtr());
    }

    if (use_thermal_diffusion) {
	MakeThermalCoeffs(snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2);

	MakeExplicitThermal(thermal2,snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
	                    temp_diffusion_formulation);
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            thermal2[lev].setVal(0.);
        }
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_new,snew,rho_omegadot,rho_Hnuc,rho_Hext,thermal2);

    // set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::LinComb(S_cc_nph[lev],0.5,S_cc_old[lev],0,0.5,S_cc_new[lev],0,0,1,0);
    }
    AverageDown(S_cc_nph,0,1);

    // and delta_p_term = peos_new - p0_new (for RHS of projection)
    if (dpdt_factor > 0.) {
	// peos_new now holds the thermodynamic p computed from snew(rho,h,X)
	PfromRhoH(snew,snew,delta_p_term);

	// no need to compute peosbar,p0_minus_peosbar since make_w0 is not called
       
	// compute p0_cart from p0
	Put1dArrayOnCart(p0_new, p0_cart, 0, 0, bcs_f, 0);

	// compute delta_p_term = peos_new - p0_new
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Subtract(delta_p_term[lev],p0_cart[lev],0,0,1,0);
	} 
    }
    else {
	// these should have no effect if dpdt_factor <= 0
	for (int lev=0; lev<=finest_level; ++lev) {
	    delta_p_term[lev].setVal(0.);
	}
    }


    //////////////////////////////////////////////////////////////////////////////
    // STEP 7 -- redo the construction of the advective velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 7 : create MAC velocities >>>" << std::endl;
    }

    // compute unprojected MAC velocities
    AdvancePremac(umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy);

    // compute Sbar
    for (int i=0; i<Sbar.size(); ++i) {
	Sbar[i] = (1.0/(gamma1bar_nph[i]*p0_new[i]))*(p0_new[i] - p0_old[i])/dt;
    }

    // compute RHS for MAC projection, beta0*(S_cc-Sbar) + beta0*delta_chi
    is_predictor = 0;
    MakeRHCCforMacProj(macrhs,rho0_new,S_cc_nph,Sbar,beta0_nph,gamma1bar_new,p0_new, 
		       delta_p_term,delta_chi,is_predictor);

    // MAC projection
    // includes spherical option in C++ function
    MacProj(umac,macphi,macrhs,beta0_nph,is_predictor);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8 -- advect the full state through dt
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 8 : advect base >>>" << std::endl;
    }

    // no need to advect the base state density

    // copy temperature from s1 into s2 for seeding eos calls
    // temperature will be overwritten later after enthalpy advance
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
    }

    if (maestro_verbose >= 1) {
        Print() << "            :  density_advance >>>" << std::endl;
        Print() << "            :   tracer_advance >>>" << std::endl;
    }

    // advect rhoX, rho, and tracers
    DensityAdvance(2,s1,s2,sedge,sflux,scal_force,etarhoflux_dummy,umac,w0mac_dummy,rho0_pred_edge_dummy);

    // no need to compute etarho
    if (evolve_base_state) {
        // correct the base state density by "averaging"
	Average(s2, rho0_new, Rho);
	compute_cutoff_coords(rho0_new.dataPtr());
    }


    // update grav_cell_new, rho0_nph, grav_cell_nph
    if (evolve_base_state) {
	make_grav_cell(grav_cell_new.dataPtr(),
		       rho0_new.dataPtr(),
		       r_cc_loc.dataPtr(),
		       r_edge_loc.dataPtr());
	
	for(int i=0; i<beta0_nph.size(); ++i) {
	    rho0_nph[i] = 0.5*(rho0_old[i]+rho0_new[i]);
	}
	
	make_grav_cell(grav_cell_nph.dataPtr(),
		       rho0_nph.dataPtr(),
		       r_cc_loc.dataPtr(),
		       r_edge_loc.dataPtr());
    } else {
        rho0_nph = rho0_old;
        grav_cell_nph = grav_cell_old;
    }

    // base state pressure update
    if (evolve_base_state) {

	// set new p0 through HSE
	p0_new = p0_old;
	
	enforce_HSE(rho0_new.dataPtr(),
		    p0_new.dataPtr(),
		    grav_cell_new.dataPtr(),
		    r_cc_loc.dataPtr(),
		    r_edge_loc.dataPtr());
	
	for (int i=0; i<p0_nph.size(); ++i) {
	    p0_nph[i] = 0.5*(p0_old[i] + p0_new[i]);
	}
	
	// no need for psi
    }

    // base state enthalpy averaging
    if (evolve_base_state) {
	Average(s2, rhoh0_new, RhoH);
    }
    
    // base state enthalpy update
    if (maestro_verbose >= 1) {
        Print() << "            : enthalpy_advance >>>" << std::endl;
    }

    EnthalpyAdvance(2,s1,s2,sedge,sflux,scal_force,umac,w0mac_dummy,thermal1);

    //////////////////////////////////////////////////////////////////////////////
    // STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 8a: thermal conduct >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
	MakeThermalCoeffs(s2star,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2);

	ThermalConduct(s1,s2,hcoeff1,Xkcoeff1,pcoeff1,hcoeff2,Xkcoeff2,pcoeff2);
    }

    // pass temperature through for seeding the temperature update eos call
    // pi goes along for the ride
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Copy(s2[lev],s1[lev],Temp,Temp,1,ng_s);
        MultiFab::Copy(s2[lev],s1[lev],  Pi,  Pi,1,ng_s);
    }

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s2,p0_new,0);
    }
    else {
        TfromRhoH(s2,p0_new);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 9 -- react the full state and then base state through dt/2
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 9 : react state >>>" << std::endl;
    }

    React(s2,snew,rho_Hext,rho_omegadot,rho_Hnuc,p0_new,0.5*dt);

    if (evolve_base_state) {
        //compute beta0 and gamma1bar
        MakeGamma1bar(snew,gamma1bar_new,p0_new);
        make_beta0_irreg(beta0_new.dataPtr(), rho0_new.dataPtr(), p0_new.dataPtr(),
			 gamma1bar_new.dataPtr(), grav_cell_new.dataPtr(),
			 r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    }

    for(int i=0; i<beta0_nph.size(); ++i) {
        beta0_nph[i] = 0.5*(beta0_old[i]+beta0_new[i]);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 10 -- compute S^{n+1} for the final projection
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 10: make new S >>>" << std::endl;
    }

    if (use_thermal_diffusion) {
	MakeThermalCoeffs(snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2);

	MakeExplicitThermal(thermal2,snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2,p0_new,
	                    temp_diffusion_formulation);
    }

    Make_S_cc(S_cc_new,snew,rho_omegadot,rho_Hnuc,rho_Hext,thermal2);

    // define dSdt = (S_cc_new - S_cc_old) / dt
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::LinComb(dSdt[lev],-1./dt,S_cc_old[lev],0,1./dt,S_cc_new[lev],0,0,1,0);
    }

    //////////////////////////////////////////////////////////////////////////////
    // STEP 11 -- update the velocity
    //////////////////////////////////////////////////////////////////////////////

    if (maestro_verbose >= 1) {
        Print() << "<<< STEP 11: update and project new velocity >>>" << std::endl;
    }

    // Define rho at half time using the new rho from Step 8
    FillPatch(0.5*(t_old+t_new), rhohalf, sold, snew, Rho, 0, 1, Rho, bcs_s);
       
    VelocityAdvance(rhohalf,umac,w0mac_dummy,w0_force_dummy,w0_force_cart_dummy,
		    rho0_nph,grav_cell_nph,sponge);


    int proj_type;

    // set Sbar to zero
    std::fill(Sbar.begin(), Sbar.end(), 0.);
    
    // Project the new velocity field
    if (is_initIter) {

        proj_type = pressure_iters_comp;

        // rhcc_for_nodalproj needs to contain
        // (beta0^nph S^1 - beta0^n S^0 ) / dt

        Vector<MultiFab> rhcc_for_nodalproj_old(finest_level+1);
        for (int lev=0; lev<=finest_level; ++lev) {
            rhcc_for_nodalproj_old[lev].define(grids[lev], dmap[lev], 1, 1);
            MultiFab::Copy(rhcc_for_nodalproj_old[lev], rhcc_for_nodalproj[lev], 0, 0, 1, 1);
        }

        MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_new,Sbar,beta0_nph);
        
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Subtract(rhcc_for_nodalproj[lev], rhcc_for_nodalproj_old[lev], 0, 0, 1, 1);
            rhcc_for_nodalproj[lev].mult(1./dt,0,1,1);
        }

    }
    else {

        proj_type = regular_timestep_comp;

        MakeRHCCforNodalProj(rhcc_for_nodalproj,S_cc_new,Sbar,beta0_nph);

	// compute delta_p_term = peos_new - p0_new (for RHS of projection)
        if (dpdt_factor > 0.) {
	    // peos_new now holds the thermodynamic p computed from snew(rho h X)
	    PfromRhoH(snew,snew,delta_p_term);
	    
	    // no need to compute peosbar, p0_minus_peosbar since make_w0 is not called

	    // compute peosbar_cart from peosbar
	    Put1dArrayOnCart(p0_new, p0_cart, 0, 0, bcs_f, 0);

	    // compute delta_p_term = peos_new - p0_new
	    for (int lev=0; lev<=finest_level; ++lev) {
		MultiFab::Subtract(delta_p_term[lev],p0_cart[lev],0,0,1,0);
	    }
	    
	    CorrectRHCCforNodalProj(rhcc_for_nodalproj,rho0_new,beta0_nph,gamma1bar_new,
				    p0_new,delta_p_term);
        }
    }

    // call nodal projection
    NodalProj(proj_type,rhcc_for_nodalproj);

    if (!is_initIter) {
	if (!fix_base_state) { 
	    // compute tempbar by "averaging"
	    Average(snew,tempbar,Temp);
	}

	// output any runtime diagnostics
	// pass in the new time value, time+dt
	// call diag(time+dt,dt,dx,snew,rho_Hnuc2,rho_Hext,thermal2,rho_omegadot2,&
        //          rho0_new,rhoh0_new,p0_new,tempbar, &
        //          gamma1bar_new,beta0_new, &
        //          unew,w0,normal, &
        //          mla,the_bc_tower)
    }

    Print() << "\nTimestep " << istep << " ends with TIME = " << t_new
            << " DT = " << dt << std::endl;

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;
	
    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to advance time step: " << end_total << '\n';
    }

}
