
#include <Maestro.H>

using namespace amrex;

// check to see if we need to regrid, then regrid
void
Maestro::Regrid ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Regrid()",Regrid);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    Vector<Real> rho0_temp ( (max_radial_level+1)*nr_fine );
    rho0_temp.shrink_to_fit();
    
    if (spherical == 0) {
        finest_radial_level = finest_level;
	
        // look at MAESTRO/Source/varden.f90:750-1060
	// Regrid psi, etarho_cc, etarho_ec, and w0.
	// We do not regrid these in spherical since the base state array only
	// contains 1 level of refinement.
	// We do not regrid these if evolve_base_state=F since they are
	// identically zero, set this way in initialize.
        if (evolve_base_state) {

	    // We must regrid psi, etarho_cc, etarho_ec, and w0
	    // before we call init_multilevel or else we lose
	    // track of where the old valid data was 
	    
            // FIXME: may need if statement for irregularly-spaced base states

	    regrid_base_state_cc(psi.dataPtr());
	    regrid_base_state_cc(etarho_cc.dataPtr());
	    regrid_base_state_edge(etarho_ec.dataPtr());

        } else {
	    // evolve_base_state == F and spherical == 0

	    // Here we want to fill in the rho0 array so there is
	    // valid data in any new grid locations that are created
	    // during the regrid.
	    for (int i=0; i<rho0_old.size(); ++i) {
		rho0_temp[i] = rho0_old[i];
	    }
	    
	    // We will copy rho0_temp back into the rho0 array after we regrid.
	    regrid_base_state_cc(rho0_temp.dataPtr());
	    
        }
	
	// regardless of evolve_base_state, if new grids were
	// created, we need to initialize tempbar_init there, in
	// case drive_initial_convection = T
	regrid_base_state_cc(tempbar_init.dataPtr());
    }

    // regrid could add newly refine levels (if finest_level < max_level)
    // so we save the previous finest level index
    regrid(0, t_new);
    
    // Redefine numdisjointchunks, r_start_coord, r_end_coord
    TagArray();
    init_multilevel(tag_array.dataPtr(),&finest_level);
	
    if (spherical == 1) {
        MakeNormal();
    }

    if (evolve_base_state) {
        // force rho0 to be the average of rho
        Average(sold,rho0_old,Rho);
    } else {
	for (int i=0; i<rho0_old.size(); ++i) {
	    rho0_old[i] = rho0_temp[i];
	}
    }
    
    // compute cutoff coordinates
    compute_cutoff_coords(rho0_old.dataPtr());

    // make gravity
    make_grav_cell(grav_cell_old.dataPtr(),
                   rho0_old.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());

    // enforce HSE
    enforce_HSE(rho0_old.dataPtr(),
                p0_old.dataPtr(),
                grav_cell_old.dataPtr(),
                r_cc_loc.dataPtr(),
                r_edge_loc.dataPtr());

    
#ifdef AMREX_USE_CUDA
    // turn on GPU for eos calls
    Device::beginDeviceLaunchRegion();
#endif
    
    if (use_tfromp) {
        // compute full state T = T(rho,p0,X)
        TfromRhoP(sold,p0_old,0);
    } else {
        // compute full state T = T(rho,h,X)
        TfromRhoH(sold,p0_old);
    }

#ifdef AMREX_USE_CUDA
    // turn off GPU after eos calls
    Device::endDeviceLaunchRegion();
#endif
    
    // force tempbar to be the average of temp
    Average(sold,tempbar,Temp);

    // gamma1bar needs to be recomputed
    MakeGamma1bar(sold,gamma1bar_old,p0_old);

    // beta0_old needs to be recomputed
    if (use_exact_base_state) {
        make_beta0_irreg(beta0_old.dataPtr(), rho0_old.dataPtr(), p0_old.dataPtr(),
                         gamma1bar_old.dataPtr(), grav_cell_old.dataPtr(),
                         r_cc_loc.dataPtr(), r_edge_loc.dataPtr());
    } else {
        make_beta0(beta0_old.dataPtr(), rho0_old.dataPtr(), p0_old.dataPtr(),
                   gamma1bar_old.dataPtr(), grav_cell_old.dataPtr());
    }

    // wallclock time
    Real end_total = ParallelDescriptor::second() - strt_total;

    // print wallclock time
    ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
    if (maestro_verbose > 0) {
        Print() << "Time to regrid: " << end_total << '\n';
    }

}

// set tagging array to include buffer zones for multilevel
void
Maestro::TagArray ()
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::TagArray()",TagArray);

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    
    for (int lev=1; lev<=max_radial_level; ++lev) {
	
	const MultiFab& state = sold[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

        for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
        {
            const Box& tilebox  = mfi.tilebox();

            // retag refined cells to include cells in buffered regions
            retag_array(&tagval,&clearval, 
			ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
			&lev, tag_array.dataPtr());
        }
    }
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
Maestro::ErrorEst (int lev, TagBoxArray& tags, Real time, int ng)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ErrorEst()",ErrorEst);

    if (lev >= tag_err.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();

    const MultiFab& state = sold[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

        for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
        {
            const Box& tilebox  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];

            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            // set itags initially to 'untagged' everywhere
            // we define itags over the tilebox region
            tagfab.get_itags(itags, tilebox);

            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebox.loVect();
            const int*  thi     = tilebox.hiVect();

            // tag cells for refinement
	    // use row lev of tag_array to tag the correct elements
            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D(state[mfi]),
                        &tagval, &clearval,
                        ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
                        ZFILL(dx), &time,
			tag_err[lev].dataPtr(), &lev, tag_array.dataPtr());

	    // for planar refinement, we need to gather tagged entries in arrays
	    // from all processors and then re-tag tileboxes across each tagged
	    // height
	    if (spherical == 0) {
		ParallelDescriptor::ReduceIntMax(tag_array.dataPtr(),(max_radial_level+1)*nr_fine);

		tag_boxes(tptr, ARLIM_3D(tlo), ARLIM_3D(thi),
			  &tagval, &clearval,
			  ARLIM_3D(tilebox.loVect()), ARLIM_3D(tilebox.hiVect()),
			  ZFILL(dx), &time, &lev, tag_array.dataPtr());
	    }
	    
            //
            // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
            //
            tagfab.tags_and_untags(itags, tilebox);
        }
    }
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that existed before, using pre-existing fine and interpolated coarse data
// overrides the pure virtual function in AmrCore
void
Maestro::RemakeLevel (int lev, Real time, const BoxArray& ba,
                      const DistributionMapping& dm)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::RemakeLevel()",RemakeLevel);

    const int ng_s = snew[lev].nGrow();
    const int ng_u = unew[lev].nGrow();
    const int ng_S = S_cc_new[lev].nGrow();
    const int ng_g = gpi[lev].nGrow();
    const int ng_d = dSdt[lev].nGrow();

    MultiFab snew_state    (ba, dm,          Nscal, ng_s);
    MultiFab sold_state    (ba, dm,          Nscal, ng_s);
    MultiFab unew_state    (ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab uold_state    (ba, dm, AMREX_SPACEDIM, ng_u);
    MultiFab S_cc_new_state(ba, dm,              1, ng_S);
    MultiFab S_cc_old_state(ba, dm,              1, ng_S);
    MultiFab gpi_state     (ba, dm, AMREX_SPACEDIM, ng_g);
    MultiFab dSdt_state    (ba, dm,              1, ng_d);

    FillPatch(lev, time, snew_state, sold, snew, 0, 0, Nscal, 0, bcs_s);
    std::swap(snew_state, snew[lev]);
    std::swap(sold_state, sold[lev]);

    FillPatch(lev, time, unew_state, uold, unew, 0, 0, AMREX_SPACEDIM, 0, bcs_u);
    std::swap(unew_state, unew[lev]);
    std::swap(uold_state, uold[lev]);

    FillPatch(lev, time, S_cc_new_state, S_cc_old, S_cc_new, 0, 0, 1, 0, bcs_f);
    std::swap(S_cc_new_state, S_cc_new[lev]);
    std::swap(S_cc_old_state, S_cc_old[lev]);

    FillPatch(lev, time, gpi_state, gpi, gpi, 0, 0, AMREX_SPACEDIM, 0, bcs_f);
    std::swap(gpi_state, gpi[lev]);

    FillPatch(lev, time, dSdt_state, dSdt, dSdt, 0, 0, 1, 0, bcs_f);
    std::swap(dSdt_state, dSdt[lev]);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }
}

// within a call to AmrCore::regrid, this function fills in data at a level
// that did NOT exist before, using interpolated coarse data
// overrides the pure virtual function in AmrCore
void
Maestro::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                 const DistributionMapping& dm)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeNewLevelFromCoarse()",MakeNewLevelFromCoarse);

    snew[lev].define    (ba, dm,          Nscal, 0);
    sold[lev].define    (ba, dm,          Nscal, 0);
    unew[lev].define    (ba, dm, AMREX_SPACEDIM, 0);
    uold[lev].define    (ba, dm, AMREX_SPACEDIM, 0);
    S_cc_new[lev].define(ba, dm,              1, 0);
    S_cc_old[lev].define(ba, dm,              1, 0);
    gpi[lev].define     (ba, dm, AMREX_SPACEDIM, 0);
    dSdt[lev].define    (ba, dm,              1, 0);

    t_new = time;
    t_old = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    FillCoarsePatch(lev, time,     snew[lev],     sold,     snew, 0, 0,          Nscal, bcs_s);
    FillCoarsePatch(lev, time,     unew[lev],     uold,     unew, 0, 0, AMREX_SPACEDIM, bcs_u);
    FillCoarsePatch(lev, time, S_cc_new[lev], S_cc_old, S_cc_new, 0, 0,              1, bcs_f);
    FillCoarsePatch(lev, time,      gpi[lev],      gpi,      gpi, 0, 0, AMREX_SPACEDIM, bcs_f);
    FillCoarsePatch(lev, time,     dSdt[lev],     dSdt,     dSdt, 0, 0,              1, bcs_f);
}

// within a call to AmrCore::regrid, this function deletes all data
// at a level of refinement that is no longer needed
// overrides the pure virtual function in AmrCore
void
Maestro::ClearLevel (int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ClearLevel()",ClearLevel);

    sold[lev].clear();
    snew[lev].clear();

    uold[lev].clear();
    unew[lev].clear();

    S_cc_old[lev].clear();
    S_cc_new[lev].clear();

    gpi[lev].clear();

    dSdt[lev].clear();

    flux_reg_s[lev].reset(nullptr);
    flux_reg_u[lev].reset(nullptr);
}
