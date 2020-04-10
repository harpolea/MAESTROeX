
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, const Array4<Real> scal, const Array4<Real> vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    const auto tileBox = mfi.tilebox();

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });

    const auto& s0_p = s0_init;

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i,j,k,Rho) = s0_p(lev,r,Rho);
        scal(i,j,k,RhoH) = s0_p(lev,r,RhoH);
        scal(i,j,k,Temp) = s0_p(lev,r,Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i,j,k,FirstSpec+comp) = s0_p(lev,r,FirstSpec+comp);
        }
        // initialize pi to zero for now
        scal(i,j,k,Pi) = 0.0;
    });    
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
			   MultiFab& scal, MultiFab& vel)
{
    Abort("ERROR: InitLeveLDataSphr not implemented.");
}
