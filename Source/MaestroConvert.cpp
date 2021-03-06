
#include <Maestro.H>

using namespace amrex;

void
Maestro::PutInPertForm(Vector<MultiFab>& scal, 
		       const Vector<Real>& s0, 
		       int comp, int bccomp, 
                       const Vector<BCRec>& bcs,
                       bool flag)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PutInPertForm()",PutInPertForm);

    // place 1d array onto a cartesian grid
    Vector<MultiFab> s0_cart(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	s0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    
    // s0 is not edge centered
    // note that bcs parameter is not used
    Put1dArrayOnCart(s0,s0_cart,0,0);

    if (flag) {
	for (int lev = 0; lev <= finest_level; ++lev) {
	    MultiFab::Subtract(scal[lev],s0_cart[lev],0,comp,1,0);
	}
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
	    MultiFab::Add(scal[lev],s0_cart[lev],0,comp,1,0);
	}
    }

    AverageDown(scal,comp,1);
    FillPatch(t_old,scal,scal,scal,comp,comp,1,bccomp,bcs);

}

void
Maestro::ConvertRhoXToX(Vector<MultiFab>& scal,
                        bool flag)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ConvertRhoXToX()",ConvertRhoXToX);

    if (flag) {
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                MultiFab::Divide(scal[lev],scal[lev],Rho,comp,1,0);
            }
        }
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                MultiFab::Multiply(scal[lev],scal[lev],Rho,comp,1,0);
            }
        }
    }

    // average down data and fill ghost cells
    AverageDown(scal,FirstSpec,NumSpec);
    if (flag) {
        FillPatch(t_old,scal,scal,scal,FirstSpec,FirstSpec,NumSpec,0,bcs_f);
    }
    else {
        FillPatch(t_old,scal,scal,scal,FirstSpec,FirstSpec,NumSpec,FirstSpec,bcs_s);
    }


    
}

void
Maestro::ConvertRhoHToH(Vector<MultiFab>& scal,
                        bool flag)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ConvertRhoHToH()",ConvertRhoHToH);

    if (flag) {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Divide(scal[lev],scal[lev],Rho,RhoH,1,0);
        }
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Multiply(scal[lev],scal[lev],Rho,RhoH,1,0);
        }
    }

    // average down data and fill ghost cells
    AverageDown(scal,RhoH,1);
    if (flag) {
        FillPatch(t_old,scal,scal,scal,RhoH,RhoH,1,0,bcs_f);
    }
    else {
        FillPatch(t_old,scal,scal,scal,RhoH,RhoH,1,RhoH,bcs_s);
    }
}
