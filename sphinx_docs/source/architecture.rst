**********************
MAESTROeX Architecture
**********************

MAESTRO Basics
==============

MAESTRO models problems that are in tight hydrostatic equilibrium.
The fluid state is decomposed into a 1D radial base state that
describes the hydrostatic structure of the star or atmosphere, and a
2- or 3D Cartesian full state, that captures the departures from
hydrostatic equilibrium. Two basic geometries are allowed. A
*plane-parallel* geometry assumes that the domain is thin compared to
the radius of curvature of the star, and therefore the 1D base state
is perfectly aligned with the Cartesian state. A *spherical*
geometry is for modeling an entire star [1]_. Here, the 1D base state is not
aligned with the Cartesian state. Figure \ `[fig:base_state] <#fig:base_state>`__ shows
these geometries.

.. raw:: latex

   \centering

|[fig:base_state] MAESTRO geometries, showing both the
1D base state and the full Cartesian state. (Left) For multi-level
problems in planar geometry, we force a direct alignment between the
radial array cell centers and the Cartesian grid cell centers by
allowing the radial base state spacing to change with space and
time. (Right) For multi-level problems in spherical geometry, since
there is no direct alignment between the radial array cell centers
and the Cartesian grid cell centers, we choose to fix the radial
base state spacing across levels. Figure taken
from :raw-latex:`\cite{multilevel}`.|
|[fig:base_state] MAESTRO geometries, showing both the
1D base state and the full Cartesian state. (Left) For multi-level
problems in planar geometry, we force a direct alignment between the
radial array cell centers and the Cartesian grid cell centers by
allowing the radial base state spacing to change with space and
time. (Right) For multi-level problems in spherical geometry, since
there is no direct alignment between the radial array cell centers
and the Cartesian grid cell centers, we choose to fix the radial
base state spacing across levels. Figure taken
from :raw-latex:`\cite{multilevel}`.|

MAESTRO can use adaptive mesh refinement to focus resolution on
complex regions of flow. For Cartesian/plane-parallel geometries, all
cells at the same height must be at the same level of refinement.
This restriction is to allow for the base state to directly align with
the Cartesian state everywhere. For spherical geometries, there is no
such restriction (again, see Figure \ `[fig:base_state] <#fig:base_state>`__).
The MAESTRO grids are managed by the AMReX library, which is
distributed separately.

The MAESTRO ‘Ecosystem’
=======================

.. raw:: latex

   \centering

.. figure:: \archfigpath/maestro_ecosystem2
   :alt: [fig:arch:eco] The basic
   MAESTRO ‘ecosystem’. Here we see the different packages that
   contribute to building the reacting_bubble problem in MAESTRO. The
   red directories are part of most standard MAESTRO build. The
   purple lines show the directories that are pulled in through
   various Makefile variables (AMREX_HOME, NETWORK_DIR,
   EOS_DIR, and CONDUCTIVITY_DIR).

   [fig:arch:eco] The basic
   MAESTRO ‘ecosystem’. Here we see the different packages that
   contribute to building the reacting_bubble problem in MAESTRO. The
   red directories are part of most standard MAESTRO build. The
   purple lines show the directories that are pulled in through
   various Makefile variables (AMREX_HOME, NETWORK_DIR,
   EOS_DIR, and CONDUCTIVITY_DIR).

Building MAESTRO requires both the MAESTRO-specific source
files (distributed in the MAESTRO/ directory), and the
AMReX library (distributed separately, consisting of the amrex/ directory).
AMReX provides both a C++ and a Fortran framework. MAESTRO only uses the Fortran portions of AMReX. Figure \ `[fig:arch:eco] <#fig:arch:eco>`__
shows the relationship between the different packages, highlighting
what goes into defining a specific MAESTRO problem.

Problems piece together various MAESTRO directories, choosing a
reaction network, equation of state, and conductivity routine to build
an executable. Briefly, the MAESTRO sub-directories are:

-  MAESTRO/

   The main MAESTRO algorithm directory.

   Important directories under MAESTRO/ include:

   -  Docs/

      Documentation describing the basic algorithm (including this
      document).

   -  Exec/

      The various problem-setups. Each problem in MAESTRO gets it own
      sub-directory under SCIENCE/, TEST_PROBLEMS/, or
      UNIT_TESTS/. The GNUmakefile in the problem directory
      includes the instructions on how to build the executable,
      including what modules in Microphysics/ are used. Any file that
      you place in a sub-directory here takes precedence over a file of
      the same name in MAESTRO/. This allows problems to have
      custom versions of the main MAESTRO routines (e.g. initial
      conditions via initdata.f90. See § \ `3.1 <#sec:makefile>`__ and
      Chapter \ `[ch:make] <#ch:make>`__ for details on the build system).

      -  SCIENCE/

         MAESTRO problem directories for science studies. These are
         the setups that have been used for science papers in the past,
         or are the basis for current science studies.

      -  TEST_PROBLEMS/

         MAESTRO problem directories for simple test problems that have
         been used while developing MAESTRO. Many of these problems
         have appeared in the various papers describing the
         MAESTRO algorithm.

      -  UNIT_TESTS/

         Special MAESTRO problem directories that test only a single
         component of the MAESTRO algorithm. These often have their
         own main drivers (varden.f90) that setup and initialize
         some data structures and then call only a few of the
         MAESTRO routines. See Chapter \ `[chapter:unit_tests] <#chapter:unit_tests>`__ for details.

   -  Microphysics/ [2]_

      The basic microphysics routines used by MAESTRO. These are organized
      into the following sub-directories.

      -  conductivity/

         Various routines for computing the thermal conductivity used in
         the thermal diffusion part of the algorithm.

      -  EOS/

         The gamma_law_general/.

      -  networks/

         The basic general_null network that defines arbitrary
         non-reacting species.

   -  Source/

      The main MAESTRO source. Here you will find the driver routine,
      the advection routines, etc. All MAESTRO problems will compile
      this source.

   -  Util/

      Various helper routines exist in this directory. Some of these
      are externally developed.

      -  BLAS/

         Linear algebra routines.

      -  initial_models/

         Simple routines for generating toy initial models in hydrostatic equilibrium.

      -  model_parser/

         A simple Fortran module for reading in 1D initial model files.
         This is used by the initialization routines to get the initial
         model data.

      -  random/

         A random number generator.

      -  VODE/

         The VODE :raw-latex:`\cite{vode}` package for integrating ODEs. At the
         moment, this is used for integrating various reaction networks.

.. raw:: latex

   \centering

.. figure:: \archfigpath/amrex_directory2
   :alt: [fig:arch:amrex] The
   basic AMReX directory structure. The directories used by
   MAESTRO are indicated in red.

   [fig:arch:amrex] The
   basic AMReX directory structure. The directories used by
   MAESTRO are indicated in red.

The AMReX directory structure is shown in
Figure \ `[fig:arch:amrex] <#fig:arch:amrex>`__. The subset of the directories that are
used by MAESTRO are:

-  Src/

   The main AMReX source directory. In MAESTRO, we only use the
   Fortran source files. The core directories are:

   -  F_BaseLib/

      The Fortran AMReX files. This is a library for describing
      meshes consisting of a union of boxes. The AMReX modules
      define the basic datatypes used in MAESTRO. AMReX also
      provides the routines that handle the parallelization and I/O.

   -  LinearSolvers/

      The AMReX linear solvers—these are used to solve elliptic
      problems in the MAESTRO algorithm.

      -  F_MG

         The Fortran multigrid solver, with support for both
         cell-centered and node-centered data.

-  Tools/

   Various tools used for building AMReX applications. Here we use:

   -  F_mk/

      The generic Makefiles that store the compilation flags for various
      platforms. Platform/compiler-specific options are stored in the
      comps/ sub-directory.

   -  F_scripts/

      Some simple scripts that are useful for building, running,
      maintaining MAESTRO.

Finally the amrex/Tools/Postprocessing/F_Src package provides simple
Fortran-based analysis routines (e.g. extract a line from a
multidimensional dataset) that operate on AMReX datasets. These are
described in § \ `[sec:analysis] <#sec:analysis>`__. Several sub-directories with
python-based routines are also here. These are described both in
§ \ `[sec:analysis] <#sec:analysis>`__ and § \ `[sec:vis:python] <#sec:vis:python>`__.

.. _sec:adding_problems:

Adding A New Problem
====================

Different MAESTRO problems are defined in sub-directories under
Exec/ in SCIENCE, TEST_PROBLEMS, or UNIT_TESTS.
To add a problem, start by creating a new sub-directory—this is
where you will compile your problem and store all the problem-specific
files.

The minimum requirement to define a new problem would be a
GNUmakefile which describes how to build the application and an
input file which lists the runtime parameters. The problem-specific
executable is built in the problem directory by typing make.
Source files are found automatically by searching the directories
listed in the GNUmakefile. Customized versions of any source
files placed in the problem-directory override those with the same
name found elsewhere. Any unique source files (and not simply a
custom version of a file found elsewhere) needs to be listed in a file
called GPackage.mak in the problem-directory (and this needs to
be told to the build system—see below).

.. _sec:makefile:

The GNUmakefile
---------------

A basic GNUmakefile begins with:

::

      NDEBUG := t
      MPI    :=
      OMP    :=

Here, NDEBUG is true if we are building an optimized executable.
Otherwise, the debug version is built—this typically uses less
optimization and adds various runtime checks through compiler flags.
MPI and OMP are set to true if we want to use either MPI
or OpenMP for parallelization. If MPI := t, you will need to
have the MPI libraries installed, and their location may need to be
specified in MAESTRO/mk/GMakeMPI.mak.

The next line sets the compiler to be used for compilation:

::

      COMP := gfortran

The MAESTRO build system knows what options to use for various
compiler families. The COMP flag specifies which compiler to
use. Allowed values include Intel, gfortran, PGI,
PathScale, and Cray. The specific details of these
choices are defined in the MAESTRO/mk/comps/ directory.

MKVERBOSE set to true will echo the build commands to the
terminal as the are executed.

::

      MKVERBOSE := t

The next line defines where the top of the MAESTRO source tree is located.

::

      MAESTRO_TOP_DIR := ../../..

A MAESTRO application is built from several packages (the
multigrid solver, an EOS, a reaction network, etc.). The core
MAESTRO packages are always included, so a problem only needs
to define the EOS, reaction network, and conductivities to
use, as well as any extra, problem-specific files.

::

    EOS_DIR := helmholtz   
    CONDUCTIVITY_DIR := constant
    NETWORK_DIR := ignition_simple

    EXTRA_DIR := Util/random

Note that the microphysics packages are listed simply by the name of
the directory containing the specific implementation (e.g. helmholtz).
By default, the build system will look in Microphysics/EOS/ for
the EOS, Microphysics/conductivity/ for the conductivity routine,
and Microphysics/networks/ for the reaction network. To
override this default search path, you can set EOS_TOP_DIR,
CONDUCTIVITY_TOP_DIR, and NETWORK_TOP_DIR respectively.

Generally, one does not need to include the problem directory itself
in EXTRA_DIR, unless there are unique source files found there,
described in a GPackage.mak file. These variables are
interpreted by the GMaestro.mak file and used to build a master
list of packages called Fmdirs. The build system will attempt
to build all of the files listed in the various GPackage.mak
files found in the Fmdirs directories. Furthermore,
Fmdirs will be will be added to the make VPATH, which
is the list of directories to search for source files. The problem
directory will always be put first in the VPATH, so any source
files placed there override those with the same name found elsewhere
in the source tree.

Some packages (for instance, the helmholtz
EOS) require Fortran include files. The Fmincludes variable
lists all those directories that contain include files that are
inserted into the Fortran source at compile time via the include
statement. Presently, the only instance of this is with the Helmholtz
general equation of state found in Microphysics/EOS/helmholtz/. This is
automatically handled by the GMaestro.mak instructions.

Runtime parameters listed in the MAESTRO/_parameters file are
parsed at compile time and the file probin.f90 is written and
compiled. This is a Fortran module that holds the values of the
runtime parameters and makes them available to any routine. By
default, the build system looks for a file called \_parameters
in the problem directory and adds those parameters along with the
master list of MAESTRO parameters (MAESTRO/_parameters) to
the probin_module.

The final line in the GNUmakefile includes the rules to actually
build the executable.

::

      include $(MAESTRO_TOP_DIR)/GMaestro.mak

Handling Problem-Specific Source Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As mentioned above, any source files placed in the problem directory
override a files with the same name found elsewhere in the source
tree. This allows you to create a problem-specific version of any
routine. Source files that are unique to this problem (i.e. there is
no file with the same name elsewhere in the source tree) need to be
listed in a file GPackage.mak in the problem directory, and
the problem-directory needs to be explicitly listed in the EXTRA_DIR
list in the GNUmakefile.

.. _sec:def_runtime_param:

Defining Runtime Parameters
---------------------------

The runtime parameters for the core MAESTRO algorithm are listed in
MAESTRO/_parameters. That file is parsed at compile-time by
the MAESTRO/write_probin.py script (along with any
problem-specific parameters). The script outputs the probin.f90
source file. Each line in the \_parameters file has the form:
10em *data-type* 10em *value*
where *parameter* is the name of the runtime parameter,
*data-type* is one of {character, real,
integer, logical}, and the *value* specifies the default
value for the runtime parameter. Comments are indicated by a ‘
#’ character and are used to produce documentation about the
available runtime parameters. For the documentation, runtime parameters are grouped together
in the \_parameters file into categories. The category headings
are defined by comments in the \_parameters file and any comments
following that heading are placed into that category. The documentation
(Chapter `[ch:parameters] <#ch:parameters>`__) is produced by the script
MAESTRO/docs/runtime_parameters/rp.py.

At runtime, the default values for the parameters can be overridden
either through the inputs file (by adding a line of the form:
parameter = value) or through a command-line argument (taking the
form: –parameter value). The probin_module makes the
values of the runtime parameters available to the various functions
in the code (see § \ `6.7 <#sec:probin>`__).

Problem-specific runtime parameters should be defined in the
problem-directory in a file called \_parameters. This file will
be automatically found at compile time.

.. _sec:initial_models:

Preparing the Initial Model
---------------------------

MAESTRO models subsonic, non-hydrostatic flows as deviations from
a background state in hydrostatic equilibrium.
The solution in MAESTRO is broken up into a 1D base state and the 2-
or 3D full state. The job of the 1D base state in the algorithm is
to represent the hydrostatic structure. The full, Cartesian state
carries the departures from hydrostatic equilibrium. The underlying
formulation of the low Mach number equations assumes that the base
state is in hydrostatic equilibrium. At the start of a simulation,
the initial model is read in and taken as the base state. Therefore,
any initial model needs to already be in hydrostatic equilibrium.

The routines in Util/initial_models/ prepare an initial model
for MAESTRO. In general, there are two different proceduces that are
needed. The first type modify an existing 1D initial model produced
somewhere else (e.g. a 1D stellar evolution code), and map it onto a
uniform grid, at the desired resolution, using the equation of state
in MAESTRO, and using MAESTRO’s discretization of hydrostatic
equilibrium. The second type generate the initial model internally,
by integrating the condition of hydrostatic equilibrium together with
a simplifying assumption on the energy (e.g. isothermal or
isentropic). In both cases hydrostatic equilibrium is enforced as:

.. math::

   \frac{p_{i+1} - p_i}{\Delta r} = \frac{1}{2} (\rho_i + \rho_{i+1})
   g_{i+1/2}

Here, :math:`g_{i+1/2}` is the edge-centered gravitational acceleration.

The toy_atm example provides a simple approximation for a thin
(plane-parallel) convectively-unstable accreted layer on the surface
of a star. This can be used as the starting point for a more complex
model.

MAESTRO initial models are read in by the Util/model_parser
routines. This expects the initial model to contain a header giving
the number of variables and their names, followed by rows of data
giving the coordinate and data values at that coordinate. The initial
model should contain the same species data (in the form of mass fractions) as
defined in the network module used by the MAESTRO problem.

Full details on which initial model routine matches each problem and
how the initial models are used to initialize the full state data can
be found in § \ `[sec:initial_models_main] <#sec:initial_models_main>`__.

Customizing the Initialization
------------------------------

The best way to customize the initialization (e.g. add perturbations)
is to copy from one of the existing problems. The file initveldata.f90 controls the velocity field initialization and initscaldata.f90 controls the initialization of the scalars
(:math:`\rho`, :math:`\rho X_k`, :math:`\rho h`). The reacting_bubble problem is a good
starting point for plane-parallel and wdconvect is a good
starting point for full stars.

AMReX Data Structures
=====================

MAESTRO’s gridding is handled by the AMReX library, which
contains the most fundamental objects used to construct parallel
block-structured AMR applications—different
regions of the domain can have different spatial resolutions.
At each level of refinement, the region covered by that level is divided
into grids, or boxes. The entire computational domain is covered by
the coarsest (base) level of refinement, often called level :math:`\ell=0`, either by one
grid or divided into many grids.
Higher levels of refinement have cells that are finer by a “refinement ratio”
(typically 2). The grids are properly nested in the sense that the union
of grids at level :math:`\ell+1` is contained in the union of grids at level :math:`\ell`.
Furthermore, the containment is strict in the sense that, except at physical
boundaries, the level :math:`\ell` grids are large enough to guarantee that there is
a border at least :math:`n_{\rm buffer}` level :math:`\ell` cells wide surrounding each level
:math:`\ell +1` grid (grids at all levels are allowed to extend to the physical
boundaries so the proper nesting is not strict there).
For parallel computations, the boxes are spread across processors, in
a fashion designed to put roughly equal amounts of work on each
processor (load balancing).

.. raw:: latex

   \centering

.. figure:: \archfigpath/data_loc2
   :alt: [fig:dataloc] Some of the different data-centerings:
   (a) cell-centered, (b) nodal in the :math:`x`-direction, and (c) nodal in
   both the :math:`x`- and :math:`y`-directions. Note that for nodal data, the
   integer index corresponds to the lower boundary in that direction.
   In each of these centerings, the red point has the same indices: (1,2).
   Not shown is the case where data is nodal in the :math:`y`-direction only.
   :width: 6.5in

   [fig:dataloc] Some of the different data-centerings:
   (a) cell-centered, (b) nodal in the :math:`x`-direction, and (c) nodal in
   both the :math:`x`- and :math:`y`-directions. Note that for nodal data, the
   integer index corresponds to the lower boundary in that direction.
   In each of these centerings, the red point has the same indices: (1,2).
   Not shown is the case where data is nodal in the :math:`y`-direction only.

On a grid, the data can be stored at cell-centers, on a face/edge, or
on the corners. In AMReX, data that is on an edge is termed ‘nodal’
in that direction (see Figure \ `[fig:dataloc] <#fig:dataloc>`__). Data that is on the
corners is nodal in all spatial directions. In MAESTRO, the state
data (density, enthalpy, velocity, :math:`\ldots`) is generally
cell-centered. Fluxes are nodal in the direction they represent.
A few quantities are nodal in all directions (e.g. :math:`\phi` used in
the final velocity projection).

To simplify the description of the underlying AMR grid, AMReX provides a number of Fortran types. We briefly summarize the major
data types below. A more extensive introduction to AMReX is
provided by the AMReX User’s Guide, distributed with the library.

box
---

A box is simply a rectangular domain in space. Note that boxes
do not hold the state data themselves. A box has a lo
and hi index in each coordinate direction that gives the
location of the lower-left and upper-right corner with respect to
a global index space.

.. raw:: latex

   \centering

.. figure:: \archfigpath/index_grid2
   :alt: [fig:boxes] Three boxes that comprise a single level. At this
   resolution, the domain is 20\ :math:`\times`\ 18 zones. Note that the
   indexing in AMReX starts with :math:`0`.
   :width: 4in

   [fig:boxes] Three boxes that comprise a single level. At this
   resolution, the domain is 20\ :math:`\times`\ 18 zones. Note that the
   indexing in AMReX starts with :math:`0`.

The computational domain is divided into boxes. The collection of
boxes with the same resolution comprise a level.
Figure \ `[fig:boxes] <#fig:boxes>`__ shows three boxes in the same level of
refinement. The position of the boxes is with respect to the global
index space at that level. For example, box 1 in the figure has
lo = (3,7) and hi = (9,12). Note that the global indexing
is 0-based.

The global index space covers the entire domain at a given resolution.
For a simulation setup with n_cellx = 32 and n_celly =
32, the coarsest level (level 1) has :math:`32 \times 32` zones, and the
global index space will run from :math:`0, \ldots, 31` in each coordinate
direction. Level 2 will have a global index space running from :math:`0,
\ldots, 63` in each coordinate direction (corresponding to :math:`64 \times
64` zones if fully refined), and level 3 will have a global index
space running from :math:`0, \ldots, 127` in each coordinate direction
(corresponding to :math:`128\times 128` zones if fully refined).

Common Operations on a box
~~~~~~~~~~~~~~~~~~~~~~~~~~

A box declared as:

::

      type(box) :: mybox

The upper and lower bounds of the box (in terms of the global
index space) are found via:

-  lo = lwb(mybox) returns an array, lo(dm), with
   the box lower bounds

-  hi = upb(mybox) returns an array, hi(dm), with
   the box upper bounds

boxarray and ml_boxarray
------------------------

A boxarray is an array of boxes. A ml_boxarray is a collection of
boxarrays at different levels of refinement.

layout and ml_layout
--------------------

A layout is basically a boxarray that knows information about other
boxes, or box “connectivity.” It contains additional information
that is used in filling ghost cells from other fine grids or from
coarser grids. This information is stored as long as the layout
exists so that we don’t have to recompute intersections every time we
do some operation with two multifabs that have that layout, for
example.

By separating the layout from the actual data, we can allocate and
destroy data that lives on the grid as needed.

fab
---

A fab is a “Fortran Array Box”. It contains the state data in a
multidimensional array and several box-types to describe where in
the global index-space it lives:

::

      type fab
         ...
         type(box) :: bx
         type(box) :: pbx
         type(box) :: ibx
      end type fab

bx represents the box in the global index-space over which the
fab is defined, pbx represents the “physical” box in the
sense that it includes bx plus ghost cells, and ibx is the
same as bx unless the fab is nodal. As can be seen in
Figure \ `[fig:dataloc] <#fig:dataloc>`__, for the same grid nodal data requires one
more array element than cell-centered data. To address this ibx
is made by growing bx by one element along all nodal dimensions.

It’s important to note that all state data is stored in a
four-dimensional array *regardless of the problem’s
dimensionality*. The array is (nx,ny,nz,nc) in size, where
nc is the number of components, for instance representing
different fluid variables, and (nx,ny,nz) are the number of
cells in each respective spatial dimension. For 2D problems,
nz=1.

A fab would represent the data for a single box in the domain.
In MAESTRO, we don’t usually deal with fabs alone, but rather
we deal with multifabs, described next.

multifab
--------

A multifab is a collection of fabs at the same level of
refinement. This is the primary data structure that MAESTRO routine operate on. A multilevel simulation stores the
data in an array of multifabs, where the array index refers
to the refinement level.

All fabs in a given multifab have the same number of ghost cells,
but different multifabs can have different numbers of ghost cells
(or no ghost cells).

Working with multifabs
~~~~~~~~~~~~~~~~~~~~~~

To build a multifab, we need to provide a layout, the number of
components to store in the multifab  and the number of ghostcells. In
MAESTRO  the hierarchy of grids will be described by a single
ml_layout. A multifab can be declared and built at any time in a
simulation using the ml_layout, thereby allocating space at every
grid location in the simulation. The sequence to build a multifab appears as

::

      type(multifab) :: mfab(nlevs)
      ...
      do n = 1, nlevs
         call multifab_build(mfab(n), mla%la(n), nc, ng)
      enddo

Here, nc is the number of components and ng is the number
of ghostcells. The multifab is built one level at a time, using the
layout for that level taken from the ml_layout, mla.

A common operation on a multifab is to initialize it to :math:`0`
everywhere. This can be done (level-by-level) as

::

    call setval(mfab(n), ZERO, all=.true.)

where ZERO is the constant 0.0 from bl_constants_module.

The procedure for accessing the data in each grid managed by the
multifab is shown in § \ `[sec:example] <#sec:example>`__. Subroutines to add,
multiply, or divide two multifabs exist, as do subroutines to copy
from one multifab to another—see
amrex/Src/F_BaseLib/multifab.f90 for the full list of
routines that work with multifabs.

When you are done working with a multifab, its memory can be freed by
calling multifab_destroy on the multifab.

bc_tower
--------

A bc_tower holds the information about what boundary conditions are
in effect for each variable in a
MAESTRO simulation. These are interpretted by the ghost cell filling
routines. See § \ `10 <#sec:arch:bcs>`__ for more detail.

MAESTRO Data Organization
=========================

The state of the star in MAESTRO is described by both a
multidimensional state and the 1D base state. The full
multidimensional state is stored in multifabs while the base state
is simply stored in Fortran arrays. Here we describe the
major MAESTRO data-structures.

‘s’ multifabs (fluid state)
---------------------------

The fluid state (density, enthalpy, species, temperature, and tracer)
are stored together in a cell-centered multi-component multifab,
typically named sold, s1, s2, or snew
(depending on which time-level it represents). The enthalpy is stored
as :math:`(\rho h)`, and the species are stored as partial-densities :math:`(\rho
X_k)`. The tracer component is not used at present time, but can
describe an arbitrary advected quantity.

Individual state variables should be indexed using the integer keys
provided by the variables module (see §
`6.8 <#sec:variables_module>`__). For example, the integer rho_comp
will always refer to the density component of the state.

Note: the pressure is not carried as part of the ‘s’ multifabs.

‘u’ multifabs (fluid velocity)
------------------------------

The fluid velocity at time-levels :math:`n` and :math:`n+1` is stored in
a cell-centered multi-component multifab, typically named
uold or unew. Here the dm
components correspond to each coordinate direction.

umac (the MAC velocity)
-----------------------

In creating the advective fluxes, we need the time-centered velocity
through the faces of the zone—the :math:`x`-velocity on the :math:`x`-edges, the
:math:`y`-velocity on the :math:`y`-edges, etc. (see figure \ `[fig:mac] <#fig:mac>`__). This
type of velocity discretization is termed the MAC velocity (after the
“marker-and-cell” method for free boundaries in incompressible
flows :raw-latex:`\cite{harlowwelch:1965}`).

.. raw:: latex

   \centering

.. figure:: \archfigpath/mac2
   :alt: [fig:mac] The MAC grid for the velocity.
   Here the :math:`x`-velocity is on the :math:`x`-edges (shown as the
   blue points) and the :math:`y`-velocity is on the :math:`y`-edges
   (shown as the red points).
   :width: 2.5in

   [fig:mac] The MAC grid for the velocity.
   Here the :math:`x`-velocity is on the :math:`x`-edges (shown as the
   blue points) and the :math:`y`-velocity is on the :math:`y`-edges
   (shown as the red points).

|  

The MAC velocities are allocated at each level of refinement, n,
by making a multifab array where each of the dm components is
nodal in its respective direction:

::

      type(multifab) :: umac(nlevel,dm)

      do n=1,nlevel
         do comp=1,dm
            call multifab_build_edge(umac(n,comp), mla%la(n),1,1,comp)
         enddo
      enddo

Base State Arrays
-----------------

The base state is defined by :math:`\rho_0`, :math:`p_0`, and :math:`w_0`. There is no
base state composition. Other arrays are defined as needed, such as
:math:`h_0`, the base state enthalpy.

The base state arrays are 2-dimensional, with the first dimension
giving the level in the AMR hierarchy and the second the radial index
into the base state. For spherical geometries, the base state only
exists at a single level, so the first index will always be 1. The
radial index is 0-based, to be consistent with the indexing for the
Cartesian state data. For example, the base state density would be
dimensioned: rho0(nlevs,0:nr_fine-1). Here, nlevs is the
number of levels of refinement and nr_fine is the number of
cells in the radial direction at the finest level of refinement.

For multilevel, plane-parallel geometry, all grids at the same height
will have the same resolution so that the full state data is always
aligned with the base state (see Figure \ `[fig:base_state] <#fig:base_state>`__). Base
state data on coarse grids that are covered by fine grids is not
guaranteed to be valid.

For spherical problems, the base state resolution, :math:`\Delta r`, is
generally picked to be finer than the Cartesian grid resolution,
:math:`\Delta x`, i.e. \ :math:`\Delta r < \Delta x`. The ratio is controlled
by the parameter drdxfac.

Note there are no ghost cells for the base state outside of the
physical domain. For plane-parallel, multilevel simulations, there
are ghostcells at the jumps in refinement—these are filled by the
fill_code_base routine. The convention when dealing with the
base state is that we only access it inside of the valid physical
domain. Any multi-dimensional quantity that is derived using the base
state then has its ghost cells filled by the usually multifab ghost
cell routines.

MAESTRO Helper Modules
======================

A number of MAESTRO modules appear frequently throughout the source.
Below, we describe some of the more common functionality of the most
popular modules.

average_module
--------------

The average_module module provides a routine average that takes
a multilevel multifab array and averages the full Cartesian data
onto the 1D base state.

eos_module
----------

The eos_module provides the interface to the equation of
state to connect the state variables thermodynamically. It
gets the information about the fluid species from the network
module (for example, the atomic number, :math:`Z`, and atomic weight, :math:`A`,
of the nuclei).

Presently there is a single EOS that comes with MAESTRO, tt gamma_law_general,
but many more are available through the external Microphysics repo [3]_. The Microphysics EOSs share the same interface and can be compiled into MAESTRO directly.
Here are the more popular EOSs:

-  helmholtz represents a general stellar equation
   of state, consisting of nuclei (as an ideal gas), radiation,
   and electrons (with arbitrary degeneracy and degree of relativity).
   This equation of state is that described in :raw-latex:`\cite{timmes_eos}`.

   A runtime parameter, use_eos_coulomb, is defined in
   this EOS to enable/disable Coulomb corrections.

-  gamma_law_general assumes an ideal gas with a mixed
   composition and a constant ratio of specific heats, :math:`\gamma`:

   .. math:: p = \rho e (\gamma - 1) = \frac{\rho k_B T}{\mu m_p}

   where :math:`k_B` is Boltzmann’s constant and :math:`m_p` is the mass of the
   proton.
   The mean molecular weight, :math:`\mu`, is computed assuming
   electrically neutral atoms:

   .. math:: \mu = \left ( \sum_k \frac{X_k}{A_k} \right )^{-1}

   An option in the source code itself exists for treating the
   species as fully-ionized, but there is no runtime-parameter to
   make this switch.

-  multigamma is an ideal gas equation of state where each
   species can have a different value of :math:`\gamma`. This mainly affects
   how the internal energy is constructed as each species, represented
   with a mass fraction :math:`X_k` will have its contribution to the total
   specific internal energy take the form of :math:`e = p/\rho/(\gamma_k -                                               
     1)`. The main thermodynamic quantities take the form:

   .. math::

      \begin{aligned}
      p &= \frac{\rho k T}{m_u} \sum_k \frac{X_k}{A_k} \\
      e &= \frac{k T}{m_u} \sum_k \frac{1}{\gamma_k - 1} \frac{X_k}{A_k} \\
      h &= \frac{k T}{m_u} \sum_k \frac{\gamma_k}{\gamma_k - 1} \frac{X_k}{A_k}\end{aligned}

   We recognize that the usual astrophysical :math:`\bar{A}^{-1} = \sum_k                                                  
   X_k/A_k`, but now we have two other sums that involve different
   :math:`\gamma_k` weightings.

   The specific heats are constructed as usual,

   .. math::

      \begin{aligned}
      c_v &= \left . \frac{\partial e}{\partial T} \right |_\rho =
          \frac{k}{m_u} \sum_k \frac{1}{\gamma_k - 1} \frac{X_k}{A_k} \\
      c_p &= \left . \frac{\partial h}{\partial T} \right |_p =
          \frac{k}{m_u} \sum_k \frac{\gamma_k}{\gamma_k - 1} \frac{X_k}{A_k}\end{aligned}

   and it can be seen that the specific gas constant, :math:`R \equiv c_p - c_v` is
   independent of the :math:`\gamma_i`, and is simply :math:`R = k/m_u\bar{A}` giving the
   usual relation that :math:`p = R\rho T`. Furthermore, we can show

   .. math::

      \Gamma_1 \equiv \left . \frac{\partial \log p}{\partial \log \rho} \right |_s =
         \left ( \sum_k \frac{\gamma_k}{\gamma_k - 1} \frac{X_k}{A_k} \right ) \bigg /
         \left ( \sum_k \frac{1}{\gamma_k - 1} \frac{X_k}{A_k} \right ) =
      \frac{c_p}{c_v} \equiv \gamma_\mathrm{effective}

   and :math:`p = \rho e (\gamma_\mathrm{effective} - 1)`.

   This equation of state takes several runtime parameters that can set the
   :math:`\gamma_i` for a specific species:

   -  eos_gamma_default: the default :math:`\gamma` to apply for
      all species

   -  species_X_name and species_X_gamma: set the :math:`\gamma_i`
      for the species whose name is given as species_X_name to the
      value provided by species_X_gamma. Here, X can be one
      of the letters: a, b, or c, allowing us to specify
      custom :math:`\gamma_i` for up to three different species.

The thermodynamic quantities are stored in a Fortran type eos_t,
which has fields for all the thermodynamic inputs and outputs. The
type definition is brought in through eos_type_module.
 [4]_

The first argument to the eos call is an integer key that
specifies which thermodynamic variables (in addition to the mass
fractions) are used as input. EOS input options are listed
in table \ `[arch:table:eosinput] <#arch:table:eosinput>`__.

.. table:: [arch:table:eosinput] EOS input flags

   +--------------+-------------------------+
   | key          | input quantities        |
   +==============+=========================+
   | eos_input_rt | :math:`\rho`, :math:`T` |
   +--------------+-------------------------+
   | eos_input_rh | :math:`\rho`, :math:`h` |
   +--------------+-------------------------+
   | eos_input_tp | :math:`T`, :math:`p`    |
   +--------------+-------------------------+
   | eos_input_rp | :math:`\rho`, :math:`p` |
   +--------------+-------------------------+
   | eos_input_re | :math:`\rho`, :math:`e` |
   +--------------+-------------------------+
   | eos_input_ps | :math:`p`, :math:`s`    |
   +--------------+-------------------------+

fill_3d_module
--------------

The fill_3d_module provides routines that map from the 1D
base state to the full Cartesian 2- or 3D state. Variations in the
routines allow for cell-centered or edge-centered data on either the
base state or full Cartesian state.

fundamental_constants_module
----------------------------

The fundamental_constants_module provides a simple list of
various fundamental constants (e.g. Newton’s gravitational constant)
in CGS units.

geometry
--------

network
-------

The network module defines the number species advected by the
code (nspec), their ordering, and gives their basic properties
(like atomic number, :math:`Z`, and atomic mass, :math:`A`). All MAESTRO problems
require a network module, even if there are no reactions
modeled. Many different reaction modules (containing different sets
of isotopes) exist in Microphysics/networks. The particular network
used by a problem is defined in the problem’s GNUmakefile.

To find the location of a particular species (for instance, “carbon-12”)
in the allowed range of 1:nspec, you do the following query:

::

      ic12 = network_species_index("carbon-12")

If the resulting index is -1, then the requested species was not
found.

.. _sec:probin:

probin_module
-------------

probin_module provides access to the runtime parameters.
The runtime parameters appear simply as module variables. To get the
value of a parameter, one simply needs to ‘use probin_module’.
The preferred method is to add the ‘only’ clause to the
use statement and explicitly list only those parameters that
are used in the routine. Defining new runtime parameters is
described in § \ `3.2 <#sec:def_runtime_param>`__.

.. _sec:variables_module:

variables
---------

The variables module provides integer keys to index the state
multifabs and other arrays dealing with the scalar quantities. The
most commonly used keys are are list in table \ `[arch:table:variables] <#arch:table:variables>`__.

.. table:: [arch:table:variables] Common variables module keys

   +-----------+------------------------------------------------------------+
   | rho_comp  | density                                                    |
   +-----------+------------------------------------------------------------+
   | rhoh_comp | density :math:`\times` specific enthalpy, :math:`(\rho h)` |
   +-----------+------------------------------------------------------------+
   | spec_comp | first species partial density, :math:`(\rho X_1)`          |
   +-----------+------------------------------------------------------------+
   | temp_comp | temperature                                                |
   +-----------+------------------------------------------------------------+

The species indices are contiguous in the state array, spanning
spec_comp:spec_comp-1+nspec. To find a particular species, a
query can be made through the network module, such as:

::

      ic12 = network_species_index("carbon-12")

and then the fab can be indexed using spec_comp-1+ic12 to
get “carbon-12”.
The variables module also provides keys for the plotfile
variables and boundary condition types.

Other keys in the variables modules are reserved for boundary
conditions (foextrap_comp and hoextap_comp), the
projection of the pressure (press_comp), or constructing
the plotfile.

AMReX Helper Modules
====================

There are a large number of modules in amrex/ that provide
the core functionality for managing grids. Here we describe
the most popular such modules.

bl_types
--------

The main purpose of this module is to define the Fortran kind dp_t
which is used throughout the code to declare double precision variables.

bl_constants
------------

This module provides descriptive names for a number of common double precision
numbers, e.g. ONE = 1.0_dp_t. This enhances the readability of
the code.

parallel
--------

All MPI calls are wrapped by functions in the parallel module. For
serial jobs, the wrappers simply do the requested operation on processor.
By wrapping the calls, we can easily switch between serial and parallel
builds.

[sec:example] Example: Accessing State and MAC Data
===================================================

In MAESTRO, the state data is stored in a cell-centered multifab array
(the array index refers to the AMR level) and the MAC velocities are
stored in a 2D nodal multifab array (with indices referring to the AMR
level and the velocity component). Here we demonstrate a typical way
to extract the state and MAC velocity data.

All MAESTRO routines are contained in a module, to allow for compile-time
argument checking.

::

    module example_module

    contains

The main interface to our routine is called example—this will
take the multifabs containing the data and then pass them to the
work routines, example_2d or example_3d, depending on
the dimensionality.

::

      subroutine example(mla,s,umac,dx,dt)

        use bl_types
        use multifab_module
        use ml_layout_module
        use variables, only: rho_comp

Here, the bl_types and multifab_module
modules bring in the basic AMReX data types. Specifically, here,
bl_types defines dp_t which is the kind used for
declaring double precision data, and multifab_module defines
the multifab data type. The ml_layout_module defines the
datatype for a ml_layout—many routines will take an ml_layoutto
allow us to fill ghostcells. The variables module is a
MAESTRO module that provides integer keys for indexing the state
arrays. In this case the integer rho_comp refers to the
location in the state array corresponding to density.

Next we declare the subroutine arguments:

::

        type(ml_layout), intent(in   ) :: mla
        type(multifab) , intent(inout) :: s(:)
        type(multifab) , intent(inout) :: umac(:,:)
        real(kind=dp_t), intent(in   ) :: dx(:,:),dt

Here, s(:) is our multifab array that holds the state data.
with the array index in s refers to the AMR level. The MAC
velocities are held in the multifab umac, with the array
indices referring to the AMR level and the component.

Local variable declarations come next:

::

        ! Local variables
        real(kind=dp_t), pointer :: sp(:,:,:,:)
        real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
        integer :: i,n,dm,nlevs,ng_sp,ng_um
        integer :: lo(mla%dim),hi(mla%dim)

Amongst the local variables we define here are a pointer,
sp, that will point to a single fab from the
multifab s, and a pointer for each component of the MAC
velocity, ump, vmp, and wmp (for a 2D run,
we won’t use wmp). We note that regardless of the dimensionality,
these pointers are 4-dimensional: 3 spatial + 1 component.

Next we get the dimensionality and number of levels

::

        dm = mla%dim
        nlevs = mla%nlevel

Each multifab can have their own number of ghostcells, so we get
these next:

::

        ng_sp = nghost(s(1))
        ng_um = nghost(umac(1,1))

By convention, all levels in a given multifab have the same number of
ghostcells, so we use level 1 in the nghost() call. We also use
the same number of ghostcells for each component of the velocity, so
we only need to consider the first component in the nghost()
call. The ghostcells will be needed to access the data stored in the
fabs.

To access the data, we loop over all the levels, and all the boxes in
the given level.

::

        do n=1,nlevs
           do i = 1, nfabs(s(n))

nfabs(s(n)) is simply the number of boxes in level n on
the current processor. Each processor knows which fabs in its
multifabare local to that processor, and this loop will only loop
over those.

For a given box, we get the data and the bounds of the box.

::

              sp  => dataptr(s(n), i)
              ump => dataptr(umac(n,1),i)
              vmp => dataptr(umac(n,2),i)
              lo =  lwb(get_box(s(n), i))
              hi =  upb(get_box(s(n), i))

The actual data array is accessed through the dataptr function,
which takes a multifab (e.g. s(n)) and the index of the
box (i) we want. We see that the :math:`x` MAC velocity for the
current box is stored in ump and the :math:`y` MAC velocity is stored
in vmp. We don’t get the :math:`z` velocity data here, since that
would not be available for a 2D run—we defer that until we test on
the dimensionality below.

Finally, the index bounds of the box (just the data, not the ghostcells) are
stored in the dm-dimensional arrays lo and hi. These indices
refer to the current box, and hold for both the state, sp, and the MAC
velocity, ump and vmp. However, since the MAC velocity is nodal
in the component direction, the loops over the valid data will differ
slight (as we see below).

With the data extracted, we call a subroutine to operate on it. We use
different subroutines for the different dimensionalities (and many times
have a separate routine for spherical geometries).

::

              select case (dm)
              case (2)
                 call example_2d(sp(:,:,1,rho_comp),ng_sp, &
                                 ump(:,:,1,1),vmp(:,:,1,1),ng_um, &
                                 lo,hi,dx(n,:),dt)
              case (3)
                 wmp => dataptr(umac(n,3),i)
                 call example_3d(sp(:,:,:,rho_comp),ng_sp, &
                                 ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                 lo,hi,dx(n,:),dt)
              end select
           enddo    ! end loop over boxes

        enddo    ! end loop over levels

      end subroutine example

We call either the function
example_2d for two-dimensional data or example_3d
for three-dimensional data. Note that in the two-dimensional
case, we index the data as sp(:,:,1,rho_comp). Here a
‘1’ is used as the ‘z’-coordinate spatial index, since this
is a 2D problem, and the density component of the state is selected
(using the integer key rho_comp). The 3D version accesses
the data as sp(:,:,:,rho_comp)—only the component regarding
the variable is needed here. Notice that we also pass through
the number of ghostcells for each of the quantities.

This routine will be supplimented with example_2d and
example_3d, which actually operate on the data. The form of
the 2D function is:

::

      subroutine example_2d(density,ng_sp, &
                            umac,vmac,ng_um, &
                            lo,hi,dx,dt)

        use bl_constants_module
        use probin_module, only: prob_lo

        integer        , intent(in) :: lo(:),hi(:), ng_sp, ng_um
        real(kind=dp_t), intent(in) :: density(lo(1)-ng_sp:,lo(2)-ng_sp:)
        real(kind=dp_t), intent(in) ::    umac(lo(1)-ng_um:,lo(2)-ng_um:)
        real(kind=dp_t), intent(in) ::    vmac(lo(1)-ng_um:,lo(2)-ng_um:)

        real(kind=dp_t), intent(in) :: dx(:),dt

        integer         :: i, j
        real(kind=dp_t) :: x, y
        real(kind=dp_t) :: dens, u, v

        do j = lo(2), hi(2)
           y = prob_lo(2) + (dble(j) + HALF)*dx(2)

           do i = lo(1), hi(1)
              x = prob_lo(1) + (dble(i) + HALF)*dx(1)

              dens = density(i,j)

              ! compute cell-centered velocity
              u = HALF*(umac(i,j) + umac(i+1,j))
              v = HALF*(vmac(i,j) + vmac(i,j+1))

              ! operate on the data
              ! ...

           enddo
        enddo

      end subroutine example_2d

    end module example_module

In this function, the bounds of the density array take
into account the ng_sp ghostcells and the index space of the
current box. Likewise, the MAC velocities refer to the ng_um
ghostcells. The j and i loops loop over all the valid
zones. Coordinate information is computed from dx and
prob_lo which is the physical lower bound of the domain.
bl_constants_module declares useful double-precision
constants, like HALF (0.5). Here, we see how to access the
density for the current zone and compute the cell-centered velocities
from the MAC velocities. By convection, for a nodal array, the
indices refer to the *lower* interface in the nodal direction, so
for umac, umac(i,j) and umac(i+1,j) are the :math:`x` MAC
velocities on the lower and upper edge of the zone in the
:math:`x`-direction.

The three-dimensional case is similar, with the density array
declared as

::

      density(lo(1)-ng_sp:,lo(2)-ng_sp:,lo(3)-ng_sp:)

and an additional loop over the ‘z’ coordinate (from lo(3) to
hi(3)).

In this example, we looped over the valid zones. If we wished to loop
over the interfaces bounding the valid zones, in the :math:`x`-direction,
we would loop as

::

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            ! access umac(i,j)
         enddo
      enddo

Filling Ghostcells
==================

Ghostcells are filled through a variety of different routines, depending
on the objective.

-  multifab_fill_boundary fills ghost cells for two
   adjacent grids at the same level, which als includes periodic domain
   boundary ghost cells.

-  multifab_physbc fills ghostcells at the physical boundaries.

-  multifab_fill_ghost_cells is used for multilevel
   problems, and fills ghostcells in the finer grid (level n) by
   interpolating from data in the coarser grid (level n-1).
   This function, by default, will also call multifab_fill_boundary
   and multifab_physbc for both levels n and n-1 (you
   can override this behavior for speed optimization purposes).
   This call is usually preceded by a call to
   ml_cc_restriction_c which sets the level n-1 data to be
   the average of the level n data covering it.

You generally won’t see calls in the MAESTRO source code to these subroutines,
as there is now a special AMReX subroutine, ml_restrict_and_fill,
that takes an array of multifabs at different level, and in order calls:
(1) ml_cc_restriction_c, (2) multifab_fill_boundary,
(3) multifab_physbc, and (4) multifab_fill_ghost_cells.
These four subroutines are called in such a way to avoid extra
ghostcell filling, saving on communication time. You can specify the
starting component, starting boundary condition component,
the number of components, the number of ghost cells,
and whether or not you want to use the same boundary condition component
for all variables.

.. _sec:arch:bcs:

Boundary Conditions
===================

When MAESTRO is run, the boundary condition parameters are read in
from the input file and used to build the bc_tower. The
bc_tower consists of a bc_level object for each level of resolution
in the simulation. The bc_level contains 3 different descriptions of
the boundary conditions for each box in the domain at that level of
refinement: phys_bc_level_array, adv_bc_level_array,
and ell_bc_level_array. In all cases, the boundary
conditions are specified via integer values that are defined in
bc_module (part of AMReX).

Each level has a phys_bc_level_array(0:ngrids,dim,2) array,
where ngrids is the number of boxes on that level, dim is
the coordinate direction, and the last index refers to the lower (1)
or upper (2) edge of the box in the given coordinate direction. This
stores the *physical desciption* of the boundary type (outlet, inlet,
slipwall, etc.)—this description is independent of the variables
that live on the grid. The phys_bc_level_array(0,:,:) ‘box’
refers to the entire domain. If an edge of a box is not on a physical
boundary, then it is set to a default value (typically
INTERIOR). These boundary condition types are used to interpret
the actual method to fill the ghostcells for each variable, as
described in adv_bc_level_array and
ell_bc_level_array.

Whereas phys_bc_level_array provides a physical description
of the type of boundary, the array adv_bc_level_array
describes the *action* taken (e.g. reflect, extrapolate, etc.)
for each variable when filling boundaries.
adv_bc_level_array specifically describes the boundary
conditions that are in play for the advection (hyperbolic) equations.
The form of this array is
adv_bc_level_array(0:ngrids,dim,2,nvar) where the additional
component, nvar, allows for each state variable that lives on a
grid to have different boundary condition actions associated with it.
The convention is that the first dm variables in bc_level(where dm is
the dimensionality of the simulation) refer to the
velocity components, and the subsequent slots are for the other
standard variables described in the variables_module. For
instance, to reference the boundary condition for density, one would
index with dm+rho_comp. For temporary variables that are
created on the fly in the various routines in MAESTRO there may not
be a variable name in variables_module that describes the
temporary variable. In this case, the special variables
foextrap_comp and hoextrap_comp (first-order and high-order
extrapolation) are used.

ell_bc_level_array is the analog to
adv_bc_level_array for the elliptic solves in MAESTRO. This
will come into play in the multigrid portions of the code. The
actions that are used for ell_bc_level_array are either
Dirichlet or Neumann boundary condtions. For the velocity
projections, we are dealing with a pressure-like quantity, :math:`\phi`, so
the pressure boundary conditions here reflect the behavior we want for
the velocity. After the projection, it is :math:`\nabla \phi` that modifies
the velocity field. At a wall or for inflow conditions, we already
have the velocity we want at the boundary, so we want the velocity to
remain unchanged after the projection. This requires :math:`d\phi/dn=0` on
those boundaries. For outflow, we impose a condition that we do not
want the boundaries to introduce any tangental acceleration (or
shear), this is equivalent to setting :math:`\phi = 0` (then :math:`\partial
\phi/\partial t = 0`, with :math:`t` meaning ‘tangental’). This allows the
velocity to adjust as needed to the domain (see, for example,
:raw-latex:`\cite{almgrenBellSzymczak:1996}`).

The actual filling of the ghostcells according to the descriptions
contained in the bc_tower is carried out by the multifab_physbc routine. When you have an EXT_DIR
condition in multifab_physbc (an specified in the inputs file
as inlet), the advection solver (via the slope routine) and
linear solvers will then assume that the value in the ghost cells is
equal to the value that actually lies on the wall.

Multigrid
=========

MAESTRO uses the multigrid solver to enforce the divergence
constraint both on the half-time edge-centered advective velocities
(the “MAC projection”) and on the final cell-centered velocities
(the “HG projection”). For the MAC projection, since the velocity
data is edge-centered (the MAC grid), the projection is cell-centered.
For the HG projection, since the velocity data is cell-centered, the
projection is node-centered. The
multigrid solver performs a number of V-cycles until the residual
drops by 10-12 orders-of-magnitude. There are several options that
affect how the multigrid solver behaves, which we describe below.
More detail on the multigrid solvers is given in Chapter \ `[ch:mg] <#ch:mg>`__.

Multilevel and Refinement Criteria
==================================

.. _arch:sec:particles:

Particles
=========

MAESTRO has support for Lagrangian particles that are passively
advected with the velocity field. These are useful for diagnostics
and post-processing. To use particles, particles must be seeded into
the domain by writing a problem-specific init_particles.f90
routine. This routine is called at the start of the simulation. The
init_particles routines add particles at specific locations by
calling the particle_module’s add routine when a given
criteria is met by the fluid state.

When you run the code, particles are enabled by setting
use_particles = T. At the end of each timestep the locations of
all the particles are written out into a series of files called
timestamp_NN, where NN is the CPU number on which the
particle *currently* resides. Particles are always kept on the
processor containing the state data corresponding to their present
location. Several bits of associated data (density, temperature, and
mass fractions) are stored along with the particle ID and position.

Some simple python scripts allow for the plotting of the particle
positions. See § \ `[analysis:sec:particles] <#analysis:sec:particles>`__ for details.

Regression Testing
==================

There is an extensive regression test suite for AMReX that works with
MAESTRO. Full details, and a sample MAESTRO configuration file are
provided in the AMReX User’s Guide and source.

.. [1]
   Spherical geometry
   only exists for 3-d. This is a design decision—convection is 3-d.
   You can however run as an octant

.. [2]
   Note: many more compatible routines are available in the separate Microphysics git repo

.. [3]
   Microphysics is
   available at https://github.com/starkiller-astro/Microphysics. MAESTRO will
   find it via the MICROPHYSICS_HOME environment variable

.. [4]
   Note: an older interface to the EOS exists, but is
   deprecated. In this mode, the eos_old_interface module declares
   the variables that need appear in the old-style eos call
   argument list. MAESTRO routines use these module variables in the
   EOS call to avoid having to declare each quantity in each routine
   that calls the EOS. Most code has been updated to use the new interface.

.. |[fig:base_state] MAESTRO geometries, showing both the
1D base state and the full Cartesian state. (Left) For multi-level
problems in planar geometry, we force a direct alignment between the
radial array cell centers and the Cartesian grid cell centers by
allowing the radial base state spacing to change with space and
time. (Right) For multi-level problems in spherical geometry, since
there is no direct alignment between the radial array cell centers
and the Cartesian grid cell centers, we choose to fix the radial
base state spacing across levels. Figure taken
from :raw-latex:`\cite{multilevel}`.| image:: \archfigpath/base_grid
   :height: 2in
.. |[fig:base_state] MAESTRO geometries, showing both the
1D base state and the full Cartesian state. (Left) For multi-level
problems in planar geometry, we force a direct alignment between the
radial array cell centers and the Cartesian grid cell centers by
allowing the radial base state spacing to change with space and
time. (Right) For multi-level problems in spherical geometry, since
there is no direct alignment between the radial array cell centers
and the Cartesian grid cell centers, we choose to fix the radial
base state spacing across levels. Figure taken
from :raw-latex:`\cite{multilevel}`.| image:: \archfigpath/base_spherical
   :height: 2in
