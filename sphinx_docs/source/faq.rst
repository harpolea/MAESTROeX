**************************
Frequently Asked Questions
**************************

Coding
======

#. *Why is everything in its own module?*

   If a subroutine is in a Fortran module, then at compile time,
   there is argument checking. This ensures that the right number
   and data types of arguments are present. It makes the code safer.

#. *How do tags work when editing source?*

   Tags allow the editor to follow function/subroutine names and
   automatically switch you to the source code corresponding to that
   function. Using tags in MAESTRO depends on the editor:

   -  vi:

      In the build directory, type ‘make tags’. Then in
      vi, move the cursor over a function name and use
      ``^``-] to bring up the source corresponding to that
      function. Use ``^``-t to go back. (Here, ``^``-
      refers to the Control key.)

   -  emacs:

      In the build directory, type ‘make TAGS’. Then in
      emacs, move the cursor over a function name and use M-. to bring up the source corresponding to that function. Use
      M-\* to go back. (Here, M- refers to the META key.)

Compiling
=========

#. *What version of the Fortran standard is needed?*

   Some features of Fortran 2003 are used, in particular, the
   ISO_C_BINDING feature of Fortran 2003 is needed to define a long
   integer type for some MPI operations in Fortran. This is supported
   by most Fortran compilers even if they don’t support the entire
   Fortran 2003 standard.

   We also rely on the Fortran 95 standard that specifies that any
   local allocated arrays are automatically deallocated when a
   subroutine ends. Fortran 90 does not do this. Most
   MAESTRO routines rely on this Fortran 95 feature.

#. *The code doesn’t compile, but complains right away that there
   is “No rule to make target ‘fabio_c.c’, needed by ‘t/Linux.gfortran.mpi/c.depends’”*

   The environment variable AMREX_HOME needs to be the full path
   to the amrex/ directory. You cannot use ‘:math:`\sim`’ as a shortcut
   for your home directory.

#. *make issues an error like:*

   .. raw:: latex

      \small

   ::

              .../BoxLib/Tools/F_mk/GMakeMPI.mak:40: Extraneous text after `else' directive
              .../BoxLib/Tools/F_mk/GMakeMPI.mak:47: *** only one `else' per conditional.  Stop
             

   You need to use GNU make version 3.81 or later. Earlier versions do
   not support an *else-if* clause.

Running
=======

#. *How do we turn off all the initial projections to look at the
   initial velocity field as specified in initdata, instead of as
   modified by the velocity constraint?*

   ::

           max_step  = 1
           init_iter = 0
           init_divu_iter = 0
           do_initial_projection = F

   max_step
   init_iter
   init_divu_iter
   do_initial_projection

#. *MAESTRO crashes because the multigrid algorithm fails to
   converge—how do I get around this?*

   Setting general convergence criteria for multigrid is as much
   art as science.
   First, it is important to determine if the multigrd solver is
   close to convergence and just dancing around near the desired
   tolerance, but never reaching it, or if it is no where near
   convergence. For the latter, it may be that the multigrid
   solver was fed bad data and the problem arose in one of the earlier
   steps. To get more detail information from the multigrid solver,
   set mg_verbose to a positive integer from 1-4 (the higher
   the number the more information you receive.

   If the multigrid solver is failing during one of the initial
   “divu” iterations, it may be because the velocity is initially
   zero, so there is no velocity magnitude to use as a reference for
   convergence, and that (:math:`S - \bar{S}`) is very small (or zero). In
   this case, it is usually a good idea to perturb the initial state
   slightly, so the righthand side is non-zero.

   The tolerances used for the various multigrid solves in the code
   can be overridden on a problem-by-problem basis by putting a
   copy of MAESTRO/Source/mg_eps.f90 into the problem directory
   and resetting the tolerances. The role of each of these tolerances
   is described in MAESTRO/docs/mg/.

#. *Why do the initial projection and “divu” iters sometimes
   have a harder time converging than the multigrid solves in the main algorithm?*

   The initial projection and “divu” solve sets the density to :math:`1`
   (see § \ `[sec:flow:initialization] <#sec:flow:initialization>`__), so the coefficients in the
   elliptic solve are :math:`O(\beta_0) \sim O(\rho)`. But in the main
   algorithm, the coefficients are :math:`O(\beta_0/\rho) \sim O(1)`. Since
   :math:`\rho` can vary a lot, the variation in the coefficients in the
   initial projection and “divu” solve present a harded linear system
   to solve.

#. *How can I obtain profiling infomation for my run?*

   The code is already instrumented with timers. Simply compile with
   PROF=TRUE in the GNUmakefile, or equvalently do
   make PROF=t. A file containing a summary of the timings will
   be output in the run directory.

   An alternate way to get single-processor timings, when using the GCC
   compilers is to add -pg to the compilation flags for both
   gfortran and gcc. This can be accomplished by setting:

   ::

           GPROF := t

   in your GNUmakefile. Upon completion, a file
   names gmon.out will be produced. This can be processed by
   gprof running

   ::

           gprof exec-name

   where *exec-name* is the name of the executable. More detailed
   line-by-line timings can be obtained by using the -l argument
   to gprof.

#. *How can I force MAESTRO to abort?*

   In the output directory, do ‘touch .abort_maestro’. This
   will cause the code to write out a final checkpoint file, free up
   any allocated memory, and gracefully exit. Be sure to remove the
   .abort_maestro file before restarting the code in the
   same directory.

Debugging
=========

#. *How can we dump out a variable to a plotfile from any point in the
   code?*

   ::

           use fabio_module

           call fabio_ml_multifab_write_d(uold,mla%mba%rr(:,1),"a_uold")
           call fabio_ml_multifab_write_d(umac(:,1),mla%mba%rr(:,1),"a_umacx")

#. *How can I print out a multifab’s contents from within the code?*

   There is a print method in multifab_module. This can
   be simply called as

   ::

         call print(a)
         

   where a is a multifab (single-level).

#. *How can I debug a parallel (MPI) job with gdb?*

   If you only need to use a few processors, the following command will work:

   ::

       mpiexec -n 4 xterm -e gdb ./main.Linux.gfortran.mpi.exe 

   where the executable needs to be created with the “-g” flag to
   the compiler. This will pop up multiple xterms with gdb running
   in each. You need to then issue:

   ::

       run inputs

   where inputs is the desired inputs file *in each* xterm.

#. *How can I get more information about floating point exceptions?*

   AMReX can intercept floating point exceptions and provide a helpful
   backtrace that shows you where they were generated. See §
   `[ch:makefiles:special] <#ch:makefiles:special>`__.

I/O
===

#. *How can I tell from a plotfile what runtime parameters were
   used for its run? or when it was created?*

   In each plotfile directory, there is a file called job_info
   (e.g. plt00000/job_info) that lists the build directory and
   date, as well as the value of every runtime parameter for the run.

#. *How can I force the code to output a plotfile / checkpoint
   file at the next step?*

   In the output directory (where the code is running) do ‘touch
   .dump_plotfile’. This will create an empty file called
   .dump_plotfile. At the end of each step, if the code finds
   that file, it will output a plotfile. Simply delete the file to
   restore the code to its normal plotfile behavior.

   Similarly, creating the file .dump_checkpoint will force the
   output of a checkpoint file.

Algorithm
=========

#. *Why is MAESTRO so “hard” to use (e.g. as compared to a
   compressible code)?*

   There are several complexities to the algorithm that don’t have
   straightforward compressible counterparts. These mainly involve the
   role of the base state and the constraint equation.

   Care must be taken to setup an initial model/initial base state that
   respects the thermodynamics in MAESTRO and is in hydrostatic equilibrium.
   Best results are attained when the model is processed with the MAESTRO EOS and reset into HSE, as is done in the initial_model routines.
   Because MAESTRO builds off of the base state, any flaws in that initial
   state will influence the subsequent behavior of the algorithm.

   The constraint equation brings another complexity not seen in compressible
   codes—information is instantly communicated
   across the grid. In compressible codes you can track down a problem by
   watching where it starts from and watching it move one cell per dt. In
   MAESTRO things can go wrong in multiple places without it being obvious
   where the root problem is.

#. *In the final projection in the algorithm, we project
   :math:`U^{n+1}`, using a time-centered :math:`\beta_0`, a time-centered :math:`\rho_0`, but
   an “:math:`n+1`”-centered :math:`S`. Why then is the resulting :math:`\phi` (which then
   defines :math:`\pi`) is at “:math:`n+1/2`”?*

   The short answer to this question is that you should think of this
   as really projecting :math:`(U^{n+1} - U^n)` and the right hand side as having
   :math:`(S^{n+1} - S^n)`. This is because the pressure enters the dynamic equations as
   :math:`(U^{n+1} - U^n) = \ldots + \frac{1}{\rho^{n+1/2}} \nabla \pi^{n+1/2}`.
   (We approximate :math:`\pi^{n+1/2}` by :math:`\pi^{n-1/2}` then do the projection to fix the
   :math:`\pi` as well as the :math:`U`.)

   So everything is in fact time-centered.

#. *Why is :math:`\gammabar` computed as the average of the full state
   :math:`\Gamma_1` instead of computed from the base state density and
   pressure via the equation of state?*

   The primary reason is that there is no base state composition. The
   base state density is simply the average of the full state density,
   and the base state pressure is the pressure required for hydrostatic
   equilibrium. There is no thermodynamic relationship enforced between
   these base state quantities.

#. *Can I run a full star in 2-d axisymmetric geometry?*

   No. This is a design decision. There is no support for axisymmetric
   coordinates in MAESTRO. Spherical problems must be run in 3-d.

#. *Why did we switch all the equations over to the
   :math:`\tilde{\Ub}` form instead of just working with :math:`\Ub`?*

   This is basically a numerical discretization issue. Whenever the base
   state aligns with the grid, you should be able to show that you get
   exactly the same answer each way.

   When you do a spherical star on a 3d Cartesian grid, though, the :math:`w_0`
   is defined on the radial mesh and the :math:`\tilde{\Ub}` on the Cartesian
   mesh, and the :math:`w_0` part never experiences the Cartesian projection,
   for example. So there are differences in exactly how the :math:`w_0` component
   appears (projected on the Cartesian mesh vs. interpolated from the
   radial mesh)—we made the decision at the time to separate the
   components for that reason.

#. *Why does “checkerboarding” appear in the velocity field,
   especially in regions where the flow is stagnant?*

   Checkerboarding can arise from the projection—it doesn’t see that
   mode (because it is an approximate projection) so it is unable to
   remove it. This allows the pattern to slowly build up. There are
   filtering techniques that can be used to remove these modes, but
   they are not implemented in MAESTRO.

Analysis
========

#. *I want to open a plotfile, derive a new quantity from
   the data stored there, and write out a new plotfile with this derived
   data. How do I do this?*

   One implementation of this can be found in
   amrex/Tools/Postprocessing/F_Src/tutorial/fwrite2d.f90. This reads in
   the plotfile data using the plotfile_module that the
   data_processing routines rely on, but then builds a multifab
   and writes the data out to a plotfile using the AMReX write
   routines.
