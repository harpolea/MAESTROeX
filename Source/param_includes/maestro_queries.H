
// This file is automatically created by parse_maestro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

pp.query("do_reflux", do_reflux);
pp.query("maestro_verbose", maestro_verbose);
pp.query("probin_file", probin_file);
pp.query("model_file", model_file);
pp.query("perturb_model", perturb_model);
pp.query("print_init_hse_diag", print_init_hse_diag);
pp.query("stop_time", stop_time);
pp.query("max_step", max_step);
pp.query("cfl", cfl);
pp.query("init_shrink", init_shrink);
pp.query("small_dt", small_dt);
pp.query("max_dt_growth", max_dt_growth);
pp.query("max_dt", max_dt);
pp.query("fixed_dt", fixed_dt);
pp.query("nuclear_dt_fac", nuclear_dt_fac);
pp.query("use_soundspeed_firstdt", use_soundspeed_firstdt);
pp.query("use_divu_firstdt", use_divu_firstdt);
pp.query("spherical", spherical);
pp.query("octant", octant);
pp.query("do_2d_planar_octant", do_2d_planar_octant);
pp.query("regrid_int", regrid_int);
pp.query("amr_buf_width", amr_buf_width);
pp.query("drdxfac", drdxfac);
pp.query("minwidth", minwidth);
pp.query("min_eff", min_eff);
pp.query("use_tpert_in_tagging", use_tpert_in_tagging);
pp.query("plot_int", plot_int);
pp.query("plot_deltat", plot_deltat);
pp.query("chk_int", chk_int);
pp.query("plot_h_with_use_tfromp", plot_h_with_use_tfromp);
pp.query("plot_spec", plot_spec);
pp.query("plot_omegadot", plot_omegadot);
pp.query("plot_Hext", plot_Hext);
pp.query("plot_Hnuc", plot_Hnuc);
pp.query("plot_eta", plot_eta);
pp.query("plot_trac", plot_trac);
pp.query("plot_base", plot_base);
pp.query("plot_gpi", plot_gpi);
pp.query("plot_cs", plot_cs);
pp.query("plot_base_name", plot_base_name);
pp.query("check_base_name", check_base_name);
pp.query("diag_buf_size", diag_buf_size);
pp.query("plot_ad_excess", plot_ad_excess);
pp.query("plot_processors", plot_processors);
pp.query("plot_pidivu", plot_pidivu);
pp.query("init_iter", init_iter);
pp.query("init_divu_iter", init_divu_iter);
pp.query("restart_file", restart_file);
pp.query("restart_into_finer", restart_into_finer);
pp.query("do_initial_projection", do_initial_projection);
pp.query("mg_verbose", mg_verbose);
pp.query("cg_verbose", cg_verbose);
pp.query("mg_cycle_type", mg_cycle_type);
pp.query("hg_cycle_type", hg_cycle_type);
pp.query("hg_bottom_solver", hg_bottom_solver);
pp.query("mg_bottom_solver", mg_bottom_solver);
pp.query("max_mg_bottom_nlevels", max_mg_bottom_nlevels);
pp.query("mg_bottom_nu", mg_bottom_nu);
pp.query("mg_nu_1", mg_nu_1);
pp.query("mg_nu_2", mg_nu_2);
pp.query("hg_dense_stencil", hg_dense_stencil);
pp.query("do_sponge", do_sponge);
pp.query("sponge_kappa", sponge_kappa);
pp.query("sponge_center_density", sponge_center_density);
pp.query("sponge_start_factor", sponge_start_factor);
pp.query("plot_sponge_fdamp", plot_sponge_fdamp);
pp.query("anelastic_cutoff", anelastic_cutoff);
pp.query("base_cutoff_density", base_cutoff_density);
pp.query("burning_cutoff_density", burning_cutoff_density);
pp.query("buoyancy_cutoff_factor", buoyancy_cutoff_factor);
pp.query("dpdt_factor", dpdt_factor);
pp.query("do_planar_invsq_grav", do_planar_invsq_grav);
pp.query("planar_invsq_mass", planar_invsq_mass);
pp.query("evolve_base_state", evolve_base_state);
pp.query("use_exact_base_state", use_exact_base_state);
pp.query("fix_base_state", fix_base_state);
pp.query("do_smallscale", do_smallscale);
pp.query("do_eos_h_above_cutoff", do_eos_h_above_cutoff);
pp.query("enthalpy_pred_type", enthalpy_pred_type);
pp.query("species_pred_type", species_pred_type);
pp.query("use_delta_gamma1_term", use_delta_gamma1_term);
pp.query("use_etarho", use_etarho);
pp.query("slope_order", slope_order);
pp.query("grav_const", grav_const);
pp.query("ppm_type", ppm_type);
pp.query("bds_type", bds_type);
pp.query("ppm_trace_forces", ppm_trace_forces);
pp.query("beta0_type", beta0_type);
pp.query("use_linear_grav_in_beta0", use_linear_grav_in_beta0);
pp.query("rotational_frequency", rotational_frequency);
pp.query("co_latitude", co_latitude);
pp.query("rotation_radius", rotation_radius);
pp.query("mach_max_abort", mach_max_abort);
pp.query("drive_initial_convection", drive_initial_convection);
pp.query("stop_initial_convection", stop_initial_convection);
pp.query("restart_with_vel_field", restart_with_vel_field);
pp.query("use_alt_energy_fix", use_alt_energy_fix);
pp.query("use_thermal_diffusion", use_thermal_diffusion);
pp.query("temp_diffusion_formulation", temp_diffusion_formulation);
pp.query("thermal_diffusion_type", thermal_diffusion_type);
pp.query("limit_conductivity", limit_conductivity);
pp.query("do_burning", do_burning);
pp.query("burner_threshold_species", burner_threshold_species);
pp.query("burner_threshold_cutoff", burner_threshold_cutoff);
pp.query("do_subgrid_burning", do_subgrid_burning);
pp.query("reaction_sum_tol", reaction_sum_tol);
pp.query("small_temp", small_temp);
pp.query("small_dens", small_dens);
pp.query("use_tfromp", use_tfromp);
pp.query("use_eos_e_instead_of_h", use_eos_e_instead_of_h);
pp.query("use_pprime_in_tfromp", use_pprime_in_tfromp);
pp.query("s0_interp_type", s0_interp_type);
pp.query("w0_interp_type", w0_interp_type);
pp.query("s0mac_interp_type", s0mac_interp_type);
pp.query("w0mac_interp_type", w0mac_interp_type);
pp.query("use_particles", use_particles);
pp.query("store_particle_vels", store_particle_vels);
pp.query("do_heating", do_heating);
pp.query("sdc_iters", sdc_iters);
pp.query("sdc_couple_mac_velocity", sdc_couple_mac_velocity);
pp.query("eps_init_proj_cart", eps_init_proj_cart);
pp.query("eps_init_proj_sph", eps_init_proj_sph);
pp.query("eps_divu_cart", eps_divu_cart);
pp.query("eps_divu_sph", eps_divu_sph);
pp.query("divu_iter_factor", divu_iter_factor);
pp.query("divu_level_factor", divu_level_factor);
pp.query("eps_mac", eps_mac);
pp.query("eps_mac_max", eps_mac_max);
pp.query("mac_level_factor", mac_level_factor);
pp.query("eps_mac_bottom", eps_mac_bottom);
pp.query("eps_hg", eps_hg);
pp.query("eps_hg_max", eps_hg_max);
pp.query("hg_level_factor", hg_level_factor);
pp.query("eps_hg_bottom", eps_hg_bottom);
