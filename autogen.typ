STRUCT(CatalysisRule)
	VAR(u32, a, 0, Binary())
	VAR(u32, b, 0, Binary())
END

STRUCT(BindSite)
	VAR(pose, mol_from_sit, identity(), InputFloat())
	VAR(u8, flavor, 0, InputInt())
	VAR(u8, gender, 0, InputInt())
END

STRUCT(MoleculeType)
	BIN(BindSite, sites, 6, Struct())
	VAR(u16, bond_max_count, 0, InputInt())
	VAR(float, rad, 1.f, InputFloat())
	VAR(float, mass, 1.f, InputFloat())
	VAR(float, moi, ((2.f/5.f)*1.f*1.f*1.f), InputFloat()) // sphere moi is 2/5mr^2
	//VAR(float, cor, 1.f, InputFloat()) // coeficient of restitution
END

STRUCT(Molecule)
	VAR(s32, body_idx, -1, InputInt())	
	VAR(u8, type, 0, InputInt())
	VAR(u8, sites_taken, 0, Binary())
	VAR(pose, bod_from_mol, identity(), InputFloat())
	VAR(s32, l_bond, -1, InputInt())
	VAR(s32, r_bond, -1, InputInt())
END

STRUCT(Body)
	VAR(pose, wol_from_bod, identity(), InputFloat())
	VAR(vec2, bod_from_com, vec2(0.f), InputFloat())
	VAR(vec2, dp, vec2(0.f), InputFloat())
	VAR(float, dr, 0.f, InputFloat())	
	VAR(float, mass, 0.f, InputFloat())
	VAR(float, moi, 0.f, InputFloat())
	VAR(float, vibration, 0.f, InputFloat())
	VAR(s32, some_mol_idx, -1, InputInt())
	VAR(u32, types_inside, 0, Binary())
	VAR(u32, mol_count, 0, InputInt())
	VAR(u8, type, (u8)-1, InputInt())
	VAR(bool, valid, true, Checkbox())
END

STRUCT(BodyTemplate)
	BUNCH(Molecule, mols, Struct())
END

STRUCT(BodyHist)
	VAR(s32, body_idx, -1, InputInt())
	BIN(float, body_type_masses, BODY_TYPE_COUNT_MAX, Struct())
END

STRUCT(HeatPumpParams)
	VAR(bool, enable, false, Checkbox())
	VAR(bool, cycle_enable, false, Checkbox())
	VAR(float, avg_energy_min, 100.f, InputFloat())
	VAR(float, avg_energy_max, 100.f, InputFloat())
	VAR(float, permeability, 0.01f, InputFloat())
	VAR(float, frequency, 0.5f, InputFloat())
END

STRUCT(SimStats)
	VAR(u64, step_count, 0, InputInt())
	VAR(float, mass_sum, 0.f, InputFloat())
	VAR(vec2, lin_momentum_sum, vec2(0.f), InputFloat())
	VAR(float, lin_kinetic_sum, 0.f, InputFloat())
	VAR(float, rot_kinetic_sum, 0.f, InputFloat())
	VAR(float, rot_momentum_sum, 0.f, InputFloat())
	VAR(float, vibration_sum, 0.f, InputFloat())
	VAR(s32, bod_count, 0, InputInt())
	VAR(s32, mol_count, 0, InputInt())
	BIN(float, speed_hist, 72, Struct())
	VAR(float, bucket_count_max, 0.f, InputFloat())
	VAR(float, speed_range, 250.f, InputFloat())
	BIN(s32, compound_counts, BODY_TYPE_COUNT_MAX, Struct())
	BUNCH(BodyHist, shell_guts, Struct())
	VAR(s32, free_capsid_count, 0, InputInt())
	VAR(s32, bound_capid_count, 0, InputInt())
	VAR(s32, shelled_capsid_count, 0, InputInt())
END


STRUCT(SimState)
	VAR(float, dt, (1.f/60.f), InputFloat())
	VAR(float, time, 0.f, InputFloat())
	VAR(u64, step_count, 0, InputInt())
	VAR(float, bound_rad, 200.f, InputFloat())
	VAR(HeatPumpParams, heat_pump, HeatPumpParams(), Struct())
	BUNCH(Molecule, molecules, Struct())
	BUNCH(Body, bodies, Struct())
	BUNCH(s32, free_bodies, Struct())
	VAR(u32, coll_hash_salt, 0x1337c0de, InputInt())
	VAR(u32, rng_state, 0xdeadbeef, InputInt())
	VAR(float, capsid_break_threshold, 100.f, InputFloat())
	VAR(float, dissapated_lin_kinetic_energy, 0.f, InputFloat())
	VAR(float, dissapated_rot_kinetic_energy, 0.f, InputFloat())
	VAR(vec2, lin_bound_momentum, vec2(0.f), InputFloat())
	VAR(float, rot_bound_momentum, 0.f, InputFloat())
	VAR(SimStats, init_stats, SimStats(), Struct())
	// runtime compound description, SoA for speed?
	BIN(MoleculeType, mol_types, BODY_TYPE_COUNT_MAX, Struct()) 
	BIN(CatalysisRule, catalysis_rules, BODY_TYPE_COUNT_MAX, Struct())
	BIN(vec3, body_cols, BODY_TYPE_COUNT_MAX, Struct())
	BIN(u32, mol_types_inside, BODY_TYPE_COUNT_MAX, Struct())
	BIN(float, catalysis_chances, BODY_TYPE_COUNT_MAX, InputFloat())
END

STRUCT(SimVis)
	VAR(bool, show_cells, false, Checkbox())
	VAR(bool, show_collisions, false, Checkbox())
	VAR(bool, show_velocity, false, Checkbox())
	VAR(bool, color_body_idx, false, Checkbox())
	VAR(bool, show_sites, false, Checkbox())
	VAR(bool, show_bonds, false, Checkbox())
	VAR(bool, show_com, false, Checkbox())
	VAR(bool, show_mol_rotation, false, Checkbox())
	VAR(bool, show_rings, false, Checkbox())
	VAR(bool, color_mol_idx, false, Checkbox())
	BIN(bool, filter_type, BODY_TYPE_COUNT_MAX, Struct())
	VAR(pose, wol_from_cam, identity(), InputFloat())
	VAR(vec3, bg_col, vec3(0.f), Color()) 
	VAR(vec3, fg_col, vec3(1.f), Color()) 
	BIN(vec3, type_cols, BODY_TYPE_COUNT_MAX, Color())
	VAR(bool, show_center, false, Checkbox())
	VAR(bool, show_shell_type, false, Checkbox())
END

// user compound description, AoS for convenience
STRUCT(CompoundDesc)
	VAR(s32, comp_class, 0, InputInt())
	VAR(vec3, col, vec3(1.f), Color())
	// mass, rad, etc

	// catalyst specific
	VAR(u8, catalyzes_compound, NO_TYPE, InputInt())
	VAR(float, catalysis_chance, 1.f, SliderFloat())

	// substrate specific
	BIN(u8, mol_types, 6, Struct())

	// capsid specific
	VAR(u16, bond_max_count, 24, InputInt()) // move to mol desc?
END

STRUCT(EnvParams)
	VAR(float, bound_rad, 100.f, InputFloat())
	VAR(HeatPumpParams, heat_pump, HeatPumpParams(), Struct())
END

STRUCT(TypeParams)
	BIN(CompoundDesc, comps, BODY_TYPE_COUNT_MAX, Struct())
	VAR(float, capsid_break_threshold, 200.f, InputFloat())
END

STRUCT(AutogenInitParams)
	BIN(int, compound_count, BODY_TYPE_COUNT_MAX, Struct())
END
STRUCT(InitGenParams)	
	VAR(float, spacing, 2.f, InputFloat())
	VAR(float, heat, 1.f, InputFloat())
	BIN(float, compound_distro, BODY_TYPE_COUNT_MAX, Struct())
	VAR(bool, one_starter, false, Checkbox())
	VAR(s32, starter_comp, 0, InputInt())
	VAR(bool, enable_seed_autogen, true, Checkbox())
	VAR(AutogenInitParams, seed_autogen, AutogenInitParams(), Struct())
END

STRUCT(SimParams)
	VAR(EnvParams, env, EnvParams(), Struct())
	VAR(TypeParams, type, TypeParams(), Struct())
	VAR(InitGenParams, init, InitGenParams(), Struct())
END

STRUCT(Sim)
	VAR(SimParams, sim_params, SimParams(), Struct())
	VAR(SimState, sim_state, SimState(), Struct())
	VAR(SimVis, sim_vis, SimVis(), Struct())
END

STRUCT(SimCtrl)
	VAR(bool, run, false, Checkbox())
	VAR(bool, step, false, Checkbox())
	VAR(s32, substeps, 4, InputInt())
	VAR(s32, supersteps, 1, InputInt())
	VAR(bool, record_data, false, Checkbox())
	VAR(float, time_of_last_datum, 0.f, InputFloat())
	VAR(float, time_of_last_state, 0.f, InputFloat())
	BIN(s8, record_name, 64, Struct())
END