#include "nest.h"
#include "wrap/draw_wrap.h"
#include "core/log.h"
#include "core/vec2.h"
#include "core/pose.h"
#include "core/container.h"
#include "core/draw_stack.h"
#include "core/draw_shapes.h"
#include "core/draw_prim.h"
#include "core/random.h"
#include "core/cpu_timer.h"
#include "core/dirent.h"
#include "core/color.h"
#include "core/dir.h"

/*

bugs:
- !!! molecules can "tunnel" into shells, due to projection
- crash on edit and continue, inside TimerTree. very consistent, should be trackable?
- bounds collision check is using stale wol_from_mol cache
- if you "run" while in the editor there are overlapping molecules, there's a crash because sim_vis isn't initted properly. need to genfinalize on run? nasty
- rotational momentum is constantly growing?? is it every projection? seems like even when just moving...linearized motion?? what does it mean that the linear momentum of a system grows without bound? *something* is running away right? watch if dr is drifting? or if zw is?
- breaking apart a ring causes energy loss/gain. odd, need to double check all breakapart capsin logic. something with centers of mass? doesn't happen on catalysts, so maybe an artifact of 1:1 mass ratio of broken pieces?

wrong:
- bounds collision is using a stale position cache, big speedup but technically  wrong

warning:
- should really investigate backtracking or something, in the event of total interpenetration...(which can cause 0 length direction normal)
- collision pair hash gives many false positives (on collision of hash). it's salted every frame though so hopefully the collisions change each frame, and thus resolve themselves.
- wol_from_mol returns wol_from_body when mass is 1.0, because it assumes that also means identity transform

todo:
- change time to s64 all over the place, where seconds are needed replace with double simTime() calls that convert internal counter to seconds
- rename CompoundDesc to Type? it appears like that's the meaning everywhere, but I've forgotten if they are supposed to be subtly different?
- change from Heat to Avg Energy in InitGenParams, and update code to reflect this (be consistent with heat pump)
- use simulation step + dt to infer time, since it's a fixed timestep
- save heat pump state and logic into the sim, so that it's run from within the sim rather than via "control"
  - init heat pump cycle time such that the current temperature setting matches the initial temperature (no discontinuity)
  - disalow control over heatpump other than initial conditions (or rather don't make that the preferred UI)

opt:
- maybe i can collapse bod_from_com such that wol_from_com is always the body pose, and then we just have mol_from_com for molecule positions. not sure why i didn't do it this way? it would reduce *bod_from_com multiplies in various places.. and simplify integration...
- why do 2-piece bodies cost so much more to collide? i guess because you're constantly taxing the hash from being next to each other.
- basically re-computing wol_from_mol costs 1ms per frame each time we do it over all things. so caching would def help there i think, especially when it's simple

ideas:
- it's not clear if it's better to store the rotation as a complex angle or as a number. need to actually finish the work, and then look at all the cases to see how often we convert back and forth. for example the dr variable seems to make sense to be a real rather than complex, because you need to timestep scale it (which involves trig if we use a complex angle, negating the usefulness). also for rotational momentum it seems we'd want a single number for omega as well?
  but for poses, it seems valuable to keep it complex, since pose composition and pose*vec are reasonably common operations in the code (and have to be done cpu side)
- if a body is just a collection of circles or bodies, maybe we need to look at body flex and other such things for actual binding sites to "line up"
  in that case, it's a matter of aggregating the forces onto each circle, applying them, and then iterating verlet constraints on the linkages
  verlet linkages on spheres seems like a very good model though, for atoms and shapes.
  it also seems like it might be easier than compounding or parenting shapes on the fly...
- save the state of the editor on each "run" 
- add R to rewind back to whatever state was hit during last run
- for something like the "petri dish" case, where lots of mols are jammed towards the edge, i can make "jammed" or "penetrated with bounds" (or any other static obj) have artificially more mass, when it comes to projecting out from them
  in other words, molecules that are in contact with the edge become essentially "the edge" and perhaps that flag can propagate from any molecule that touches an edge molecule. so they temporarily "get to share in the infinity of mass" for the purposes or projection and bouncing.
- run sim forward, until a ring forms, containing stuff. mark which molecules make that ring up. rewind, and highlight those molecules to watch how they get to that situation. look for patterns. maybe visualize their paths.
- can try to make capsid molecules longer arcs, so fewer are needed to assemble a ring, but the ring can still be large enough to capture catalyst (substrate would still be at least the size of a capsid arc, so no help there)
- hover over types to make them pulse color - good for spotting lone catalyst
*/

#include "autogen.h"
#include "physics.h"

#include "core/type_magic.h"

#define TYPE() TOOL
#define VALUES(X) \
X(add) \
X(move) \
X(velocity) \
X(COUNT)
DECLARE_ENUM(TYPE(), VALUES)
#undef VALUES
#undef TYPE


enum {CLASS_INERT, CLASS_CAPSID, CLASS_CATALYST, CLASS_COMPOUND};
enum {SET_NONE, SET_VEL, SET_ROT};
enum {TOOL_ADD, TOOL_INSPECT};

struct SimEditor {
	pose cursor_from_body = identity();
	//int mol_type = 0;
	//int mol_idx = NO_IDX; // later, for making compounds
	int compound_type = 0;
	int bond_count = 0;
	int bod_idx = NO_IDX;
	int inspect_mol_idx = NO_IDX;
	int inspect_bod_idx = NO_IDX;
	int want_to = TOOL_INSPECT;
	int setting = SET_NONE;
	pose world_from_cursor_press = identity();
};
static SimEditor edit;

struct PartVis {
	bool is_colliding;
};

// generation
static SimParams sim_params;
static Bunch<BodyTemplate> gen_templates;

// state
static SimState sim_state;
static const float avg_energy_min = 0.f;
static const float avg_energy_max = 2500.f;

// control 
static SimCtrl sim_ctrl;

// visualization
static SimVis sim_vis;
static Bunch<PartVis> mols_vis;
static pose wol_from_cam = vec4(0,0,1,0);

// timing
static TimerTree timer;
static TimerResult timer_results[TIMER_TREE_MAX_TIMERS];

// io
static char sim_filename_out[64];
vec2 debug_c = vec2(0.f);
float debug_r = 5.f;

static void genBodyTemplates(Bunch<BodyTemplate>* templates, const TypeParams* type, const MoleculeType* mol_types) {
	templates->clear();
	for (int i = 0; i < type->comps.count; ++i) {
		const CompoundDesc& D = type->comps[i];
		BodyTemplate& T = templates->push();
		if (D.comp_class == CLASS_COMPOUND) {
			assert(D.mol_types.count == 2); //#DODGY for now just support 2, but can change this later
			if (D.mol_types[0] != NO_TYPE && D.mol_types[1] != NO_TYPE) {			
				Molecule A;
				A.type = D.mol_types[0];
				A.bod_from_mol = vec4(vec2(-mol_types[A.type].rad,0.f), 1.f, 0.f);
				for (int s = 0; s < mol_types[A.type].sites.count; ++s) fillSite(A, s, -1); //#DODGY bonds with no left/right..cause they aren't in the world yet so no indecies known...
				Molecule B;
				B.type = D.mol_types[1];
				B.bod_from_mol = vec4(vec2(+mol_types[B.type].rad,0.f), 1.f, 0.f);
				for (int s = 0; s < mol_types[B.type].sites.count; ++s) fillSite(B, s, -1);
				T.mols.push(A);
				T.mols.push(B);
			}
		} else {
			Molecule& M = T.mols.push();
			M.type = (u8)i;
		}
	}
}

static void simClear(SimState* S, const EnvParams* env, const TypeParams* type) {
	S->rng_state = 0xdeadbeef;
	S->coll_hash_salt = 0x1337c0de;
	S->step_count = 0;
	S->time = 0.f;

	S->molecules.clear();
	S->bodies.clear();
	S->free_bodies.clear();
	S->mol_types.clear();
	S->catalysis_rules.clear();
	S->body_cols.clear();
	S->mol_types_inside.clear();
	S->catalysis_chances.clear();

	S->bound_rad = env->bound_rad;
	S->heat_pump = env->heat_pump;
	S->capsid_break_threshold = type->capsid_break_threshold;


	S->mol_types.clear();
	for (int i = 0; i < type->comps.count; ++i) {
		const CompoundDesc& desc = type->comps[i];
		if (desc.comp_class == CLASS_CAPSID) {
			MoleculeType& capsid = S->mol_types.push();
			capsid.bond_max_count = desc.bond_max_count;
		
			float bind_ang = TAU/desc.bond_max_count;
			{
				BindSite site;
				site.mol_from_sit = trans(vec2(2.f*capsid.rad,0.f),bind_ang,1.f);
				site.mol_from_sit.zw(normalize(site.mol_from_sit.zw()));
				site.flavor = 0;
				site.gender = GENDER_MALE;
				capsid.sites.push(site);
			}
			{
				BindSite site;
				//site.mol_from_sit = trans(rotation(-bind_ang)*vec2(-2.f*MOL_RAD,0.f),0.f,1.f);
				site.mol_from_sit = rotation(-bind_ang) * trans(vec2(-2.f*capsid.rad,0.f), 0.f, 1.f);
				site.mol_from_sit.zw(normalize(site.mol_from_sit.zw()));
				site.flavor = 0;
				site.gender = GENDER_FEMALE;
				capsid.sites.push(site);
			}
		} else {
			// for catalysts and inert, push a type
			// for compounds, this is just a "dummy type" to make the indexing work out
			/* MoleculeType& M = */ S->mol_types.push();
		}
	}

	S->catalysis_rules.clear();
	S->mol_types_inside.clear();
	S->body_cols.clear();
	S->catalysis_chances.clear();

	for (int i = 0; i < type->comps.count; ++i) {
		const CompoundDesc& D = type->comps[i];
		u32 types_inside = 0;
		if (D.comp_class == CLASS_COMPOUND) {
			assert(D.mol_types.count == 2); //#DODGY for now just support 2, but can change this later
			if (D.mol_types[0] != NO_TYPE && D.mol_types[1] != NO_TYPE) {			
				types_inside |= (1<<D.mol_types[0]) | (1<<D.mol_types[1]);
			} else {
				assert(false);
			}
		} else {
			if (D.comp_class == CLASS_CATALYST && D.catalyzes_compound != NO_TYPE) {
				assert(D.catalyzes_compound < type->comps.count); // input here shouldn't be malformed
				CatalysisRule r;
				r.a = 1<<u8(i);
				const CompoundDesc& B = type->comps[D.catalyzes_compound];
				r.b = 0;
				for (int c = 0; c < B.mol_types.count; ++c) {
					r.b |= 1<<B.mol_types[c];
				}
				S->catalysis_rules.push(r);
			}
			types_inside |= 1<<u8(i);
		}
		S->mol_types_inside.push(types_inside);
		S->body_cols.push(D.col);
		S->catalysis_chances.push(D.catalysis_chance);
	}
}

void computeStats(SimStats& stats) {
	stats.step_count = sim_state.step_count;
	stats.compound_counts.pushi(0, sim_state.body_cols.count);
	stats.speed_hist.pushi(0.f, stats.speed_hist.maxcount);
	const Bunch<Body>& bodies = sim_state.bodies;

	u8 capsid_types_inside = 0;
	s16 bond_max_count = 0;
	for (int i = 0; i < sim_state.mol_types.count; ++i) {
		if (sim_state.mol_types[i].bond_max_count > 0) {
			capsid_types_inside = 1<<i;
			bond_max_count = sim_state.mol_types[i].bond_max_count;
			break;
		}
	}
	for (int i = 0; i < bodies.count; ++i) {
		if (!bodies[i].valid) continue;
		const Body& B = bodies[i];
		stats.mass_sum += B.mass;
		stats.lin_momentum_sum += B.dp * B.mass;
		float vsq = dot(B.dp,B.dp);
		stats.lin_kinetic_sum += 0.5f * B.mass * vsq; 
		stats.rot_kinetic_sum += 0.5f * B.moi * B.dr * B.dr;
		stats.vibration_sum += B.vibration;
		vec2 wol_from_com = B.wol_from_bod * B.bod_from_com;
		stats.rot_momentum_sum += det(wol_from_com, B.dp*B.mass) + B.moi*B.dr;
		if (isnan(stats.lin_kinetic_sum) || isnan(stats.rot_kinetic_sum)) {
//			int break_here = 0;
		}
		//assert(B.vibration >= 0.f);
		//assert(!isnan(stats.lin_kinetic_sum) && !isnan(stats.rot_kinetic_sum));
		float v = sqrt(vsq);
		int bucket = min(stats.speed_hist.maxcount-1, int(v/stats.speed_range*stats.speed_hist.maxcount));
		stats.speed_hist[bucket] += 1.f;
		stats.bucket_count_max = max(stats.speed_hist[bucket], stats.bucket_count_max);
		stats.bod_count += 1;
		stats.compound_counts[B.type] += B.mol_count;
		if (B.types_inside == capsid_types_inside) {
			if (B.mol_count == 1) {
				stats.free_capsid_count += B.mol_count;
			} else if (B.mol_count == bond_max_count) {
				stats.shelled_capsid_count += B.mol_count;
			} else {
				stats.bound_capid_count += B.mol_count;
			}
		}
	}
	stats.mol_count = (int)sim_state.molecules.count;
}

static void genInitialConditions(SimState* S, const InitGenParams params, const Bunch<BodyTemplate>* templates) {
	float exclude_radius = 0.0f;
	if (params.enable_seed_autogen) {
		// figure out maximum bond count for whatever ring configuration we have
		u16 bond_count = 0;
		float mol_rad = 0.0f;
		u8 type = NO_TYPE;

		for (int i = 0; i < templates->count; ++i) {
			for (int m = 0; m < (*templates)[i].mols.count; ++m) {
				u8 t = (*templates)[i].mols[m].type;
				int c = S->mol_types[t].bond_max_count;
				if (c > bond_count) {
					bond_count = c;
					mol_rad = S->mol_types[t].rad;
					type = t;
				}
			}
		}
		float t = TAU/bond_count;
		float L = mol_rad * 2.0f;
		float D = (L*0.5f) / tan(t*0.5f);
		float shell_radius = sqrt(sqr(L*0.5f) + sqr(D));
		exclude_radius = shell_radius + mol_rad;

		BodyTemplate T;
		for (int i = 0; i < bond_count; ++i) {
			float a = i * t;
			Molecule A;
			A.type = type;
			vec2 r = quat2d(a);
			float external_angle = TAU/4.0f - (TAU/2.0f - TAU/4.0f - t * 0.5f);
			A.bod_from_mol = vec4(r*shell_radius, rotate(r, quat2d(TAU/4.0f + external_angle)));
			T.mols.push(A);
		}
		shellAdd(S, T, identity());

		int idx = 1;
		for (int c = 0; c < params.seed_autogen.compound_count.count; ++c) {
			for (int i = 0; i < params.seed_autogen.compound_count[c]; ++i) {
				float r = 1.25f * sqrtf(float(idx)) * mol_rad;
				float t = idx*GOLDEN_ANGLE;
				vec2 p = quat2d(t)*r;
				pose wol_from_bod = trans(p, randFloat(0.f,TAU), 1.f);
				int bod_idx = compoundAdd(S, (*templates)[c], wol_from_bod);
				Body& B = bodyGet(&sim_state, bod_idx);
				B.dp = normalize(vec2(randFloat(-0.5f,+0.5f), randFloat(-0.5f,+0.5f))) * randFloat(0.f,params.heat);
				B.dr = randFloat(-params.heat,params.heat);
				B.type = c;
				idx += 1;
			}
		}
	}

	int b = (int)S->bound_rad;
	for (int x = -b; x <= b; ++x) {
		for (int y = -b; y <= b; ++y) {
			vec2 p = vec2((float)x,(float)y)*params.spacing;
			//vec2 p = (vec2(x,y) - vec2(b)*0.5f) * 2.f * 2.f;// * 2.f;
			float len = length(p);
			if (len < (S->bound_rad - MOL_RAD*2.f) && len > (exclude_radius + MOL_RAD*2.f)) {
				//pose wol_from_bod = trans(p, 0.f, 1.f); //trans(vec2(x,y)*2.f*2.f + vec2(randFloat(-0.5f,+0.5f), randFloat(-0.5f,+0.5f)), 0.f, 1.f);
				pose wol_from_bod = trans(p, randFloat(0.f,TAU), 1.f);
				u8 template_idx = (u8)roulette(params.compound_distro.ptr, params.compound_distro.count);
				int bod_idx = compoundAdd(S, (*templates)[template_idx], wol_from_bod);
				Body& B = bodyGet(&sim_state, bod_idx);
				B.dp = normalize(vec2(randFloat(-0.5f,+0.5f), randFloat(-0.5f,+0.5f))) * randFloat(0.f,params.heat);
				B.dr = randFloat(-params.heat,params.heat);
				B.type = template_idx;
			}
		}
	}

	if (params.one_starter) { 
		pose wol_from_bod = trans(vec2(0.f), randFloat(0.f, TAU), 1.f);
		//pose wol_from_bod = trans(vec2(-params.bound_rad*0.5, 0.f), randFloat(0.f, TAU), 1.f);
		u8 template_idx = (u8)params.starter_comp;//1;
		int bod_idx = compoundAdd(S, (*templates)[template_idx], wol_from_bod);
		Body& B = bodyGet(&sim_state, bod_idx);
		B.dp = normalize(vec2(randFloat(-0.5f, +0.5f), randFloat(-0.5f, +0.5f))) * randFloat(0.f, params.heat);
		B.dr = randFloat(-params.heat, params.heat);
		B.type = template_idx;
	}
}

CatalysisRule catalysisRule(u32 a, u32 b) {
	CatalysisRule r;
	r.a = a;
	r.b = b;
	return r;
}
static void genFinalize() {
	simFinalize(&sim_state);
	mols_vis.clear();
	mols_vis.pushi(sim_state.molecules.count);
	sim_ctrl.time_of_last_datum = 0.f;
	sim_ctrl.time_of_last_state = 0.f;
}

//#WARNING for editor use only, to keep re-setted states in sync with vibration
static void genZeroVibration() {
	Bunch<Body>& bodies = sim_state.bodies;
	for (int i = 0; i < bodies.count; ++i) {
		if (!bodies[i].valid) continue;
		bodies[i].vibration = 0.f;
	}
}
static void genStat() {
	sim_state.init_stats = SimStats();
	computeStats(sim_state.init_stats);
	sim_state.dissapated_lin_kinetic_energy = 0.f;
	sim_state.dissapated_rot_kinetic_energy = 0.f;
	sim_state.lin_bound_momentum = vec2(0.f);
	sim_state.rot_bound_momentum = 0.f;
}

static void genTypesAndTemplatesAndDistro(SimParams* params) {
	genBodyTemplates(&gen_templates, &params->type, sim_state.mol_types.ptr); // #DODGY referring to currently running sim to make templates? are templates necessarily locked to this? is this guaranteed? (i guess so since sim clear sets up the mol_types?)
	if (params->init.compound_distro.count != gen_templates.count) {
		params->init.compound_distro.clear();
		params->init.compound_distro.pushi(1.f/gen_templates.count, (int)gen_templates.count);
	}
	if (params->init.seed_autogen.compound_count.count != gen_templates.count) {
		params->init.seed_autogen.compound_count.clear();
		params->init.seed_autogen.compound_count.pushi(0, (int)gen_templates.count);
	}
	sim_vis.filter_type.pushi(true, max(0, (int)params->type.comps.count-sim_vis.filter_type.count));
	sim_vis.type_cols.pushi(vec3(1,0,1),  max(0, (int)params->type.comps.count-sim_vis.type_cols.count));
}

static void genSim(const SimParams& params) {
	simClear(&sim_state, &params.env, &params.type);
	genBodyTemplates(&gen_templates, &params.type, sim_state.mol_types.ptr);
	genInitialConditions(&sim_state, params.init, &gen_templates);
	genFinalize();
	genStat();
}

static void genBlankSim(SimParams* params) {
	simClear(&sim_state, &params->env, &params->type);
	genTypesAndTemplatesAndDistro(params);
	genFinalize();
	genStat();
}

struct NameAndDate {
	char* name;
	time_t date;
};
int NameAndDateCompare(const void * da, const void * db) {
	NameAndDate* a = (NameAndDate*)da;
	NameAndDate* b = (NameAndDate*)db;
	return (a->date > b->date) ? -1 : +1;
}
#pragma warning(push)
#pragma warning(disable:4706) // assignment within conditional expression
static void guiSaveLoadState(const char* ext) {
	{ // save
		gui::Separator();
		bool save = false;
		save |= gui::Button("save");
		gui::SameLine();
		save |= gui::InputText("",sim_filename_out,sizeof(sim_filename_out),ImGuiInputTextFlags_EnterReturnsTrue);
		if (save) {
			BlockTimer T("save");
			Sim S;
			S.sim_state = sim_state;
			S.sim_params = sim_params;
			S.sim_vis = sim_vis;
			SaveToFile(S, "Sim", "sim", TempStr("%s%s%s", PATH, sim_filename_out,ext));
			//SaveToFile(sim_state, "SimState", "sim_state", tempstr);
		}
	}
	{ // load

		gui::Separator();
		gui::Text("load:");
		DIR *dirp = opendir(PATH);
		int extlen = (int)strlen(ext);
		if (dirp) {
			struct dirent *entry;
			Bunch<NameAndDate> strings;
			while ((entry = readdir(dirp))) {
				if (entry->d_namlen > extlen && strstr(entry->d_name, ext)) {
					char* str = (char*)malloc(entry->d_namlen+1);
					memcpy(str, entry->d_name, entry->d_namlen-extlen);
					str[entry->d_namlen-extlen] = 0;
					struct _stat buf;
					_stat(TempStr("%s%s", PATH, entry->d_name), &buf);
					NameAndDate n;
					n.name = str;
					n.date = buf.st_mtime;
					strings.push(n);
				}
			}
			qsort(strings.ptr, strings.count, sizeof(NameAndDate), NameAndDateCompare);
			for (int i = 0; i < strings.count; ++i) {
				if (gui::Button(strings[i].name)) {
					BlockTimer T("load");
					Sim S;
					//if (LoadFromFile(sim_state, "SimState", "sim_state", tempstr)) {
					if (LoadFromFile(S, "Sim", "sim", TempStr("%s%s%s", PATH, strings[i].name, ext))) {
						sim_state = S.sim_state;
						sim_params = S.sim_params;
						sim_vis = S.sim_vis;
						genBodyTemplates(&gen_templates, &sim_params.type, sim_state.mol_types.ptr);
						genFinalize();
					}
				}
			}
			closedir(dirp);
		}  
	}
}
#pragma warning(pop)

struct SaveLoadState {
	char full_filename[2048];
	char display_filename[256];
};
void SaveLoadCallback(const FileDesc* desc, void* data) {
	SaveLoadState* state = (SaveLoadState*)data;
	if (desc->full_path_len < sizeof(SaveLoadState::full_filename) && desc->name_len < sizeof(SaveLoadState::display_filename) && ImGui::MenuItem(desc->name))  {
		strcpy(state->full_filename, desc->full_path);
		strcpy(state->display_filename, desc->name);
	}
}
/*
static void guiSaveLoadFile(const char* name, bool* is_open, const char* dir, const char* ext, bool is_save) {
	static char full_filename[2048];
	static char filename[256];
	if (gui::Begin(name, is_open)) {
		gui::Text(is_save ? "Save As:" : "Load:");
		gui::SameLine();
		gui::InputText("##inputtext", filename, sizeof(filename));
		gui::SameLine();
		gui::Button("ok");
		gui::SameLine();
		if (gui::Button("cancel")) *is_open = false;
		dirScan(dir, ext, SaveLoadCallback, full_filename);
	} gui::End();

	gui::ShowTestWindow();
}
*/
static void DragDistro(float* distro, int count, vec4* cols = NULL) {
	int active_i = -1;
	float sum_inactive = 0.f;
	for (int i = 0; i < count; ++i) {
		gui::PushID(i);
		static char buff[16];
		sprintf_s(buff, sizeof(buff), "%6.2f%%###", distro[i] * 100.f);
		//gui::PushStyleColor(ImGuiCol_FrameBg, ImColor::HSV(i/7.0f, 0.5f, 0.5f));
		//gui::PushStyleColor(ImGuiCol_FrameBgHovered, ImColor::HSV(i/7.0f, 0.6f, 0.5f));
		//gui::PushStyleColor(ImGuiCol_FrameBgActive, ImColor::HSV(i/7.0f, 0.7f, 0.5f));
		if (cols) {
			//gui::ColorConvertRGBtoHSV(cols[i].x,cols[i].y,cols[i].z,
			gui::PushStyleColor(ImGuiCol_SliderGrab, cols[i]);
			gui::PushStyleColor(ImGuiCol_SliderGrabActive, cols[i] * 1.3f);
		}
		gui::DragFloat(buff, &distro[i], 0.005f, 0.f, 1.f, "");
		if (cols) gui::PopStyleColor(2);
		if (gui::IsItemActive()) active_i = i;
		else sum_inactive += distro[i];
		gui::PopID();
	}
	if (active_i != -1) {
		float remainder = 1.f - distro[active_i];
		for (int i = 0; i < count; ++i) {
			if (i != active_i) {
				if (sum_inactive == 0.f) {
					float portion = 1.f / (count - 1);
					distro[i] = remainder * portion;
				} else {
					float portion = distro[i] / sum_inactive;
					distro[i] = remainder * portion;
				}
			}
		}
	}
}
static void guiGen(const Input& in) {
	if (gui::Begin("Initial Conditions")) {
		if (gui::Button("generate", vec2(0,30)) || (!gui::IsAnyItemActive() && in.key.press['G'])) {
			genSim(sim_params);
		}
		{ // types
			TypeParams& params = sim_params.type;
			bool t = false;

			gui::Separator();
			//gui::Text("Compounds:");
			gui::CollapsingHeader("Compounds:", ImGuiTreeNodeFlags_Leaf);
			//gui::PushStyleColor(ImGuiCol_FrameBg, vec4(1,0,0,1));
			//gui::Text("Compounds:");
			//gui::PopStyleColor();
			gui::Spacing();
			for (int i = 0; i < params.comps.count; ++i) {
				gui::PushID(i);
				CompoundDesc& desc = params.comps[i];
				gui::Text("%c:", char('A'+i));
				gui::SameLine();
				t |= gui::ColorEdit3("##color", (float*)&desc.col, ImGuiColorEditFlags_NoOptions | ImGuiColorEditFlags_NoInputs);
				gui::SameLine();
				gui::PushItemWidth(100.f);
				t |= gui::Combo("##class", &desc.comp_class, "inert\0capsid\0catalyst\0substrate\0");
				gui::PopItemWidth();
				if (desc.comp_class == CLASS_CATALYST) {
					gui::SameLine();
					char choices[64]; // 32 choices + 0 separators
					int choice_idx = 0;
					for (int c = 0; c < params.comps.count; ++c) {
						//if (c != i) { //#DODGY had to get rid of this to keep the combobox indecies ordered right. ugly
							choices[choice_idx] = char('A'+c);
							choices[choice_idx+1] = 0;
							choice_idx+=2; 
						//}
					}
					choices[choice_idx] = 0;
					gui::Text("catalyzes:");
					gui::SameLine();
					gui::PushItemWidth(50.f);
					bool pushcol = false;
					if (desc.catalyzes_compound>=0 && desc.catalyzes_compound<params.comps.count) {
						pushcol = true;
						gui::PushStyleColor(ImGuiCol_FrameBg,vec4(params.comps[desc.catalyzes_compound].col, 1.f));
					}
					int choice = desc.catalyzes_compound;
					t |= gui::Combo("##catalyzes", &choice,choices);
					desc.catalyzes_compound = (u8)choice;
					if (pushcol)
						gui::PopStyleColor();
					gui::PopItemWidth();
					gui::SameLine();
					gui::PushItemWidth(80.f);
					gui::InputFloat("chance", &desc.catalysis_chance, 0.f, 0.f, 6);
					gui::PopItemWidth();
				}
				if (desc.comp_class == CLASS_COMPOUND) {
					gui::SameLine();
					gui::Text(" contains:");
					if (desc.mol_types.count != 2) {
						desc.mol_types.clear();
						desc.mol_types.pushi(NO_TYPE,2);
					}
					for (int p = 0; p < 2; ++p) {
						gui::PushID(p);
						char choices[64]; // 32 choices + 0 separators
						int choice_idx = 0;
						for (int c = 0; c < params.comps.count; ++c) {
							//if (c != i) { //#DODGY had to get rid of this to keep the combobox indecies ordered right. ugly
								choices[choice_idx] = char('A'+c);
								choices[choice_idx+1] = 0;
								choice_idx+=2; 
							//}
						}
						choices[choice_idx] = 0;
						gui::SameLine();
						gui::PushItemWidth(50.f);
						
						bool pushcol = false;
						if (desc.mol_types[p] != NO_TYPE && desc.mol_types[p]<params.comps.count) {
							pushcol = true;
							gui::PushStyleColor(ImGuiCol_FrameBg,vec4(params.comps[desc.mol_types[p]].col, 1.f));
						}
						int choice = desc.mol_types[p];
						t |= gui::Combo("##contains",&choice,choices);
						desc.mol_types[p] = (u8)choice;
						if (pushcol)
							gui::PopStyleColor();
						gui::PopItemWidth();
						gui::PopID();
					}
				}
				if (desc.comp_class == CLASS_CAPSID) {
					gui::SameLine();
					int bond_max_count = desc.bond_max_count;
					t |= gui::InputInt("max bonds", &bond_max_count);
					desc.bond_max_count = max((u16)3, (u16)bond_max_count); 

				}
				gui::PopID();
			}
			if (!params.comps.full() && gui::Button("add")) {
				params.comps.push();
				t = true;
			}
			if (!params.comps.empty()) {
				gui::SameLine();
				if (gui::Button("remove")) {
					params.comps.pop();
					t = true;
				}
			}
			gui::Separator();

			for (int i = 0; i < params.comps.count; ++i) {
				CompoundDesc& D = params.comps[i];
				if (D.comp_class == CLASS_CATALYST && D.catalyzes_compound >= params.comps.count)
					D.catalyzes_compound = NO_TYPE;
				for (int d = 0; d < D.mol_types.count; ++d)
					if (D.mol_types[d] >= params.comps.count)
						D.mol_types[d] = NO_TYPE;
			}

			if (t) {
				genTypesAndTemplatesAndDistro(&sim_params);
			}
		}

		{ // initial conditions
			gui::InputFloat("bound radius", &sim_params.env.bound_rad);
			gui::Checkbox("heat pump", &sim_params.env.heat_pump.enable);
			if (sim_params.env.heat_pump.enable) {
				gui::SliderFloat("permeability", &sim_params.env.heat_pump.permeability, 0.001f, 1.f, "%.4f");
				gui::Checkbox("cycle enable", &sim_params.env.heat_pump.cycle_enable);
				if (sim_params.env.heat_pump.cycle_enable) {
					gui::SliderFloat("frequency", &sim_params.env.heat_pump.frequency, 0.001f, 0.5f, "%.4f"); 
					gui::SliderFloat("min avg energy", &sim_params.env.heat_pump.avg_energy_min, avg_energy_min, avg_energy_max, "%4.1f");
					gui::SliderFloat("max avg energy", &sim_params.env.heat_pump.avg_energy_max, avg_energy_min, avg_energy_max, "%4.1f");
				} else {
					gui::SliderFloat("avg energy", &sim_params.env.heat_pump.avg_energy_min, avg_energy_min, avg_energy_max, "%4.1f");
					sim_params.env.heat_pump.avg_energy_max = sim_params.env.heat_pump.avg_energy_min;
				}
			}
			gui::Spacing();
			gui::SliderFloat("spacing", &sim_params.init.spacing, 2.f, 8.f);
			gui::SliderFloat("heat", &sim_params.init.heat, 1.f, 80.f);
		}

		{ // compound distro
			gui::Spacing();
			vec4 cols[32];
			for (int i = 0; i < sim_params.type.comps.count; ++i) 
				cols[i] = vec4(sim_params.type.comps[i].col,1.f);
			DragDistro(sim_params.init.compound_distro.ptr, sim_params.init.compound_distro.count, cols);
		}
		gui::SliderFloat("capsid break threshold", &sim_params.type.capsid_break_threshold, 1.f, 1000.f, "%3.3f", 2.f);
		gui::Checkbox("one starter", &sim_params.init.one_starter);
		gui::SameLine();
		gui::InputInt("compound", &sim_params.init.starter_comp);
		sim_params.init.starter_comp = max(0, min(sim_params.type.comps.count-1, sim_params.init.starter_comp));
		gui::Checkbox("enable seed autogen", &sim_params.init.enable_seed_autogen);
		if (sim_params.init.enable_seed_autogen) {
			for (int i = 0; i < sim_params.type.comps.count; ++i) {
				gui::ColorButton("colbutton",vec4(sim_params.type.comps[i].col,1.f));	
				gui::SameLine();
				gui::PushID(i);
				TempStr str("%d", i);
				gui::InputInt(str, &sim_params.init.seed_autogen.compound_count[i]);
				gui::PopID();
			}
		}
		gui::Spacing();
		if (gui::Button("save default")) SaveToFile(sim_params, "SimParams", "sim_params", PATH"sim_params.def");
		gui::SameLine();
		if (gui::Button("load default")) {
			LoadFromFile(sim_params, "SimParams", "sim_params", PATH"sim_params.def");
			genTypesAndTemplatesAndDistro(&sim_params);
		}

		guiSaveLoadState(".run");

		static bool save_open = false;
		//if (gui::Button("Save")) save_open = true;
		//if (save_open) guiSaveLoadFile("Save Initial Conditions", &save_open, PATH, ".inc", true);
	} gui::End();
}

extern vec3 func_srgb_from_jch(vec3 lch);
static vec3 part_col(int idx) {
	//vec3 lch(0.6f, 0.4f, WangHashFloat(idx));
	vec3 lch(WangHashFloat(idx^0x8423efab, 0.2f,0.8f), WangHashFloat(idx^0x1234abcd, 0.2f,0.8f), WangHashFloat(idx));
	return func_srgb_from_jch(lch);
}

void guiInspector(const char* title, int mol_idx, int bod_idx, bool highlight) {
	gui::Text(title);
	//gui::Separator();
	if (mol_idx != NO_IDX && mol_idx >= 0 && mol_idx < sim_state.molecules.count) {
		Molecule& M = molGet(&sim_state, mol_idx);
		Gui("molecule", M, Struct());
		if (bod_idx == NO_IDX)
			bod_idx = M.body_idx;
		Body& B = bodyGet(&sim_state, bod_idx);
		pose wol_from_mol = WolFromMol(M,B);
		if (highlight) {
			primColor(0.8f,1.0f,0.3f);
			primCircle(wol_from_mol, 10.f);
		}
	}
	//gui::Spacing();
	if (bod_idx != NO_IDX && bod_idx >= 0 && bod_idx < sim_state.bodies.count) {
		Body& B = bodyGet(&sim_state, bod_idx);
		Gui("body", B, Struct());
		if (highlight) {
			float r; vec2 c;
			shellRadius(&sim_state, B, 1.0, &r, &c);
			primColor(0.8f, 0.5f, 0.1f);
			primCircle(c, r);
		}
	}
}

void guiEditor(Input& in, pose world_from_cursor) {

	if (gui::Begin("Editor")) {

		//drawColor(1,0,0);
		//drawCircleF(world_from_cursor.x,world_from_cursor.y,2.f, 16.f);

		if (gui::Button("add")) {
			edit.want_to = TOOL_ADD;
		}
		if (gui::Button("blank slate")) {
			genBlankSim(&sim_params);
		}
		gui::SliderInt("type", &edit.compound_type, 0, (int)gen_templates.count-1);
		//for (int i = 1; i < TOOL_COUNT; ++i) {
		//	if (i != 1) gui::SameLine();
		//	gui::RadioButton(TOOL_NAME(i), &edit.tool, i);
		//}

		bool updating = sim_ctrl.run || sim_ctrl.step;
		if (edit.want_to == TOOL_ADD && edit.bod_idx == NO_IDX && !gui::IsMouseHoveringAnyWindow() && !updating) {
			BodyTemplate T;
			if (edit.bond_count == 0) {
				T = gen_templates[edit.compound_type];	
			} else {
				for (int i = 0; i < edit.bond_count; ++i) {
					// ?
				}
			}
			edit.bod_idx = compoundAdd(&sim_state, T, identity());
			edit.setting = SET_NONE;
		}

		bool commit = false;
		bool delete_mol = false;

		if (in.mouse.press[MOUSE_BUTTON_1] || in.mouse.press[MOUSE_BUTTON_2])
			edit.world_from_cursor_press = world_from_cursor;

		if (gui::IsMouseHoveringAnyWindow()) delete_mol |= true;
		if (updating) delete_mol = true;
		if (!edit.want_to == TOOL_ADD) delete_mol = true;

		if (edit.bod_idx != NO_IDX) {
			Body& B = bodyGet(&sim_state, edit.bod_idx);
			if (edit.setting == SET_NONE) {
				B.wol_from_bod = world_from_cursor * edit.cursor_from_body;
				B.wol_from_bod.zw(normalize(B.wol_from_bod.zw()));
			} else if (edit.setting == SET_VEL) {
				B.dp = (world_from_cursor.xy() - edit.world_from_cursor_press.xy());
				if (in.key.press['R']) edit.setting = SET_ROT;
			} else if (edit.setting == SET_ROT) {
				B.wol_from_bod.zw(normalize(world_from_cursor.xy() - edit.world_from_cursor_press.xy()));
				if (in.key.press['V']) edit.setting = SET_VEL;
			}

			if (in.mouse.press[MOUSE_BUTTON_1]) {
				edit.setting = SET_VEL;
			} else if (in.mouse.release[MOUSE_BUTTON_1]) {
				commit |= true;
				//mol->dr = quat2d(randFloat(-0.5f,0.5f));
			} else if (in.mouse.press[MOUSE_BUTTON_2]) {
				edit.want_to = TOOL_INSPECT;
				delete_mol |= true;
			}
		}

		if (delete_mol) {
			if (edit.bod_idx != NO_IDX) {
				// #DODGY this is supremely dangerous! removing molecules invalidates all indecies past this one (for bonds, and others?) in the data
				// so this only  "works" if the molecule we are deleting is the one we just added! (so it's at the end of the array) (and if it hasn't somehow bonded to something earlier in the array already)
				// don't do this in the general case!
				for (int m = 0; m < sim_state.molecules.count; ) {
					if (sim_state.molecules[m].body_idx == edit.bod_idx) {
						sim_state.molecules.remove(m);
					} else {
						m++;
					}
				}
				bodyRem(&sim_state, edit.bod_idx); // #WRONG this pushes a free body, but we actually want to delete for real in this case
				genFinalize();
				genStat();
				genZeroVibration(); // have to zero vibration, because initial state has changed and it will be out of sync with stats
				edit.bod_idx = NO_IDX;
			}
		} else if (commit) {
			if (edit.bod_idx != NO_IDX) {
				edit.bod_idx = NO_IDX;
				genFinalize();
				genStat();
				genZeroVibration(); // have to zero vibration, because initial state has changed and it will be out of sync with stats
			}
		}

		guiSaveLoadState(".edit");
	} gui::End();

	if (gui::Begin("Inspector")) {
		gui::InputInt("mol idx", &edit.inspect_mol_idx);
		gui::InputInt("bod idx", &edit.inspect_bod_idx);
		if (edit.want_to == TOOL_INSPECT) {
			int hover_mol_idx = NO_IDX;
			if (!gui::IsMouseHoveringAnyWindow()) {
				float closest = FLT_MAX;
				for (int i = 0; i < sim_state.molecules.count; ++i) {
					const Molecule& M = sim_state.molecules[i];
					const Body& B = bodyGet(&sim_state, M.body_idx);
					float rad = molTypeGet(&sim_state, M.type).rad;
					vec2 diff = world_from_cursor.xy() - WolFromMol(M, B).xy(); // #NONGEN ignores object scale..
					float d = dot(diff,diff);
					if (d < (rad*rad) && d < closest) {
						hover_mol_idx = i;
						closest = d;
					}
				}
				if (in.mouse.press[MOUSE_BUTTON_1]) {
					edit.inspect_mol_idx = hover_mol_idx;

					if (edit.inspect_mol_idx == NO_IDX)
						edit.inspect_bod_idx = NO_IDX;
				}
			}
			if (edit.inspect_mol_idx != NO_IDX && edit.inspect_mol_idx >= 0 && edit.inspect_mol_idx < sim_state.molecules.count) 
				edit.inspect_bod_idx = molGet(&sim_state, edit.inspect_mol_idx).body_idx;
			
			gui::Columns(2, "somecolumns");	
			guiInspector("clicked",edit.inspect_mol_idx,edit.inspect_bod_idx, true);
			gui::NextColumn();
			guiInspector("hovered",hover_mol_idx, -1, false);
			gui::Columns(1);
		}
	} gui::End();
}

static void guiControl(const Input& in) {
	if (gui::Begin("Sim Control")) {
		gui::Checkbox("run", &sim_ctrl.run);
		if (!gui::IsAnyItemActive() && in.key.press[KEY_SPACE]) sim_ctrl.run = !sim_ctrl.run;
		gui::PushButtonRepeat(true);
		if (gui::Button("step") || (!gui::IsAnyItemActive() && in.key.down['S'])) {
			sim_ctrl.step = true;
		}
		gui::InputInt("substeps", &sim_ctrl.substeps);
		sim_ctrl.substeps = max(1, min(16, sim_ctrl.substeps));
		gui::InputInt("supersteps", &sim_ctrl.supersteps);
		sim_ctrl.supersteps = max(1, min(1000, sim_ctrl.supersteps));
		gui::PopButtonRepeat();
		gui::Checkbox("record data", &sim_ctrl.record_data);
		if (sim_ctrl.record_data) {
			gui::InputText("record name", (char*)sim_ctrl.record_name.ptr, sim_ctrl.record_name.maxcount*sizeof(char));
		}
	} gui::End();
}

static void guiStats(const SimStats& stats) {
	u32 t_sec = u32(stats.step_count * (1.f/60.f)); // #DODGY
	u32 t_min = u32(t_sec / 60.f);
	gui::Text("Step %lu, Time %02dm:%02ds", stats.step_count, t_min, t_sec - t_min*60);
	float total_count = 0;
	float counts[32];
	vec4 cols[32];
	for (int i = 0; i < stats.compound_counts.count; ++i) {
		//max_count = max(max_count, stats.compound_counts[i]);
		total_count += stats.compound_counts[i];
		cols[i] = vec4(sim_state.body_cols[i],1.f);
	}
	for (int i = 0; i < stats.compound_counts.count; ++i)
		counts[i] = stats.compound_counts[i] / total_count;

	for (int i = 0; i < stats.compound_counts.count; ++i) {
		static char buff[64];
		sprintf_s(buff, sizeof(buff), "%6.2f%% (%d)###", counts[i]*100.f, stats.compound_counts[i]);
		gui::PushID(i);
		gui::PushStyleColor(ImGuiCol_SliderGrab, cols[i]);
		gui::PushStyleColor(ImGuiCol_SliderGrabActive, cols[i]*1.3f);
		gui::DragFloat(buff, &counts[i], 0.005f, 0.f, 1.f, "");
		gui::PopStyleColor(2);
		gui::PopID();
	}

	gui::Spacing();
	gui::PlotHistogram("speed distro", stats.speed_hist.ptr, stats.speed_hist.maxcount, 0, 0, 0.f, stats.bucket_count_max, vec2(400,200));
	//gui::Text("sum mass - %f", mass_sum);
	gui::Text("%d molecules", stats.mol_count);
	gui::Text("%d bodies", stats.bod_count);
	gui::Text("sum linear momentum:          (%f,%f) |%f|", stats.lin_momentum_sum.x,stats.lin_momentum_sum.y, length(stats.lin_momentum_sum));
	gui::Text("sum rotational momentum:       %f", stats.rot_momentum_sum);
	gui::Text("sum linear kinetic energy:     %f", stats.lin_kinetic_sum);
	gui::Text("sum rotational kinetic energy: %f", stats.rot_kinetic_sum);
	gui::Text("sum kinetic energy:            %f", stats.lin_kinetic_sum + stats.rot_kinetic_sum);
	gui::Text("sum vibration:                 %f", stats.vibration_sum);
	gui::Text("avg energy:                    %f", (stats.lin_kinetic_sum+stats.rot_kinetic_sum+stats.vibration_sum) / stats.mass_sum);
	//gui::Text("mols %dKb\n", sim_state.molecules.bytes() / 1024);
	//gui::Text("bodies %dKb\n", sim_state.bodies.bytes() / 1024);
	//gui::Text("hash %dMb\n", (cells.bytes() / (1024*1024)));
		
	
}
void drawStats() {

	static float speed_range = 250;

	static SimStats stats_prev;
	SimStats stats;
	stats.speed_range = speed_range;
	computeStats(stats);

	StatsReq req;
	req.shell_contents = true;
	req.stats = &stats;
	simStatsOnly(&sim_state, &timer, req);

	const SimParams& SP = sim_params; //#DODGY need to pull this from saved state, not current editor state

	int catalyst_type_count_max = 0;
	for (u8 d = 0; d < SP.type.comps.count; ++d) {
		const CompoundDesc& D = SP.type.comps[d];
		if (D.catalyzes_compound != NO_TYPE) {
			catalyst_type_count_max += 1;
		}
	}

	int shell_status_counts[SHELL_STATUS_count] = { 0 };
	if (gui::Begin("COUNT")) {
		for (int i = 0; i < stats.shell_guts.count; ++i) {
			gui::PushID(i);

			const BodyHist& hist = stats.shell_guts[i];
			u8 hi = 0;
			gui::Text("body_idx %d", hist.body_idx);
			gui::Indent();
			gui::Text("Contents:");
			gui::PushID("hist");
			for (u8 h = 0; h < (u8)hist.body_type_masses.count; ++h) {
				gui::PushID(h);
				float c = hist.body_type_masses[h];	
				if (c > 0) {
					if (hi == 0) gui::SameLine(); // for break after text
					if (hi != 0) gui::SameLine();
					gui::ColorButton("colbutton",vec4(simBodyCol(&sim_state, h), 1.f));	
					gui::SameLine();
					gui::Text("%02.0f", c);
					++hi;
				}
				gui::PopID();
			}
			gui::PopID();
			gui::Unindent();


			int catalyst_active_count[BODY_TYPE_COUNT_MAX] = {0};
			int catalyst_inactive_count[BODY_TYPE_COUNT_MAX] = {0};
			float total_mass = 0.f;
			for (u8 h = 0; h < (u8)hist.body_type_masses.count; ++h) {
				const CompoundDesc& D = SP.type.comps[h];
				float compound_mass = 0.f;
				if (D.comp_class == CLASS_COMPOUND) {
					for (u8 ti = 0; ti < D.mol_types.count; ++ti) {
						if (D.mol_types[ti] != NO_TYPE) {
							compound_mass += molTypeGet(&sim_state, D.mol_types[ti]).mass;
						}
					}
				} else {
					compound_mass += molTypeGet(&sim_state, h).mass;
				}
				total_mass += hist.body_type_masses[h];
				int c = round_nearest(hist.body_type_masses[h] / compound_mass);
				if (D.comp_class == CLASS_CATALYST) { 
					// count active catalysts
					catalyst_active_count[h] += c;
				} else {
					// if this is a substrate, add it's count to each catalyst type that it contains
					for (u8 ti = 0; ti < D.mol_types.count; ++ti) {
						if (D.mol_types[ti] != NO_TYPE) {
							u8 t = D.mol_types[ti];
							if (SP.type.comps[t].comp_class == CLASS_CATALYST) {
								catalyst_inactive_count[t] += c;
							}
						}
					}
				}
			}

			hi = 0;
			gui::Indent();
			gui::Text("Active:  ");
			gui::PushID("active");
			for (int h = 0; h < BODY_TYPE_COUNT_MAX; ++h) {
				int c = catalyst_active_count[h];
				if (c > 0) {
					if (hi == 0) gui::SameLine(); // for break after text
					if (hi != 0) gui::SameLine();
					gui::ColorButton("colbutton",vec4(simBodyCol(&sim_state, h), 1.f));	
					gui::SameLine();
					gui::Text("%02d", c);
					++hi;
				}
			}
			gui::PopID();
			gui::Unindent();

			hi = 0;
			gui::Indent();
			gui::Text("Inactive:");
			gui::PushID("inactive");
			for (int h = 0; h < BODY_TYPE_COUNT_MAX; ++h) {
				int c = catalyst_inactive_count[h];
				if (c > 0) {
					if (hi == 0) gui::SameLine(); // for break after text
					if (hi != 0) gui::SameLine();
					gui::ColorButton("colbutton",vec4(simBodyCol(&sim_state, h), 1.f));	
					gui::SameLine();
					gui::Text("%02d", c);
					++hi;
				}
			}
			gui::PopID();
			gui::Unindent();

			gui::Indent();
			gui::Text("Status:");
			gui::SameLine();

			bool any_active = false;
			bool any_present = false;
			bool all_present = true;
			bool all_active = true;
			for (int a = 0; a < typeGetCount(&sim_state); ++a) {
				if (SP.type.comps[a].comp_class == CLASS_CATALYST) {
					if (catalyst_active_count[a] > 0)
						any_active = true;
					if (catalyst_active_count[a] > 0 || catalyst_inactive_count[a] > 0)
						any_present = true;
					if (catalyst_active_count[a] == 0)
						all_active = false;
					if (catalyst_active_count[a] == 0 && catalyst_inactive_count[a] == 0)
						all_present = false;
				}
			}
			int status;
			if (all_active) {
				status = SHELL_STATUS_ripe;
			} else if (any_active && all_present) {
				status = SHELL_STATUS_fertile;
			} else if (all_present) {
				status = SHELL_STATUS_dormant;
			} else if (any_present) {
				status = SHELL_STATUS_inert;
			} else if (total_mass > 0.f) {
				status = SHELL_STATUS_capsid;
			} else {
				status = SHELL_STATUS_empty;
			}
			shell_status_counts[status] += 1;
			gui::Text("%s", SHELL_STATUS_NAME(status));
			gui::PushID("status");
			gui::SameLine();
			gui::ColorButton("colbutton",vec4(simShellStatusCol(status), 1.f));	
			gui::PopID();
			gui::Unindent();

			gui::PopID();

			if (sim_vis.show_shell_type) {
				const Body& B = bodyGet(&sim_state, hist.body_idx);
				float r; vec2 c;
				shellRadius(&sim_state, B, 1.0, &r, &c);
				primColor(vec4(simShellStatusCol(status), 1.f)*0.4f);
				primCircle(c, r);
			}
		}
	}
	gui::End();

	if (gui::Begin("Shell Status")) {
		for (int s = 0; s < SHELL_STATUS_count; ++s) {
			gui::PushStyleColor(ImGuiCol_Text, vec4(simShellStatusCol(s), 1.f));
			gui::Text("%s: %d", SHELL_STATUS_NAME(s), shell_status_counts[s]);
			gui::PopStyleColor(1);
		}
		gui::Spacing();
		gui::Checkbox("highlight", &sim_vis.show_shell_type);
	} gui::End();

	
	stats_prev = stats;

	vec2 lin_total_momentum = stats.lin_momentum_sum +  sim_state.lin_bound_momentum;
	float rot_total_momentum = stats.rot_momentum_sum +  sim_state.rot_bound_momentum;

	if (gui::Begin("Stats")) {
		//gui::SetWindowFontScale(2.f);
		gui::Text("==== Current Stats =====");
		guiStats(stats);
		gui::SliderFloat("speed range", &speed_range, 0.f, 200.f);
		gui::Spacing();
		gui::Text("==== Initial Stats =====");
		guiStats(sim_state.init_stats);
		gui::Text("===== Losses =====");
		float lin_loss = sim_state.init_stats.lin_kinetic_sum - stats.lin_kinetic_sum;
		gui::Text("linear kinetic energy:           %f", lin_loss);
		gui::Text("dissapated linear kinetic energy %f", sim_state.dissapated_lin_kinetic_energy);
		gui::Text("                           delta %f", lin_loss -  sim_state.dissapated_lin_kinetic_energy);
		float rot_loss = sim_state.init_stats.rot_kinetic_sum - stats.rot_kinetic_sum;
		gui::Text("rotational kinetic energy:           %f", rot_loss);
		gui::Text("dissapated rotational kinetic energy %f",  sim_state.dissapated_rot_kinetic_energy);
		gui::Text("                               delta %f", rot_loss -  sim_state.dissapated_rot_kinetic_energy);
		gui::Text("total energy:            %f", lin_loss + rot_loss);
		gui::Text("vibration track:         %f", (lin_loss + rot_loss) - stats.vibration_sum);
		float acc = (sim_state.init_stats.lin_kinetic_sum+sim_state.init_stats.rot_kinetic_sum) - (stats.lin_kinetic_sum+stats.rot_kinetic_sum+stats.vibration_sum);
		gui::Text("accounting %f (%f%%)", acc, acc / (sim_state.init_stats.lin_kinetic_sum+sim_state.init_stats.rot_kinetic_sum)*100.f);
		gui::Text("dissapated total energy: %f",  sim_state.dissapated_lin_kinetic_energy +  sim_state.dissapated_rot_kinetic_energy);
		float delta_loss = (lin_loss+rot_loss) - ( sim_state.dissapated_lin_kinetic_energy+ sim_state.dissapated_rot_kinetic_energy);
		gui::Text("                   delta %f", delta_loss);
		gui::Text("        percent of total %f", delta_loss / (sim_state.init_stats.lin_kinetic_sum + sim_state.init_stats.rot_kinetic_sum) * 100.f);
		gui::Text("     diff with vibration %f", ( sim_state.dissapated_lin_kinetic_energy+ sim_state.dissapated_rot_kinetic_energy) - stats.vibration_sum);
		gui::Text("sum lin momentum (with bounds): (%f,%f) |%f|", lin_total_momentum.x,lin_total_momentum.y, length(lin_total_momentum));
		vec2 delta_lin_momentum = sim_state.init_stats.lin_momentum_sum - lin_total_momentum;
		gui::Text("                        delta: (%f,%f) |%f|", delta_lin_momentum.x,delta_lin_momentum.y, length(delta_lin_momentum));
		gui::Text("                      percent: (%f,%f) |%f|", delta_lin_momentum.x/sim_state.init_stats.lin_momentum_sum.x*100.f,delta_lin_momentum.y/sim_state.init_stats.lin_momentum_sum.y*100.f, length(delta_lin_momentum)/length(sim_state.init_stats.lin_momentum_sum)*100.f);
		float delta_rot_momentum = sim_state.init_stats.rot_momentum_sum - rot_total_momentum;
		gui::Text("sum rot momentum (with bounds): %f", rot_total_momentum);
		gui::Text("                         delta: %f", delta_rot_momentum);
		gui::Text("                       percent: %f", delta_rot_momentum/sim_state.init_stats.rot_momentum_sum*100.f);
	} gui::End();

	if (gui::Begin("Thermometer")) {
		float total_energy_sum = stats.lin_kinetic_sum + stats.rot_kinetic_sum + stats.vibration_sum;
		float total_energy_avg = stats.mass_sum > 0.f ? total_energy_sum/stats.mass_sum : 0.f;
		gui::PushStyleColor(ImGuiCol_FrameBg, vec4(1.f,0.9f,0.9f,1.f)*0.7f);
		gui::PushStyleColor(ImGuiCol_SliderGrab, vec4(1.f,0.3f,0.3f,1.f));
		gui::PushStyleColor(ImGuiCol_SliderGrabActive, vec4(1.f,0.3f,0.3f,1.f)*1.3f);
		gui::VSliderFloat("##temp", vec2(45,500), &total_energy_avg, avg_energy_min, avg_energy_max, "%4.1f");
		gui::PopStyleColor(3);

		SimState* S = &sim_state;
		if (S->heat_pump.enable) {
			gui::SameLine();
			gui::VSliderFloat("##ctrlmin", vec2(45,500), &S->heat_pump.avg_energy_min, avg_energy_min, avg_energy_max, "%4.1f");
			if (S->heat_pump.cycle_enable) {
				gui::SameLine();
				gui::VSliderFloat("##ctrlmax", vec2(45,500), &S->heat_pump.avg_energy_max, avg_energy_min, avg_energy_max, "%4.1f");
			} else {
				S->heat_pump.avg_energy_max = S->heat_pump.avg_energy_min;
			}
		}

		gui::Checkbox("heat pump", &S->heat_pump.enable);
		gui::SliderFloat("permeability", &S->heat_pump.permeability, 0.001f, 1.f, "%.4f");
		gui::Checkbox("cycling", &S->heat_pump.cycle_enable);
		if (S->heat_pump.cycle_enable) {
			gui::SliderFloat("cycle frequency", &S->heat_pump.frequency, 0.001f, 0.5f, "%.4f"); 
		}

		//gui::Text("%4.1f", total_energy_avg);

		//gui::SameLine();
		//gui::Text("%4.1f", sim_ctrl.target_avg_energy);
		//gui::Spacing();
		//gui::Text("linear:     %6.2f%%", stats.lin_kinetic_sum / total_energy_sum*100.f);
		//gui::Text("rotational: %6.2f%%", stats.rot_kinetic_sum / total_energy_sum*100.f);
		//gui::Text("vibration:  %6.2f%%", stats.vibration_sum / total_energy_sum*100.f);
	} gui::End();

	if (sim_vis.show_velocity) {
		primColor(0.8f,0.3f,0.1f);
		primCircle(vec2(0.f), MOL_RAD*0.25f);
		primLine(vec2(0.f),stats.lin_momentum_sum,MOL_RAD*0.1f);
		primWedge(pose(vec2(stats.lin_momentum_sum), normalize(stats.lin_momentum_sum)), vec2(1.f,0.5f)*MOL_RAD*0.5f);	
		primColor(0.1f,0.3f,0.8f);
		primLine(vec2(0.f),sim_state.lin_bound_momentum,MOL_RAD*0.1f);
		primWedge(pose(sim_state.lin_bound_momentum, normalize(sim_state.lin_bound_momentum)), vec2(1.f,0.5f)*MOL_RAD*0.5f);	
		
		primColor(0.2f,0.8f,0.1f);
		primLine(vec2(0.f),lin_total_momentum,MOL_RAD*0.1f);
		primWedge(pose(lin_total_momentum, normalize(lin_total_momentum)), vec2(1.f,0.5f)*MOL_RAD*0.5f);	
	}

	// save into a folder
	// sample every sim second
	// save state every 30 sim minutes
	// wall time, sim time, heater temp, sim temp, substrate A, substrate B, catalyst A, catalyst B, free capsid, bound capsid, shelled capsid, total shells, ripe, fertile, dormant, inert, capsid, empty

	SimState* S = &sim_state;
	const char* record_path = "data/autogen/record1/";
	const char* data_pathfile = "data/autogen/record1/record.csv";
	if (sim_ctrl.record_data && sim_ctrl.run || sim_ctrl.step) {
		if (S->time == 0.f) {
			FILE* f = fopen(data_pathfile, "w");
			fprintf(f, "sim time,heater avg energy,bodies avg energy,catalyst A,catalyst B,substrate,free capsid,bound capsid,shelled capsid,ripe,fertile,dormant,inert,capsid,empty\n");
			fclose(f);
		}
		if (S->time == 0.f || (S->time - sim_ctrl.time_of_last_datum) > 1.f) {
			FILE* f = fopen(data_pathfile, "a");
			u8 capsid_mol_type = NO_TYPE;
			for (int i = 0; i < sim_state.mol_types.count; ++i) {
				if (sim_state.mol_types[i].bond_max_count > 0) {
					capsid_mol_type = i;
					break;
				}
			}
			u8 capsid_substrate_compound = NO_TYPE;
			for (int i = 0; i < sim_state.mol_types_inside.count; ++i) {
				u32 ti = sim_state.mol_types_inside[i];
				if ((ti&(1<<capsid_mol_type)) && (ti&(~(1<<capsid_mol_type)))) { // if it has capsid, and also something other than it
					capsid_substrate_compound = i;
					break;
				}
			}
			u8 catalyst_compounds[2] = {NO_TYPE,NO_TYPE};
			for (int i = 0; i < sim_state.mol_types_inside.count; ++i) {
				if (sim_state.catalysis_chances[i] != 1.f) { // #DODGY! 1.f is just a bad sentinel value...
					if (catalyst_compounds[0] == NO_TYPE) catalyst_compounds[0] = i;
					else if (catalyst_compounds[1] == NO_TYPE) catalyst_compounds[1] = i;
					else break;
				}
			}
			fprintf(f, "%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
				S->time,
				simHeatPumpAvgEnergy(S),
				(stats.lin_kinetic_sum+stats.rot_kinetic_sum+stats.vibration_sum) / stats.mass_sum,
				catalyst_compounds[0] != NO_TYPE ? stats.compound_counts[catalyst_compounds[0]] : 0,
				catalyst_compounds[1] != NO_TYPE ? stats.compound_counts[catalyst_compounds[1]] : 0,
				capsid_substrate_compound != NO_TYPE ? stats.compound_counts[capsid_substrate_compound]/2 : 0, //there are 2 molecules in each substrate
				stats.free_capsid_count,
				stats.bound_capid_count,
				stats.shelled_capsid_count,
				shell_status_counts[SHELL_STATUS_ripe],
				shell_status_counts[SHELL_STATUS_fertile],
				shell_status_counts[SHELL_STATUS_dormant],
				shell_status_counts[SHELL_STATUS_inert],
				shell_status_counts[SHELL_STATUS_capsid],
				shell_status_counts[SHELL_STATUS_empty]);
			sim_ctrl.time_of_last_datum = S->time;
			fclose(f);
		}
		if (S->time == 0.f || (S->time - sim_ctrl.time_of_last_state) > (60.f*30.f)) {
			char pathfile[128];
			sprintf(pathfile, "%sstate_at_%8.0f.run", record_path, S->time);
			Sim saveS;
			saveS.sim_state = sim_state;
			saveS.sim_params = sim_params;
			saveS.sim_vis = sim_vis;
			SaveToFile(saveS, "Sim", "sim", pathfile);
			sim_ctrl.time_of_last_state = S->time;
		}
	}
}
static void guiColorscheme() {
	if (gui::Begin("Color Scheme")) {
		gui::Text("Compounds:");
		//TextFramed("compounds");
		for (int i = 0; i < sim_state.body_cols.count; ++i) {
			gui::PushID(i);
			gui::ColorEdit3("##color", (float*)&sim_state.body_cols[i], ImGuiColorEditFlags_NoOptions);
			gui::PopID();
		}
		gui::Text("Environment:");
		//TextFramed("environment");
		gui::ColorEdit3("background", (float*)&sim_vis.bg_col);
		gui::ColorEdit3("foreground", (float*)&sim_vis.fg_col);
	} gui::End();
};

static void Autoload() {
	if (0 /*load all*/) {

	} else {
		if (!LoadFromFile(sim_params, "SimParams", "sim_params", PATH"sim_params.last"))
			LoadFromFile(sim_params, "SimParams", "sim_params", PATH"sim_params.def");
		if (!LoadFromFile(sim_vis, "SimVis", "sim_vis", PATH"sim_vis.last")) {
			LoadFromFile(sim_vis, "SimVis", "sim_vis", PATH"sim_vis.def");
			wol_from_cam = sim_vis.wol_from_cam;
		}

		//LoadFromFile(sim_state, "SimState", "sim_state", PATH"sim_state.last");
	}
}
static void Autosave() {
	SaveToFile(sim_params, "SimParams", "sim_params", PATH"sim_params.last");
	SaveToFile(sim_vis, "SimVis", "sim_vis", PATH"sim_vis.last");
	//SaveToFile(sim_state, "SimState", "sim_state", PATH"sim_state.last");
}

LAY_EGG(autogen)
{
	if (args.start) {
		primInit();
		Autoload();
		genBlankSim(&sim_params);
		return;
	}
	if (args.stop) {
		Autosave();
		return;
	}
	float aspect = float(args.app.res_x)/float(args.app.res_y);

	//drawClear(0.8f, 0.8f, 0.9f);
	//primColor(0.2f, 0.55f, 0.75f);

	drawViewport(0, 0, args.app.res_x, args.app.res_y);
	drawClear(sim_vis.bg_col);

	// #WRONG have to draw bounds first, cause no layers yet
	primColor(sim_vis.fg_col);
	primCircle(0.f,0.f,sim_state.bound_rad);
	if (sim_vis.show_center) {	
		primColor(sim_vis.fg_col*0.8f);
		primCircle(0.f,0.f,200.f);
	}

	pose wol_from_cam_overview = trans(vec2(0,0), 0.f, sim_state.bound_rad);
	pose cam_from_cursor = trans(mulcomp(vec2(args.in.mouse.x*2.f-1.f,args.in.mouse.y*2.f-1.f),vec2(aspect,1.f)), 0.f, 1.f);
	pose world_from_cursor = wol_from_cam * cam_from_cursor;

	if (!gui::IsMouseHoveringAnyWindow() && !gui::IsAnyItemActive()) {
		if (args.in.key.press['0']) sim_vis.wol_from_cam = wol_from_cam_overview;
		for (int i = 0; i < 9; ++i) {
			if (args.in.key.press['1'+i]) {
				vec2 xy = sim_vis.wol_from_cam.xy();
				if (edit.inspect_mol_idx != NO_IDX) {
					Molecule& M = molGet(&sim_state, edit.inspect_mol_idx);
					Body& B = bodyGet(&sim_state, M.body_idx);
					pose wol_from_mol = WolFromMol(M,B);
					xy = wol_from_mol.xy();
				} else {

					xy = world_from_cursor.xy();
				}
				sim_vis.wol_from_cam = trans(xy, 0.f, sim_state.bound_rad*((i+1)*0.1f));
			}
		}
	} 
	wol_from_cam = lerp(wol_from_cam, sim_vis.wol_from_cam, 0.1f);

	pose cam_from_world = ~wol_from_cam;

	guiEditor(args.in, world_from_cursor);
	guiControl(args.in);

	timer.stop();
	timer.reset();
	timer.start("loop");

	timer.start("stats");
	drawStats();
	timer.stop();

	if (sim_ctrl.run || sim_ctrl.step) {
		timer.start("update");
		for (int s = 0; s < sim_ctrl.supersteps; ++s) {
			timer.start("vis reset");
			for (int i = 0; i < mols_vis.count; ++i) {
				mols_vis[i].is_colliding = false;
			}
			timer.stop();
			timer.start("step");
			StatsReq sreq;
			//sreq.shell_contents = true;
			simUpdate(&sim_state, sim_state.dt, sim_ctrl.substeps, &timer, sreq);
			timer.stop();
		}
		timer.stop();
		sim_ctrl.step = false;
	} else {
		//StatsReq sreq;
		//sreq.shell_contents = true;
		//simStatsOnly(&sim_state,  &timer, sreq);
	} 
	
	timer.start("draw");
	
	if (sim_vis.show_cells) {
		drawCells();
	}

	const Bunch<Molecule>& mols = sim_state.molecules;
	const Bunch<Body>& bods = sim_state.bodies;

	timer.start("submit draw");
	for (int i = 0; i < mols.count; ++i) {
		const Molecule& M = mols[i];
		const Body& B = bodyGet(&sim_state, M.body_idx);
		float rad = molTypeGet(&sim_state, M.type).rad;
		if (sim_vis.filter_type.count > M.type && !sim_vis.filter_type[M.type]) continue;

		if (sim_vis.show_collisions && mols_vis[i].is_colliding)
			primColor(1,0,0);
		else if (sim_vis.color_body_idx)
			primColor(M.body_idx == NO_IDX ? vec3(1.f,0.f,1.f) : part_col(M.body_idx));
		else if (sim_vis.color_mol_idx)
			primColor(part_col(i));
		else {
			primColor(simBodyCol(&sim_state, B.type));
		}

		pose wol_from_mol = WolFromMol(M,B);
		primCircle(wol_from_mol, rad);

	}
	// draw the vis on top
	// #OPT gross but best option until we have layers
	for (int i = 0; i < mols.count; ++i) {
		if (sim_vis.show_bonds || sim_vis.show_sites || sim_vis.show_mol_rotation) {
			const Molecule& M = mols[i];
			const Body& B = bodyGet(&sim_state, M.body_idx);
			float rad = molTypeGet(&sim_state, M.type).rad;
			pose wol_from_mol = WolFromMol(M,B);

			if (sim_vis.show_mol_rotation) {
				primColor(1.f);
				primLine(wol_from_mol*vec2(-rad,0.f),wol_from_mol*vec2(+rad,0.f), rad*0.1f); 
			}
			if (sim_vis.show_sites) {
				const MoleculeType& type = molTypeGet(&sim_state, M.type);
				for (int s = 0; s < type.sites.count; ++s) {
					if (siteIsFull(M,s)) continue;
					pose wol_from_sit = wol_from_mol * type.sites[s].mol_from_sit;
					if (type.sites[s].gender == GENDER_BOTH) primColor(0.4f,1.f,0.4f);
					if (type.sites[s].gender == GENDER_FEMALE) primColor(1.0f,0.4f,0.4f);
					if (type.sites[s].gender == GENDER_MALE) primColor(0.f,0.4f,1.f);
					primCircle(wol_from_sit, rad*1.1f);
					primColor(0,0,0);
					primCircle(wol_from_sit, rad*BIND_POS_THRESH);
				}
			}
			if (sim_vis.show_bonds) {
				int bond_idx = getBond(M, 0);
				if (bond_idx != NO_IDX) {
					const Molecule& O = molGet(&sim_state, bond_idx);
					//primColor(0);
					primColor(simBodyCol(&sim_state, B.type));
					//primLine(wol_from_mol.xy(), WolFromMol(O, bodyGet(&sim_state, O.body_idx)).xy(), rad*2.f);//*0.1f);
					primColor(vec3(1.0f));
					primLine(wol_from_mol.xy(), WolFromMol(O, bodyGet(&sim_state, O.body_idx)).xy(), rad*0.1f);
				}
			}
		}
	}
	if (sim_vis.show_velocity) {
		for (int i = 0; i < bods.count; ++i) {
			const Body& B = bods[i];
			if (!B.valid) continue;
			vec2 wol_from_com = B.wol_from_bod * B.bod_from_com;
			if (sim_vis.show_velocity) {
				primColor(0,0,0);
				primCircle(wol_from_com, MOL_RAD*0.25f);
				primLine(wol_from_com,wol_from_com+B.dp,MOL_RAD*0.1f);
				primWedge(pose(vec2(wol_from_com+B.dp), normalize(B.dp)), vec2(1.f,0.5f)*MOL_RAD*0.5f);
			}
		}
	}
	if (sim_vis.show_com) {// || sim_vis.show_rings) {
		for (int i = 0; i < bods.count; ++i) {
			if (!bods[i].valid) continue;
			const Body& B = bods[i];
			vec2 p = B.wol_from_bod * B.bod_from_com;
			if (sim_vis.show_com) {
				primColor(1,1,1);
				primCircle(p, 0.1f*MOL_RAD);
			}
			/*
			if (sim_vis.show_rings) {
				const Molecule& M = molGet(&sim_state, B.some_mol_idx);
				if (molTypeGet(&sim_state, M.type).bond_max_count > 0 && B.mol_count > (molTypeGet(&sim_state, M.type).bond_max_count/2u)) {
					float r; vec2 c;
					shellRadius(&sim_state, B, 1.0, &r, &c);
					primColor(0.5f, 0.5f, 0.5f, 0.5f);
					primCircle(c, r);
				}
			}
			*/
		}
	}
	timer.stop();
	/*
	primColor(0.3f,0.9f,0.1f,0.2f);
	primCircle(debug_c, debug_r);
	gui::Begin("COUNT");
	gui::SliderFloat("rad", &debug_r, 0.f, 10.f);
	if (args.in.mouse.down[MOUSE_BUTTON_2])
		debug_c = world_from_cursor.xy();
	gui::End();
	*/

	timer.start("drawcall");
	primDraw(args.app.res_x, args.app.res_y, cam_from_world);
	primClear();
	timer.stop();

	timer.stop();

	if (gui::Begin("timer")) {
		int count = timer.getResults(timer_results, TIMER_TREE_MAX_TIMERS);
		for (int i = 0; i < count; ++i) {
			for (int k = 0; k < timer_results[i].depth; ++k) gui::Indent();
			gui::Text("%s : %3.2f", timer_results[i].name, timer_results[i].time);
			for (int k = 0; k < timer_results[i].depth; ++k) gui::Unindent();
		}
	} gui::End();

	guiGen(args.in);
	Gui("Sim State", sim_state);
	Gui("Sim Vis", sim_vis);
	guiColorscheme();
	if (gui::Begin("Sim Vis")) {
		if (gui::Button("save default")) SaveToFile(sim_vis, "SimVis", "sim_vis", PATH"sim_vis.def");
		gui::SameLine();
		if (gui::Button("load default")) LoadFromFile(sim_vis, "SimVis", "sim_vis", PATH"sim_vis.def");
	} gui::End();

}
  