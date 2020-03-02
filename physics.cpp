
#include "core/vec2.h"
#include "core/pose.h"
#include "core/cpu_timer.h"
#include <unordered_set> // cleanup

// #DODGY
#include "imgui/imgui.h"

#include "autogen.h"

#include "physics.h"

struct CollPair {
	int a,b;
	bool operator==(const CollPair& rhs) const { return a==rhs.a && b==rhs.b; }
};
struct SitePair {
	int a,b;
};


vec3 simBodyCol(SimState* S, u8 type) {
	if (type != u8(-1) && !S->body_cols.empty())
		return S->body_cols[type];
	else
		return vec3(1.f,0.f,1.f);
}
vec3 simShellStatusCol(int status) {
	switch (status) {
	case SHELL_STATUS_empty: return vec3(1.f);
	case SHELL_STATUS_capsid: return vec3(0.4f,0.6f,1.f);
	case SHELL_STATUS_inert: return vec3(0.1f,1.0f,0.9f);
	case SHELL_STATUS_dormant: return vec3(0.2f,1.0f,0.1f);
	case SHELL_STATUS_fertile: return vec3(1.f,1.f,0.1f);
	case SHELL_STATUS_ripe: return vec3(1.0,0.2,0.1f);
	case SHELL_STATUS_count: return vec3(0.f);
	}
	return vec3(1.f,0.f,1.f);
}
u8 simTypeOfCompound(SimState* S, u32 types_inside) {
	for (u8 i = 0; i < (u8)S->mol_types_inside.count; ++i)
		if (S->mol_types_inside[i] == types_inside)
			return i;
	return (u8)-1;
}
float simHeatPumpAvgEnergy(SimState* S) {
	return lerp(S->heat_pump.avg_energy_min, S->heat_pump.avg_energy_max, (1.f + sin(S->time*S->heat_pump.frequency*TAU))*0.5);
}

const Molecule& molGet(const SimState* S, int idx) {
	return S->molecules[idx];
}
Molecule& molGet(SimState* S, int idx) {
	return S->molecules[idx];
}
const Body& bodyGet(const SimState* S, int idx) {
	return S->bodies[idx];
}
Body& bodyGet(SimState* S, int idx) {
	return S->bodies[idx];
}
u8 typeGetCount(const SimState* S) {
	return (u8)S->mol_types.count;
}
const MoleculeType& molTypeGet(const SimState* S, u8 idx) {
	return S->mol_types[idx];
}
pose WolFromMol(const Molecule& M, const Body& B) {
	if (M.bod_from_mol == identity()) return B.wol_from_bod; // this shaves 1ms from the loop, when it's all single molecules 
	else return B.wol_from_bod * M.bod_from_mol;
	//return B.wol_from_bod * M.bod_from_mol;
}
vec2 WolFromMolPos(const Molecule& M, const Body& B) {
	return B.wol_from_bod * M.bod_from_mol.xy();
}
void shellRadius(const SimState* S, const Body& B, float mol_radii_to_add, float* r, vec2* c) {
	//#OPT it is faster if these were inlined, so if the caller has these already accessed we shoud really pass them in
	const Molecule& M = molGet(S, B.some_mol_idx);
	const MoleculeType& T_m = molTypeGet(S, M.type);

	int mol_idx_l = B.some_mol_idx;
	int mol_idx_r = getBond(M, 0);
	if (mol_idx_r == NO_IDX) {
		mol_idx_r = mol_idx_l;
		mol_idx_l = getBond(M, 1);
	}
	if (mol_idx_l != NO_IDX && mol_idx_r != NO_IDX) {
		pose wol_from_mol_l = WolFromMol(molGet(S, mol_idx_l), B);
		pose wol_from_mol_r = WolFromMol(molGet(S, mol_idx_r), B);
		vec2 a = wol_from_mol_l.xy();
		vec2 b = wol_from_mol_r.xy();
		float t = TAU/T_m.bond_max_count;
		vec2 a_to_b = b-a;
		float L = length(a_to_b);
		float D = (L*0.5f) / tan(t*0.5f);
		*c = (b+a)*0.5f + normalize(perp(a_to_b))*D;
		*r = sqrt(sqr(L*0.5f) + sqr(D)) + T_m.rad * mol_radii_to_add;
	} else {
		*c = WolFromMolPos(M, B);
		*r = T_m.rad;
	}
}
void bodyAddMass(Body& B, vec2 bod_from_mass_pos, float mass, float moi, int mass_mol_count) {
	vec2 new_com = (B.bod_from_com*B.mass + bod_from_mass_pos*mass) / (B.mass + mass);
	B.moi = (B.moi + lengthsqr(B.bod_from_com - new_com)*B.mass) + (moi + lengthsqr(bod_from_mass_pos - new_com)*mass);
	B.mass = B.mass + mass;
	B.bod_from_com = new_com;
	B.mol_count += mass_mol_count;
}
void bodyAddFree(SimState* S) {
	S->bodies.push(); 
	S->bodies.top().valid = false;
	S->free_bodies.push((int)S->bodies.count-1); 
}
int bodyAdd(SimState* S) {
	if (S->free_bodies.count > 0) {
		int idx = S->free_bodies.top();
		S->free_bodies.pop();
		S->bodies[idx] = Body(); // also sets valid true
		return idx;
	} else {
		assert(false); // we should not be changing the array size after the initial generat
		return -1;
	}
}
void bodyRem(SimState* S, int idx) {
	bodyGet(S, idx).valid = false;
	S->free_bodies.push(idx);
}



/////////////////////////////////

int getBond(const Molecule& A, int site_idx) {
	switch(site_idx) {
		case 0: return A.r_bond;
		case 1: return A.l_bond;
		default: return -1;
	}
}
void setBond(Molecule& A, int site_idx, int other_idx) {
	switch(site_idx) {
		case 0: A.r_bond = other_idx; break;
		case 1: A.l_bond = other_idx; break;
		default: break;
	}
}
void fillSite(Molecule& A, int site_idx, int other_idx) {
	A.sites_taken |= (1<<site_idx);
	setBond(A, site_idx, other_idx);
}
void emptySite(Molecule& A, int site_idx) {
	A.sites_taken &= ~(1<<site_idx);
	setBond(A, site_idx, NO_IDX);
}
int siteIsFull(const Molecule& A, int site_idx) {
	return A.sites_taken & (1<<site_idx);
}
float computeLinearKineticEnergy(const Body* B) {
	return 0.5f * B->mass * lengthsqr(B->dp);
}
float computeRotationalKineticEnergy(const Body* B) {
	return 0.5f * B->moi * sqr(B->dr);
}
float computeRotationalMomentum(const Body* B) {
	vec2 wol_from_com = B->wol_from_bod * B->bod_from_com;
	return det(wol_from_com, B->dp*B->mass) + B->moi*B->dr;
}

float simRandFloat(SimState* S) {
	S->rng_state = WangHash(S->rng_state);
	return float(S->rng_state % 16777216)*(1.0f / 16777215.0f);
}

static bool molsOverlap(vec2 p0, float r0, vec2 p1, float r1) {
	float r = r0+r1;
	vec2 d = p0-p1;
	return dot(d,d) < (r*r);
}
static bool molsCollide(const Molecule& A, const Molecule& B, vec2 A_pos, vec2 B_pos, float A_rad, float B_rad) {
	return (A.body_idx != B.body_idx) && molsOverlap(A_pos, A_rad, B_pos, B_rad); // #OPT? are the type gets slow? cache in mol?
}
static bool siteIsNear(vec2 p_site, vec2 p_mol) {
	const float thresh = (BIND_POS_THRESH*MOL_RAD) * (BIND_POS_THRESH*MOL_RAD);
	const vec2 d = p_mol-p_site;
	return dot(d,d) < thresh;
}
static bool sitesBind(const BindSite& A, const BindSite& B) {
	return A.flavor == B.flavor && (
	      (A.gender == GENDER_BOTH && B.gender == GENDER_BOTH) ||
	      (A.gender == GENDER_MALE && B.gender == GENDER_FEMALE) ||
	      (A.gender == GENDER_FEMALE && B.gender == GENDER_MALE));
}
static bool bodyBreaks(u32 A_types_inside, u32 B_types_inside, const CatalysisRule* rules_start, const CatalysisRule* rules_end) {
	for (rules_start; rules_start != rules_end; rules_start++)
		if (A_types_inside == rules_start->a && B_types_inside == rules_start->b)
			return true;
	return false;
}
static bool molsBind(const SimState* S, const Molecule& A, const Molecule& B, pose wol_from_mol_A, pose wol_from_mol_B, int* a_site_idx, int* b_site_idx) {
	if (A.type != B.type) return false;
	const MoleculeType& TA = molTypeGet(S, A.type);
	const MoleculeType& TB = molTypeGet(S, B.type);
	for (int a = 0; a < TA.sites.count; ++a) {
		const BindSite& SA = TA.sites[a];
		if (siteIsFull(A,a)) continue;
		pose world_from_sa = wol_from_mol_A * SA.mol_from_sit;
		if (siteIsNear(world_from_sa.xy(),wol_from_mol_B.xy())) {
			for (int b = 0; b < TB.sites.count; ++b) {
				const BindSite& SB = TB.sites[b];
				if (siteIsFull(B,b)) continue;
				if (sitesBind(SA,SB)) {
					*a_site_idx = a;
					*b_site_idx = b;
					return true;
				}
			}
		}
	}
	return false;
}
static void collProject(vec2 p1, float m1, float r1, vec2 p2, float m2, float r2, vec2* proj1, vec2* proj2, vec2* n, float* inv_mass_sum_pass) {
    vec2 d = p2-p1;
	float len = length(d);
	vec2 n01 = len > 0.f ? d/len : vec2(0.f,1.f);
    vec2 n10 = -n01;
	float pen = (r1+r2) - len;
   	float inv_mass_sum = 1.f / (m1+m2);
    float amount = m2 * inv_mass_sum; 
	*proj1 = n10*pen*(amount);
    *proj2 = n01*pen*(1-amount);
	// pass on to next step
	*n = n10;
	*inv_mass_sum_pass = inv_mass_sum;
}

static void collResponse(float restitution, vec2 n, vec2 c, vec2 A_p, vec2 A_v, float A_w, float A_inv_m, float A_inv_moi, vec2 B_p, vec2 B_v, float B_w, float B_inv_m, float B_inv_moi, vec2* A_dv, float* A_dw, vec2* B_dv, float* B_dw) {	
	vec2 A_r = c - A_p;
	float A_r_cross_n = det(A_r,n);
	vec2 B_r = c - B_p;
	float B_r_cross_n = det(B_r,n);
	vec2 vr = (A_v + perp(A_r)*A_w) - (B_v + perp(B_r)*B_w); 
	float vrn = dot(vr,n);
	if (vrn < 0.f) {
		float j = (-vrn * (restitution + 1.f)) / (A_inv_m + B_inv_m + (sqr(A_r_cross_n)*A_inv_moi) + (sqr(B_r_cross_n)*B_inv_moi));
		*A_dv =  n * j * A_inv_m;
		*A_dw =  A_r_cross_n * j * A_inv_moi;
		*B_dv = -n * j * B_inv_m;
		*B_dw = -B_r_cross_n * j * B_inv_moi;
	} else {
		*A_dv = vec2(0.f);
		*A_dw = 0.f;
		*B_dv = vec2(0.f);
		*B_dw = 0.f;
	}
}

////////////////////////////

// hash
static Bunch<CollPair> coll_pairs;
static Bunch<int> cells;
static int cells_width = 0;
static int cells_height = 0;
static int cells_stride = 0;
static Bunch<bool> coll_hash;
static std::unordered_set<u64> coll_set;
static u32 coll_hash_mask = 0;
static Bunch<vec2> wol_from_mol_cache;
static Bunch<vec2i> coll_hash_pos_cache;
typedef Bin<int, 4> CheckIdxs;
static Bunch<CheckIdxs> coll_check_idxs;
static float world_from_cell_scale = 2.f;
static float cell_from_world_scale = 1.f / world_from_cell_scale;
static vec2 world_from_cell_pos = vec2(0.f);
static vec2 cell_from_world_pos = -world_from_cell_pos;

static void cells_clear() {
	memset(cells.ptr, 0xff, cells.bytes());
}

static bool in_range(int x, int y) {
	return x >= 0 && x < cells_width && y >= 0 && y < cells_height;
}

void setupCollHash(float bound_rad) {
	world_from_cell_scale = 2.f;
	cell_from_world_scale = 1.f / world_from_cell_scale;

	int b = int((bound_rad*2.f)*cell_from_world_scale); // collision hash's width/height is 2*radius
	cells_width = cells_height = b;
	cells_stride = b*CELL_COUNT_MAX;
	cells.setgarbage(b*cells_stride);
	cells_clear();

	world_from_cell_pos = (vec2(0.5f)-vec2(cells_width*0.5f,cells_height*0.5f))*world_from_cell_scale;
	cell_from_world_pos = - world_from_cell_pos;
}

/////////////////////////////////////

static void clear_cache(SimState* S, TimerTree* timer) {
	const Bunch<Molecule>& mols = S->molecules;

	timer->start("cache");
	wol_from_mol_cache.setgarbage(mols.count);
	for (int a = 0; a < mols.count; ++a) {	
		const Molecule& M = molGet(S, a);
		Body& B = bodyGet(S, M.body_idx);
		const vec2 wol_from_mol = WolFromMolPos(M,B);
		wol_from_mol_cache[a] = wol_from_mol;
		// stats
		B.types_inside |= (1<<M.type);
		B.mol_count += 1;
		if (B.some_mol_idx == NO_IDX) B.some_mol_idx = a; // if the body doesn't have 'some' molecule set, then set this one
	}
	timer->stop();
}
static void fill_hash_and_find_overalapped_pairs(SimState* S, TimerTree* timer) {
	const Bunch<Molecule>& mols = S->molecules;

	timer->start("compute pairs");

	timer->start("clear hash");
	coll_pairs.clear();
	coll_pairs.reserve(mols.count * 2);
	cells_clear();
	timer->stop();

	timer->start("write hash and pairs");
	
	//-- try doing raw pointer increment over all molecules?
	//-- whats the main cost here?? is it the overlap check? or the mem writes?
	//-- does in range check matter? it happens 4 times per mol...
	for (int a = 0; a < mols.count; ++a) {
		
#if 0
		
		Body& B = bodyGet(S, M.body_idx);
		const vec2 wol_from_mol = WolFromMolPos(M,B);
		wol_from_mol_cache[a] = wol_from_mol;
		// stats
		B.types_inside |= (1<<M.type);
		B.mol_count += 1;
		if (B.some_mol_idx == NO_IDX) B.some_mol_idx = a; // if the body doesn't have 'some' molecule set, then set this one
#else
		const vec2 wol_from_mol = wol_from_mol_cache[a];
#endif
		//#OPT can hoist this above and parallelize at least that...
		const vec2 cell_pos = (wol_from_mol + cell_from_world_pos) * cell_from_world_scale;
		const int x = (int)cell_pos.x;
		const int y = (int)cell_pos.y;
		const int idx = x*CELL_COUNT_MAX + y*cells_stride;
		const vec2i poss[] = {vec2i(x,y),vec2i(x+1,y),vec2i(x+1,y+1),vec2i(x,y+1)};
		const int idxs[] = {idx, idx+CELL_COUNT_MAX, idx+CELL_COUNT_MAX+cells_stride, idx+cells_stride};
		const Molecule& M = molGet(S, a);
		for (int i = 0; i < 4; ++i) {
			if (in_range(poss[i].x,poss[i].y)) {
				int* C = cells.ptr+idxs[i];
				int count = 0;
				while (count < CELL_COUNT_MAX) {
					const int b = *C;
					if (b >= 0) {
						const Molecule& B_mol = molGet(S, b);
						if (M.body_idx != B_mol.body_idx) {
							//u64 h = u64(a) | (u64(b) << 32);
							//if (coll_set.find(h) == coll_set.end()) {

								// add to hash regardless of whether it collides or not, to avoid rechecking either way
								//coll_set.insert(h);

								// if the object is in the collision cells, then it's wol_from_mol cache has already been updated and it's valid
								if (molsOverlap(wol_from_mol, molTypeGet(S, M.type).rad, wol_from_mol_cache[b], molTypeGet(S, B_mol.type).rad)) {
								//if (molsOverlap(wol_from_mol, rads[M.type], wol_from_mol_cache[b], rads[B_mol.type])) {	
									coll_pairs.push({a,b});
								}
							//}
						}

						++C;
						++count;
					} else {
						break;
					}
				}
				if (count < CELL_COUNT_MAX)
					*C = a;
			}
		}
	}
	timer->stop();
	
#if 1
	if (gui::Begin("Debug Hash Counts")) {
		timer->start("hash histo");
		int counts[CELL_COUNT_MAX+1];
		memset(counts, 0, sizeof(counts));
		for (int i = 0; i < cells.count/CELL_COUNT_MAX; ++i) {
			int c = 0;
			for (int k = 0; k < CELL_COUNT_MAX; ++k) {
				if (cells[i+k] >= 0)
					++c;
				else
					break;
			}
			counts[c] += 1;
		}
		timer->stop();
	
		int total = 0;
		for (int i = 0; i < CELL_COUNT_MAX+1; ++i)
			total += counts[i];
		for (int i = 0; i < CELL_COUNT_MAX+1; ++i) 
			gui::Text("%d: %d %6.2f%%", i, counts[i], counts[i]/float(total));
		gui::Text("total %d", total);
	} gui::End();
#endif

#if 0
	for (int a = 0; a < mols.count; ++a) {
		if (mols[a].rad == 0.f) continue;
		for (int b = a+1; b < mols.count; ++b) {
			if (mols[b].rad == 0.f) continue;
			if (molCollide(mols[a],mols[b])) {
				coll_pairs.push({a,b});
			}
		}
	}
#else

	
	//log("coll pairs %d\n", coll_pairs.count);
#if 0
	
	timer->start("check");
	int same = 0;
	for (int a = 0; a < coll_pairs.count; ++a) {
		for (int b = a+1; b < coll_pairs.count; ++b) {
			if (coll_pairs[a] == coll_pairs[b])
				same += 1;
		}
	}
	timer->stop();
	log("%d same\n", same);
#endif

#endif
	timer->stop();
}
static void collide(SimState* S, TimerTree* timer) {
	timer->start("collide");

	const CatalysisRule* catalysis_rules_start = S->catalysis_rules.ptr;
	const CatalysisRule* catalysis_rules_end = S->catalysis_rules.ptr + S->catalysis_rules.count;

	for (int i = 0; i < coll_pairs.count; ++i) {
		CollPair P = coll_pairs[i];

		Molecule* A_mol = &molGet(S, P.a);
		Molecule* B_mol = &molGet(S, P.b);
		Body* A_bod = &bodyGet(S, A_mol->body_idx);
		Body* B_bod = &bodyGet(S, B_mol->body_idx);
		pose wol_from_mol_a = WolFromMol(*A_mol,*A_bod);
		pose wol_from_mol_b = WolFromMol(*B_mol,*B_bod);
		float A_rad = molTypeGet(S, A_mol->type).rad;
		float B_rad = molTypeGet(S, B_mol->type).rad;

		// do mols still overlap? they may have moved from earlier collisions in the loop
		// position may have changed (due to projection)
		// body index may have changed (due to merge or split)
		if (!molsCollide(*A_mol,*B_mol, wol_from_mol_a.xy(), wol_from_mol_b.xy(), A_rad, B_rad)) continue; //#OPT should pass radii in here, since we've already accessed them above

		vec2 projA,projB;
		vec2 coll_n;
		float inv_mass_sum_pass;
		collProject(wol_from_mol_a.xy(), A_bod->mass, A_rad,
					wol_from_mol_b.xy(), B_bod->mass, B_rad,
					&projA,&projB,&coll_n,&inv_mass_sum_pass);

		A_bod->wol_from_bod.xy(A_bod->wol_from_bod.xy() + projA);
		B_bod->wol_from_bod.xy(B_bod->wol_from_bod.xy() + projB);

		// update molecule world poses post collision projection
		wol_from_mol_a = WolFromMol(*A_mol,*A_bod);
		wol_from_mol_b = WolFromMol(*B_mol,*B_bod);

		// compute true collision point (collision normal is still correct because projection is along it)
		vec2 coll_c = (wol_from_mol_a.xy()+wol_from_mol_b.xy())*0.5f;

		//if (P.a == edit.inspect_mol_idx || P.b == edit.inspect_mol_idx) {
		//	int break_here = 0;
		//}

		SitePair sites[] = { {-1,-1}, {-1,-1} };
		int site_choice = 0;
		bool ab_bond = molsBind(S, *A_mol,*B_mol,wol_from_mol_a,wol_from_mol_b,&sites[0].a,&sites[0].b);
		bool ba_bond = molsBind(S, *B_mol,*A_mol,wol_from_mol_b,wol_from_mol_a,&sites[1].a,&sites[1].b);

		// don't allow smaller molecules to bind larger ones
		if (ab_bond && (A_bod->mass < B_bod->mass)) ab_bond = false;
		if (ba_bond && (B_bod->mass < A_bod->mass)) ba_bond = false;

		// don't allow to exceed maximum bond size
		if (int(A_bod->mol_count + B_bod->mol_count) > molTypeGet(S, A_mol->type).bond_max_count) ab_bond = ba_bond =  false; // #DODGY until we have a more principled way of dealing with max ring size?
		
		// don't allow "large" structures to bond together, only small-small or small-large
		bool A_is_large = A_bod->mol_count >= 3;
		bool B_is_large = B_bod->mol_count >= 3;
		if (A_is_large && B_is_large) ab_bond = ba_bond = false;

		if (ab_bond || ba_bond) {
			// want to bond smaller object to the larger one, so that the smaller object "pops in" rather than the large one
			bool do_swap = false;
			if (ab_bond && ba_bond) {
				if (B_bod->mass > A_bod->mass)
					do_swap = true;
			} else {
				if (ba_bond)
					do_swap = true;
			}
			if (do_swap) {
				swap(A_mol, B_mol);
				swap(A_bod, B_bod);
				swap(wol_from_mol_a, wol_from_mol_b);
				swap(A_rad, B_rad);
				swap(P.a,P.b);
				site_choice = 1;
			}

			int A_site_idx = sites[site_choice].a;
			int B_site_idx = sites[site_choice].b;

			/// CAREFUL! Make sure you only use the above swapped variables

			pose bod_a_from_wol = ~A_bod->wol_from_bod;
			bod_a_from_wol.zw(normalize(bod_a_from_wol.zw()));

			// take the sites
			fillSite(*A_mol, A_site_idx, P.b);
			fillSite(*B_mol, B_site_idx, P.a);

			//float rot_moment_before = computeSumRotMomentum();
			float sum_lin_kin_energy_pre = computeLinearKineticEnergy(A_bod) + computeLinearKineticEnergy(B_bod);
			float sum_rot_kin_energy_pre = computeRotationalKineticEnergy(A_bod) + computeRotationalKineticEnergy(B_bod);

			// compute response
			vec2 A_dv,B_dv;
			float A_dw,B_dw;
			vec2 A_com = A_bod->wol_from_bod * A_bod->bod_from_com;
			vec2 B_com = B_bod->wol_from_bod * B_bod->bod_from_com;
			collResponse(0.f, coll_n, coll_c,
			             A_com, A_bod->dp, A_bod->dr, 1.f/A_bod->mass, 1.f/A_bod->moi, //#OPT memoize inv mass properties
			             B_com, B_bod->dp, B_bod->dr, 1.f/B_bod->mass, 1.f/B_bod->moi,
			             &A_dv,&A_dw, &B_dv,&B_dw);

			// align first, so that mass transfer has the right centers of mass
			pose wol_from_sit = wol_from_mol_a * molTypeGet(S, A_mol->type).sites[A_site_idx].mol_from_sit;
			pose wol_from_wol = wol_from_sit * ~wol_from_mol_b;
			B_bod->wol_from_bod = wol_from_wol * B_bod->wol_from_bod;
			B_bod->wol_from_bod.zw(normalize(B_bod->wol_from_bod.zw())); // #HEAVYHANDED #OPT
			/// CAREFUL! world_from_mol_b is now invalid!
			// #DEBUG catch usage
			wol_from_mol_b = vec4(FLT_MAX); 

			// transfer mass of B into A
			float A_mass = A_bod->mass;
			float B_mass = B_bod->mass;
			float A_moi = A_bod->moi;
			float B_moi = B_bod->moi;
			vec2 wol_from_com_b = B_bod->wol_from_bod * B_bod->bod_from_com;
			vec2 bod_a_from_com_b = bod_a_from_wol * wol_from_com_b;
			bodyAddMass(*A_bod, bod_a_from_com_b, B_bod->mass, B_bod->moi, B_bod->mol_count);

			// transfer velocity
			A_bod->dp = ((A_bod->dp + A_dv)*A_mass + (B_bod->dp + B_dv)*B_mass) / (A_bod->mass); // A_bod->mass is now the sum of the two, from above calc
			A_bod->dr = ((A_bod->dr + A_dw)*A_moi + (B_bod->dr + B_dw)*B_moi) / (A_bod->moi); // #WRONG ?
			assert(!isnan(A_bod->dp) && !isnan(A_bod->dr));

			float sum_lin_kin_energy_post = computeLinearKineticEnergy(A_bod);
			float lin_kin_energy_lost = sum_lin_kin_energy_pre - sum_lin_kin_energy_post;
			S->dissapated_lin_kinetic_energy += lin_kin_energy_lost;
			float sum_rot_kin_energy_post = computeRotationalKineticEnergy(A_bod);
			float rot_kin_energy_lost = sum_rot_kin_energy_pre - sum_rot_kin_energy_post;
			S->dissapated_rot_kinetic_energy += rot_kin_energy_lost;

			// capture lost energy into vibration, as well as the old body's vibration
			A_bod->vibration += B_bod->vibration + // old body's vibration energy
			                    (lin_kin_energy_lost + rot_kin_energy_lost);
			//assert((lin_kin_energy_lost + rot_kin_energy_lost) >= 0.f);

			//float rot_moment_after = computeSumRotMomentum();
			//log("before %f after %f (%f/%d)\n", rot_moment_before, rot_moment_after, rot_moment_after-rot_moment_before,rot_moment_before == rot_moment_after);

			// transfer all mols from B to A
			int B_bod_idx = B_mol->body_idx; // store this, since it's about to be overwritten

			int mol_idx = P.b;
			int B_edge_idx = mol_idx;
			while (mol_idx != NO_IDX) {
				Molecule& M = molGet(S, mol_idx);
				pose wol_from_mol = WolFromMol(M, *B_bod); // we know that the B_body is the same
				M.bod_from_mol = bod_a_from_wol * wol_from_mol; // recompute pose in A frame
				M.bod_from_mol.zw(normalize(M.bod_from_mol.zw())); // #HEAVYHANDED #OPT
				M.body_idx = A_mol->body_idx;
				B_edge_idx = mol_idx;
				mol_idx = getBond(M, 1-B_site_idx); // walk in the opposite direction of A
			}

			if ((int)A_bod->mol_count == molTypeGet(S, A_mol->type).bond_max_count) {
				// if we have completed a ring, we have to join the other side as well
				// from the above loop, we found the edgemost B molecule, so now do the same for A
				int A_edge_idx = P.a;
				int m_idx = P.a;
				while (m_idx != NO_IDX) {
					const Molecule& M = molGet(S, m_idx);
					A_edge_idx = m_idx;
					m_idx = getBond(M, 1-A_site_idx);
				}
				Molecule& EA_mol = molGet(S, A_edge_idx);
				Molecule& EB_mol = molGet(S, B_edge_idx);
				fillSite(EA_mol, 1-A_site_idx, B_edge_idx);
				fillSite(EB_mol, 1-B_site_idx, A_edge_idx);
			}

			// erase body B
			bodyRem(S, B_bod_idx);
		} else {
			vec2 A_dv,B_dv;
			float A_dw,B_dw;
			vec2 A_com = A_bod->wol_from_bod * A_bod->bod_from_com;
			vec2 B_com = B_bod->wol_from_bod * B_bod->bod_from_com;
			//float A_cor = A_bod->types_inside == 1 ? molTypeGet(S, A_mol->type).cor : 1.f;
			//float B_cor = B_bod->types_inside == 1 ? molTypeGet(S, B_mol->type).cor : 1.f;
			// commented out the coeficient of restitution math A_cor*B_cor
			//float rot_moment_before = computeSumRotMomentum();
			//assert(A_bod->vibration >= 0.f && B_bod->vibration >= 0.f);
#define DO_VIBRATION 1	
#if DO_VIBRATION==1
			float vibration_total = A_bod->vibration + B_bod->vibration;// + 50.f;
			float sum_lin_kin_energy_pre = computeLinearKineticEnergy(A_bod) + computeLinearKineticEnergy(B_bod);
			float sum_rot_kin_energy_pre = computeRotationalKineticEnergy(A_bod) + computeRotationalKineticEnergy(B_bod);
			float sum_energy_pre = sum_lin_kin_energy_pre+sum_rot_kin_energy_pre;
			float desired_sum_energy_post = max(0.f, vibration_total + sum_energy_pre);
			float cor = sqrt(desired_sum_energy_post / sum_energy_pre); // cor is equal to the square root of the ratio of post and pre energies (or the ratio of post and pre speed)
#else
			float cor = 1.0f;
#endif
			vec2 A_vel_pre = A_bod->dp;
			vec2 B_vel_pre = B_bod->dp;
			collResponse(cor, coll_n, coll_c,
			             A_com, A_bod->dp, A_bod->dr, 1.f/A_bod->mass, 1.f/A_bod->moi, //#OPT memoize inv mass properties
			             B_com, B_bod->dp, B_bod->dr, 1.f/B_bod->mass, 1.f/B_bod->moi,
			             &A_dv,&A_dw, &B_dv,&B_dw);
			A_bod->dp += A_dv;
			A_bod->dr += A_dw;
			B_bod->dp += B_dv;
			B_bod->dr += B_dw;
			assert(!isnan(A_bod->dp) && !isnan(A_bod->dr) && !isnan(B_bod->dp) && !isnan(B_bod->dr));

#if DO_VIBRATION==1
			//float cor_vel_a = length(A_bod->dp) / length(A_vel_pre);
			//float cor_vel_b = length(B_bod->dp) / length(B_vel_pre);
			float sum_lin_kin_energy_post = computeLinearKineticEnergy(A_bod) + computeLinearKineticEnergy(B_bod);
			float sum_rot_kin_energy_post = computeRotationalKineticEnergy(A_bod) + computeRotationalKineticEnergy(B_bod);
			float sum_energy_post = sum_lin_kin_energy_post+sum_rot_kin_energy_post;
			float desired_diff = desired_sum_energy_post - sum_energy_post;
			//float final_diff = sum_energy_post - sum_energy_pre;
			A_bod->vibration = desired_diff*0.5f;
			B_bod->vibration = desired_diff*0.5f;
#endif

			//float rot_moment_after = computeSumRotMomentum();
			//log("before %f after %f (%f/%d)\n", rot_moment_before, rot_moment_after, rot_moment_after-rot_moment_before,rot_moment_before == rot_moment_after);

			int mol_idxs[] = {P.a,P.b};
			Molecule* break_mols[] = {A_mol,B_mol};
			Body* break_bodies[] = {A_bod,B_bod};
			int bod_idxs[] = {A_mol->body_idx,B_mol->body_idx};
			float impact_str[] = {lengthsqr(A_dv) + A_dw*A_dw, lengthsqr(B_dv) + B_dw*B_dw}; 
			for (int s = 0; s < 2; ++s) {
				if (molTypeGet(S, break_mols[s]->type).bond_max_count > 0 && break_bodies[s]->types_inside == ((u32)1<<(u32)break_mols[s]->type)) { // if it's a shell forming molecule and that's all the body contains
					bool breaks = impact_str[s]> S->capsid_break_threshold;
					breaks &= break_bodies[s]->mol_count > 1; // bodies of 1 molecule can't break
					if (breaks) {
						const int arity = 2;
						int new_body_idxs[arity];
						Body* new_bodies[arity];
						int bond_idxs[arity];
						int mol_counts[arity];
						for (int a = 0; a < arity; ++a) {
							new_body_idxs[a] = bodyAdd(S);
							new_bodies[a] = &bodyGet(S, new_body_idxs[a]);
							new_bodies[a]->wol_from_bod = break_bodies[s]->wol_from_bod; // place new bodies wherever the old one was (if we change this, change bod_from_mol below
							bond_idxs[a] = getBond(*break_mols[s],a); // start from neighbors of hit molecule
							mol_counts[a] = 0;
						}

						while (true) {
							bool no_advance = true;
							for (int a = 0; a < arity; ++a) { // advance all bonds
								if (bond_idxs[a] != NO_IDX) { // if not end of line...
									Molecule& N_mol = molGet(S, bond_idxs[a]);
									if (N_mol.body_idx == break_mols[s]->body_idx) { // .. and if this molecule still hasn't been claimed by a new body
										const MoleculeType& T = molTypeGet(S, N_mol.type);
									
										// update the body
										bodyAddMass(*new_bodies[a], N_mol.bod_from_mol.xy(), T.mass, T.moi, 1);
										new_bodies[a]->types_inside |= (1<<N_mol.type);
										new_bodies[a]->type = simTypeOfCompound(S, new_bodies[a]->types_inside);
										if (new_bodies[a]->some_mol_idx == NO_IDX) new_bodies[a]->some_mol_idx = bond_idxs[a]; // if body doesn't have 'some mol' set, then set this one

										// update the molecule
										N_mol.body_idx = new_body_idxs[a];
										//bod_from_mol remains same, because the new body is where the old one was

										bond_idxs[a] = getBond(N_mol, a); //#DODGY this continues walking in same direction, which is fine, but not for arbitrary graphs
										mol_counts[a] += 1;
										no_advance = false;
									} else { // we hit a molecule which is already taken, so break that bond (but in the opposite direction, because we've gone "into it" already
										emptySite(N_mol, 1-a);
									}
								}
							}
							if (no_advance)
								break;
						}

						// pick bond branch with the least molecules in it..
						int bond_with_least_mols = -1;
						int min_count = break_bodies[s]->mol_count;
						for (int a = 0; a < arity; ++a) {
							if (mol_counts[a] < min_count) {
								min_count = mol_counts[a];
								bond_with_least_mols = a;
							}
						}

						// ..and split the opposite bond (in other words, join the hit molecule to the piece which has the fewest molecules in it (including potentially 0, in a 2 piece case))
						{
							{
								// #CLEANUP this is a copy pasta of the above loop...
								int a = bond_with_least_mols;
								Molecule& N_mol = *break_mols[s];
								const MoleculeType& T = molTypeGet(S, N_mol.type);	

								// update the body
								bodyAddMass(*new_bodies[a], N_mol.bod_from_mol.xy(), T.mass, T.moi, 1);
								new_bodies[a]->types_inside |= (1<<N_mol.type);
								new_bodies[a]->type = simTypeOfCompound(S, new_bodies[a]->types_inside);
								if (new_bodies[a]->some_mol_idx == NO_IDX) new_bodies[a]->some_mol_idx = mol_idxs[s]; // if body doesn't have 'some mol' set, then set this one

								// update the molecule
								N_mol.body_idx = new_body_idxs[a];
								// #CLEANUP end

								// break the actual bond
								assert(arity == 2);
								int bond = getBond(*break_mols[s], 1-a);
								assert(bond != NO_IDX); // we are more than 1 molecule, and we should be adding ourself to the body of least molecules, even if that body has 0, so this must be a valid bond
								Molecule& O_mol = molGet(S, bond);
								emptySite(*break_mols[s], 1-a); // #DODGY only works with arity 2 :(
								emptySite(O_mol, a); // #DODGY again, assumes the other molecule's arity is also 2 and it faces the other way :(							
							}
							{ // compute kinematics of the new bodies
								vec2 break_body_com = break_bodies[s]->wol_from_bod * break_bodies[s]->bod_from_com;
								for (int a = 0; a < arity; ++a) {
									vec2 new_body_com = new_bodies[a]->wol_from_bod * new_bodies[a]->bod_from_com;
									vec2 r = new_body_com - break_body_com;
									new_bodies[a]->dp = break_bodies[s]->dp + perp(r)*break_bodies[s]->dr;
									new_bodies[a]->dr = break_bodies[s]->dr;
									new_bodies[a]->vibration = break_bodies[s]->vibration/arity; // each body gets a percentage of old body's vibration energy
									assert(!isnan(new_bodies[a]->dp) && !isnan(new_bodies[a]->dr));
								}
							}
						}

						// lastly, delete the old body
						bodyRem(S, bod_idxs[s]);
					}
				} else { // catalysis
					//#OPT can early out here and not check all rules, if we just knew if either type was a catalyst at all?
					bool ab_break = bodyBreaks(A_bod->types_inside, B_bod->types_inside, catalysis_rules_start, catalysis_rules_end);
					bool ba_break = bodyBreaks(B_bod->types_inside, A_bod->types_inside, catalysis_rules_start, catalysis_rules_end);
					float chance = S->catalysis_chances[ab_break ? A_bod->type : B_bod->type];
					if ((ab_break || ba_break) && (simRandFloat(S) <= chance)) {
						int break_bod_idx = ab_break ? B_mol->body_idx : A_mol->body_idx;
						Body* break_body = ab_break ? B_bod : A_bod;
						Molecule* break_mol = ab_break ? B_mol : A_mol;
						int break_mol_idx = ab_break ? P.b : P.a;
						vec2 break_body_com = break_body->wol_from_bod * break_body->bod_from_com;

						//#NONGEN this assumes all catalysts are composed of 2 pieces..
						const int num_mols = 2;
						assert(break_body->mol_count == 2);
						int proc_idxs[num_mols] = { -1, -1 };
						proc_idxs[0] = break_mol_idx;
						proc_idxs[1] = getBond(*break_mol, 0);
						if (proc_idxs[1] == NO_IDX)
							proc_idxs[1] = getBond(*break_mol, 1);

						for (int p = 0; p < num_mols; ++p) {
							int m = proc_idxs[p];
							Molecule& M = molGet(S, m);
							if (M.body_idx == break_bod_idx) {
								const MoleculeType& T = molTypeGet(S, M.type);
								pose wol_from_mol = WolFromMol(M, *break_body);
								M.bod_from_mol = identity();
								int new_body_idx = bodyAdd(S);
								Body& C = bodyGet(S, new_body_idx);
								C.wol_from_bod = wol_from_mol;
								vec2 r = wol_from_mol.xy() - break_body_com;
								C.dp = break_body->dp + perp(r)*break_body->dr;
								C.dr = break_body->dr;
								C.vibration = break_body->vibration / num_mols;
								assert(!isnan(C.dp) && !isnan(C.dr));
								C.types_inside |= (1<<M.type);
								C.type = simTypeOfCompound(S, C.types_inside);
								C.some_mol_idx = m;
								bodyAddMass(C, vec2(0.f), T.mass, T.moi, 1);
								M.body_idx = new_body_idx;
								emptySite(M, 0); //#NONGEN 
								emptySite(M, 1);
							} else {
								assert(false);
							}
						}
						bodyRem(S, break_bod_idx);
						break; // #DODGY just to get out of the "check both sides loop" :(
					}
				}
			}
		}

		//mols_vis[P.a].is_colliding = true;
		//mols_vis[P.b].is_colliding = true;
	}
	timer->stop();
}
static void collide_with_bounds(SimState* S, TimerTree* timer) {
	//#WRONG #DODGY #OPT using the wol_from_mol cache is a huge speedup, because we not only avoid the compute but we avoid the body get. but the cache is stale here, so needs fixing...but speedup is significant. 2.5ms down to 0.3ms
	timer->start("bound collision");
	for (int i = 0; i < S->molecules.count; ++i) {
		Molecule* M = S->molecules.ptr + i;
		const float rad = molTypeGet(S, M->type).rad;
		vec2 d = /* vec2(0.f) */ - wol_from_mol_cache[i];
		//vec2 d = - WolFromMolPos(*M,bodyGet(sim_state_ptr, M->body_idx)); // unoptimized, but more correct
		float len = length(d);
		if ((len+rad) > S->bound_rad) {
			Body& B = bodyGet(S, M->body_idx);
			vec2 nd = d/len;
			B.wol_from_bod.xy(B.wol_from_bod.xy() + nd * (len+rad-S->bound_rad));
			vec2 reflect = nd * (dot(-nd,B.dp)*2.f);
			B.dp += reflect;
			S->lin_bound_momentum += -reflect*B.mass;
			S->rot_bound_momentum += 0.f;
		}
	}
	timer->stop();
}
float computeSumRotMomentum(const SimState* S) {
	float rot_momentum_sum = 0.f;
	for (int i = 0; i < S->bodies.count; ++i) {
		if (!S->bodies[i].valid) continue;
		const Body& B = S->bodies[i];
		vec2 wol_from_com = B.wol_from_bod * B.bod_from_com;
		rot_momentum_sum += det(wol_from_com, B.dp*B.mass) + B.moi*B.dr;
	}
	return rot_momentum_sum;
}

#include "core/draw_prim.h"
static void tallyShellContents(SimState* S, SimStats* stats, TimerTree* timer) {
	timer->start("tally shell contents");
	(void)stats;

	const Bunch<Molecule>& mols = S->molecules;

	timer->start("clear hash");
	cells_clear();
	timer->stop();

	timer->start("write hash");
	for (int a = 0; a < mols.count; ++a) {
		const Molecule& M = molGet(S, a);
		const Body& B = bodyGet(S, M.body_idx);
		const vec2 wol_from_mol = WolFromMolPos(M,B);

		vec2 cell_pos = (wol_from_mol + cell_from_world_pos) * cell_from_world_scale;
		int x = (int)cell_pos.x;
		int y = (int)cell_pos.y;
		int idx = x*CELL_COUNT_MAX + y*cells_stride;
		vec2i poss[] = {vec2i(x,y),vec2i(x+1,y),vec2i(x+1,y+1),vec2i(x,y+1)};
		int idxs[] = {idx, idx+CELL_COUNT_MAX, idx+CELL_COUNT_MAX+cells_stride, idx+cells_stride};
		
		for (int i = 0; i < 4; ++i) {
			if (in_range(poss[i].x,poss[i].y)) {
				int* C = cells.ptr+idxs[i];
				int count = 0;
				while (count < CELL_COUNT_MAX) {
					const int b = *C;
					if (b >= 0) {
						const Molecule& B_mol = molGet(S, b);
						if (M.body_idx == B_mol.body_idx) {
							count = -1; // don't write same body again, since we only care about counting bodies
							break;
						}
						++C;
						++count;
					} else {
						break;
					}
				}
				if (count >= 0 && count < CELL_COUNT_MAX)
					*C = a;
			}
		}
	}
	timer->stop();

	timer->start("count");
	stats->shell_guts.clear();

	const Bunch<Body>& bods = S->bodies;
	u8 first_capsid_type = 0;
	for (u8 t = 0; t < typeGetCount(S); ++t) {
		if (molTypeGet(S, t).bond_max_count > 0) {
			first_capsid_type = t;
			break;
		}
	}

	for (int idx = 0; idx < bods.count; ++idx) {
		if (!bods[idx].valid) continue;
		const Body& B = bods[idx];
		vec2 p = B.wol_from_bod * B.bod_from_com;

		const Molecule& M = molGet(S, B.some_mol_idx);
		const MoleculeType& T_m = molTypeGet(S, M.type);
		if (B.types_inside == (u32)(1<<first_capsid_type) && B.mol_count == (T_m.bond_max_count)) {
			int mol_idx_l = B.some_mol_idx;
			int mol_idx_r = getBond(M, 0);
			if (mol_idx_r == NO_IDX) {
				mol_idx_r = mol_idx_l;
				mol_idx_l = getBond(M, 1);
			}
			assert(mol_idx_l != NO_IDX && mol_idx_r != NO_IDX);
			pose wol_from_mol_l = WolFromMol(molGet(S, mol_idx_l), B);
			pose wol_from_mol_r = WolFromMol(molGet(S, mol_idx_r), B);
			vec2 a = wol_from_mol_l.xy();
			vec2 b = wol_from_mol_r.xy();
			float t = TAU/T_m.bond_max_count;
			vec2 a_to_b = b-a;
			float L = length(a_to_b);
			float D = (L*0.5f) / tan(t*0.5f);
			vec2 c = (b+a)*0.5f + normalize(perp(a_to_b))*D;
			float r = sqrt(sqr(L*0.5f) + sqr(D)) - T_m.rad*1.f; // inset by 1 radii, to make search bound tighter (testing collision with bound not total overlap)

			//primColor(vec4(simBodyCol(S, B.type),1.f)*0.2f);
			//primCircle(c, r);
			
			vec2i cell_ll = vec2i(((c-vec2(r)) + cell_from_world_pos) * cell_from_world_scale);
			vec2i cell_ur = vec2i(((c+vec2(r)) + cell_from_world_pos) * cell_from_world_scale) + vec2i(1);
			static Bunch<int> bod_idx_inside;
			static Bunch<int> mol_idx_inside;
			bod_idx_inside.clear();
			mol_idx_inside.clear();
			for (int y = cell_ll.y; y <= cell_ur.y; ++y) {
				for (int x = cell_ll.x; x <= cell_ur.x; ++x) {
					int i = x*CELL_COUNT_MAX + y*cells_stride;
					for (int k = 0; k < CELL_COUNT_MAX; ++k) {
						int mol_idx = cells[i+k];
						if (mol_idx >= 0) {
							const Molecule& M_a = molGet(S, mol_idx);
							if (M_a.body_idx != idx && !bod_idx_inside.find(M_a.body_idx)) {
								const Body& B_a = bodyGet(S, M_a.body_idx);
								vec2 wol_from_mol_pos = WolFromMolPos(M_a, B_a);
								if (molsOverlap(c,r, wol_from_mol_pos, molTypeGet(S, M_a.type).rad)) {
									bod_idx_inside.push(M_a.body_idx);
									mol_idx_inside.push(mol_idx);
								}
							}
						}
						else
							break;
					}
					//primColor(vec4(0.f,0.f,0.f,0.2f));
					//primBox(world_from_cell_pos.x+x*world_from_cell_scale,world_from_cell_pos.y+y*world_from_cell_scale, 0.9f,0.9f);
				}
			}
			//int mol_count = 0;
			BodyHist& hist = stats->shell_guts.push();
			hist.body_idx = idx;
			hist.body_type_masses.pushi(0, typeGetCount(S));
			for (int i = 0; i < bod_idx_inside.count; ++i) {
				const Body& B_a = bodyGet(S, bod_idx_inside[i]);
				//mol_count += B_a.mol_count;
				hist.body_type_masses[B_a.type] += B_a.mass;
			}
		}
	}

	timer->stop();
	
	timer->stop();
}

void simFinalize(SimState* S) {
	// add enough free bodies such that we don't need to create any during runtime
	int extra_bodies = int(S->molecules.count - S->bodies.count + 1); // enough bodies to cover every single molecule, +1 more since we create bodies during splitting before we destroy the existing one
	for (int i = 0; i < extra_bodies; ++i)
		bodyAddFree(S);
	setupCollHash(S->bound_rad);
}
 
static void integrate_and_clear(SimState* S, float dt, TimerTree* timer) {
	// update positions
	timer->start("integrate and clear");
	
	/*
	if (S->heat_pump.enable) {
		float avg_bod_energy = 0.f;
		float mass_sum = 0.f;
		for (Body* B = S->bodies.ptr; B != (S->bodies.ptr+S->bodies.count); B++) {
			if (B->valid) {
				avg_bod_energy += computeLinearKineticEnergy(B) + computeRotationalKineticEnergy(B) + B->vibration;
				mass_sum += B->mass;
			}
		}
		if (mass_sum > 0.f)
			avg_bod_energy /= mass_sum;
	}
	*/

	float target_avg_energy = simHeatPumpAvgEnergy(S);

	for (Body* B = S->bodies.ptr; B != (S->bodies.ptr+S->bodies.count); B++) {
		if (B->valid)  {
			if (S->heat_pump.enable) { // pump energy into molecules
				float e = computeLinearKineticEnergy(B) + computeRotationalKineticEnergy(B) + B->vibration;
				B->vibration += (target_avg_energy - e) * S->heat_pump.permeability;
			}
			// integrate
			// #OPT if we didn't need to ofset by com, then we save another 10% or so. but the entire loop is only around 1ms so it's not too worth it now
			pose bod_from_com = vec4(B->bod_from_com,1.f,0.f);
			pose wol_from_com = B->wol_from_bod * bod_from_com;
			wol_from_com.xy(wol_from_com.xy() + B->dp*dt);
			// #HEAVYHANDED #OPT renormalization is about 10% cost, and the rotate is another 30%
			wol_from_com.zw(normalize(rotate(quat2d(B->dr*dt), wol_from_com.zw()))); 
			B->wol_from_bod = wol_from_com * ~bod_from_com;

			//clear
			B->types_inside = 0;
			B->mol_count = 0;
		}
	}
	timer->stop();
}

void simUpdate(SimState* S, float dt, int substeps, TimerTree* timer, StatsReq sreq) {
	(void)sreq;

	dt = dt / substeps;
	for (int i = 0; i < substeps; ++i) {
		timer->start("substep");
		integrate_and_clear(S, dt, timer);
		clear_cache(S, timer);
		fill_hash_and_find_overalapped_pairs(S, timer);
		collide(S, timer);
		collide_with_bounds(S, timer);
		S->time += dt;
		timer->stop();
	}
	S->step_count += 1;
	//if (sreq.shell_contents) tallyShellContents(S, sreq.stats, timer);
}

void simStatsOnly(SimState* S, TimerTree* timer, StatsReq sreq) {
	if (sreq.shell_contents) tallyShellContents(S, sreq.stats, timer);
}

////////////////////////////

int compoundAdd(SimState* S, const BodyTemplate& T, pose wol_from_bod) {
	bodyAddFree(S);
	int body_idx = bodyAdd(S);
	Body& B = bodyGet(S, body_idx);
	B.wol_from_bod = wol_from_bod;
	int prev_idx = NO_IDX;
	for (int i = 0; i < T.mols.count; ++i) {
		Molecule M = T.mols[i];
		M.body_idx = body_idx;
		setBond(M, 1, prev_idx);
		S->molecules.push(M);
		int curr_idx = (int)S->molecules.count-1;
		if (prev_idx != NO_IDX)
			setBond(S->molecules[prev_idx], 0, curr_idx); 
		prev_idx = curr_idx;
		bodyAddMass(B, M.bod_from_mol.xy(), molTypeGet(S, M.type).mass, molTypeGet(S, M.type).moi, 1);
		B.types_inside |= (1<<M.type);
		B.type = simTypeOfCompound(S, B.types_inside);
		if (B.some_mol_idx == NO_IDX) B.some_mol_idx = curr_idx;
	}
	return body_idx;
}
int shellAdd(SimState* S, const BodyTemplate& T, pose wol_from_bod) {
	bodyAddFree(S);
	int body_idx = bodyAdd(S);
	Body& B = bodyGet(S, body_idx);
	B.wol_from_bod = wol_from_bod;
	int prev_idx = NO_IDX;
	int first_idx = NO_IDX;
	for (int i = 0; i < T.mols.count; ++i) {
		Molecule M = T.mols[i];
		M.body_idx = body_idx;
		S->molecules.push(M);
		int curr_idx = (int)S->molecules.count-1;
		if (first_idx == NO_IDX)
			first_idx = curr_idx;
		if (prev_idx != NO_IDX) {
			s32 A_mol_idx = prev_idx;
			s32 B_mol_idx = curr_idx;
			Molecule* A_mol = &S->molecules[A_mol_idx];
			Molecule* B_mol = &S->molecules[B_mol_idx];
			pose wol_from_mol_a = WolFromMol(*A_mol, B);
			pose wol_from_mol_b = WolFromMol(*B_mol, B);
			SitePair site ={-1,-1};
			bool ab_bond = molsBind(S, *A_mol,*B_mol,wol_from_mol_a,wol_from_mol_b,&site.a,&site.b);
			assert(ab_bond);
			// take the sites
			fillSite(*A_mol, site.a, B_mol_idx);
			fillSite(*B_mol, site.b, A_mol_idx);
		}
		if (i == T.mols.count-1) {
			s32 A_mol_idx = curr_idx;
			s32 B_mol_idx = first_idx;
			Molecule* A_mol = &S->molecules[A_mol_idx];
			Molecule* B_mol = &S->molecules[B_mol_idx];
			pose wol_from_mol_a = WolFromMol(*A_mol, B);
			pose wol_from_mol_b = WolFromMol(*B_mol, B);
			SitePair site ={-1,-1};
			bool ab_bond = molsBind(S, *A_mol,*B_mol,wol_from_mol_a,wol_from_mol_b,&site.a,&site.b);
			assert(ab_bond);
			// take the sites
			fillSite(*A_mol, site.a, B_mol_idx);
			fillSite(*B_mol, site.b, A_mol_idx);
		}
		
		prev_idx = curr_idx;
		bodyAddMass(B, M.bod_from_mol.xy(), molTypeGet(S, M.type).mass, molTypeGet(S, M.type).moi, 1);
		B.types_inside |= (1<<M.type);
		B.type = simTypeOfCompound(S, B.types_inside);
		if (B.some_mol_idx == NO_IDX) B.some_mol_idx = curr_idx;
	}
	return body_idx;
}


////////////////////////////////

void drawCells() {
#if 0
	vec2 bound_dims = vec2(cells_width);
	primColor(0,0,0);
	primPush(trans(world_from_cell_pos, 0.f, world_from_cell_scale) * trans(vec2(-0.5f),0.f,1.f));
	primLine(vec2(0.f), vec2(bound_dims.x,0.f), 0.1f);
	primLine(vec2(0.f), vec2(0.f,bound_dims.y), 0.1f);
	primLine(vec2(0.f,bound_dims.y), vec2(bound_dims.x,bound_dims.y), 0.1f);
	primLine(vec2(bound_dims.x,0.f), vec2(bound_dims.x,bound_dims.y), 0.1f);
	primPop();

#if 0
	for (int i = 0; i < cells.count; ++i) {
		if (cells[i].count == 0) continue;
		int x = i % cells.width;
		int y = i / cells.width;
		if (cells[i].count >= cells[i].maxcount)
			primColor(1,0,0);
		else if (cells[i].count > cells[i].maxcount/2)
			primColor(1,1,0);
		else if (cells[i].count > 0)
			primColor(0,0,0);
			
		primBox(world_from_cell_pos.x+x*world_from_cell_scale,world_from_cell_pos.y+y*world_from_cell_scale, 0.9f,0.9f);
	}
#else
	const int* C = cells.ptr;
	const int* C_end = cells.ptr + cells.count;
	while (C != C_end) {
		const int* t = C;
		int count = 0;
		while ((*t >= 0) && (count < CELL_COUNT_MAX)) { ++t; ++count; }	
			
		if (count > 0) {
			if (count >= CELL_COUNT_MAX)
				primColor(1,0,0);
			else if (count > CELL_COUNT_MAX/2)
				primColor(1,1,0);
			else if (count > 0)
				primColor(0,0,0);

			const int i = C-cells.ptr;
			int x = (i%cells_stride)/CELL_COUNT_MAX;
			int y = i / cells_stride;
			primBox(world_from_cell_pos.x+x*world_from_cell_scale,world_from_cell_pos.y+y*world_from_cell_scale, 0.9f,0.9f);
		}
		C += CELL_COUNT_MAX;
	}
#endif
#endif
}