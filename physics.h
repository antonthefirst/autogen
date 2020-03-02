#pragma once

struct SimState;
class TimerTree;
struct Molecule;
struct Body;
struct BodyTemplate;
struct SimStats;

struct StatsReq {
	SimStats* stats;
	bool shell_contents = false;
};

vec3 simBodyCol(SimState*, u8 type);
vec3 simShellStatusCol(int status);
u8 simTypeOfCompound(SimState* sim_state_ptr, u32 types_inside);
float simHeatPumpAvgEnergy(SimState* sim_state_ptr);
Molecule& molGet(SimState* sim_state_ptr, int idx);
Body& bodyGet(SimState* sim_state_ptr, int idx);
const MoleculeType& molTypeGet(const SimState* sim_state_ptr, u8 idx);
u8 typeGetCount(const SimState* sim_state_ptr);
pose WolFromMol(const Molecule& M, const Body& B);
vec2 WolFromMolPos(const Molecule& M, const Body& B);
void shellRadius(const SimState* S, const Body& B, float mol_radii_to_add, float* radius, vec2* center);
void bodyAddMass(Body& B, vec2 bod_from_mass_pos, float mass, float moi, int mass_mol_count);
void bodyAddFree(SimState*);
int bodyAdd(SimState*);
void bodyRem(SimState*, int idx);

int compoundAdd(SimState* S, const BodyTemplate& T, pose wol_from_bod);
int shellAdd(SimState* S, const BodyTemplate& T, pose wol_from_bod);

int getBond(const Molecule& A, int site_idx);
void setBond(Molecule& A, int site_idx, int other_idx);
void fillSite(Molecule& A, int site_idx, int other_idx);
void emptySite(Molecule& A, int site_idx);
int siteIsFull(const Molecule& A, int site_idx);
float computeLinearKineticEnergy(const Body* B);
float computeRotationalKineticEnergy(const Body* B);
float computeRotationalMomentum(const Body* B);

void simFinalize(SimState* sim_state);

void simUpdate(SimState* sim_state, float dt, int substeps, TimerTree* timer, StatsReq stats_req);
void simStatsOnly(SimState* sim_state, TimerTree* timer, StatsReq stats_req);

void drawCells();