/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef CTL_THEORY_H_
#define CTL_THEORY_H_

#include "utils/System.h"
#include "core/Theory.h"
#include "ctl/DynamicKripke.h"
#include "ctl/CTLSolver.h"
#include "ctl/CTLSolverStandalone.h"
#include "DynamicKripke.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "CTLFormula.h"


#include "utils/System.h"
#include "core/Solver.h"

#include "api/Circuit.h"

#include <vector>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <queue>
#include <stack>
#include <map>
#include <stdexcept>
#include <string>


using namespace dgl;
namespace Monosat {


class CTLTheorySolver;

class CTLTheorySolver: public Theory {
public:
	enum VarType { EDGE, NODEAP, DETECTOR};

	struct Transition{
		Var v=var_Undef;
		Var outerVar=var_Undef;
		int from=-1;
		int to=-1;
	};
	struct Assignment {
		VarType type;
		bool assign :1;
		int ID:30; // e.g. EdgeID
		int AP:30; // AP, if it is an assignment to NodeAP, or -1 otherwise
		Var var;
		Assignment(VarType type, bool assign, int ID, int AP, Var v):type(type),assign(assign),ID(ID),AP(AP),var(v){

		}
	};
	double rnd_seed;
	vec<vec<int>> * strings=nullptr;
public: // FIXME was private before!
	Solver * S;
	int local_q = 0;
public:
	int id;

	vec<lbool> assigns;

	vec<Transition> edge_labels;
	DynamicKripke* g_under;
	DynamicKripke* g_over;

	vec<FSMAcceptDetector *> accepts;
	vec<FSMGeneratesDetector *> generates;
	vec<FSMTransducesDetector *> transduces;

	CTLSolver* ctl_under; //one for each over and under
	CTLSolver* ctl_over;
	CTLSolverStandalone* ctl_standalone_under; // Used for a sanity check at the end
	CTLSolverStandalone* ctl_standalone_over;

	CTLFormula* f;
	CTLFormula * original_f;//copy of original formula, without any preprocessing or transformations, for checking the solution
	int initialNode=0; // Start state
	Lit ctl_lit; // the Lit corresponding to the ctl formula

	vec<Lit> symmetryConflict;//this allocates a new symmetryConflict vector at each propagateTheory call, which is needlessly expensive.
	vec<Lit> tmpConflict; // for comparing other conflicts
	vec<Lit> processConflict;

	vec<const char *> property_symbols;
	std::map<std::string, int> property_symbol_map;

	//// Vector of CTL Formulas on the specific Kripke structure this solver deals with
	//vec<CTLFormula> formulas;
	
	vec<Assignment> trail;
	vec<int> trail_lim;
	int next_id=0;

public:
	vec<FSMDetector*> detectors;
	vec<FSMAcceptDetector*> reach_detectors;


	vec<int> marker_map;



	bool requiresPropagation = true;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {
		VarType type;
		int occursPositive :1;
		int occursNegative :1;
		int detector_node_edge :30;	//the detector this variable belongs to, or its edge number, if it is an edge variable, or the node number, if it is a node AP
		int ap; // the AP id, if it is a NodeAP. -1 Otherwise
		Var solverVar;
	};

	vec<VarData> vars;
	vec<vec<int>> nodeAPVarLookup; // Look up the solver var given a node id and a AP id. Initialized by initNodeAPVarLookup(n, a)
	int theory_index = 0;
public:
	
	double mctime = 0;
	double reachtime = 0;
	double unreachtime = 0;
	double pathtime = 0;
	double propagationtime = 0;
	long stats_propagations = 0;
	long stats_real_propagations=0;
	long stats_num_conflicts = 0;
	long stats_decisions = 0;
	long stats_num_reasons = 0;

	double reachupdatetime = 0;
	double unreachupdatetime = 0;
	double stats_initial_propagation_time = 0;
	double stats_decision_time = 0;
	double stats_reason_initial_time = 0;
	double stats_reason_time = 0;
	long num_learnt_paths = 0;
	long learnt_path_clause_length = 0;
	long num_learnt_cuts = 0;
	long learnt_cut_clause_length = 0;
	long stats_pure_skipped = 0;
	long stats_mc_calls = 0;
	long stats_propagations_skipped = 0;

	long stats_symmetry_conflicts =0 ;
	long stats_symmetry_conflict_literals = 0;
	long stats_ctl_conflicts=0;
	long stats_ctl_conflict_literals=0;

	double stats_solve_time =0;
	double stats_theory_prop_symmetry_time=0;
	double stats_theory_prop_ctl_time=0;
	double stats_theory_prop_clause_learning_time=0;

	int bit_over;
	int bit_under;

	// Compute AG EX True
	CTLFormula* fAGEXTrue;

	int processes = 0;
	int statesperprocess = 0;

	CTLTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id){
		g_under = new DynamicKripke(id);
		g_over = new DynamicKripke(id);
		ctl_under = new CTLSolver(id, *g_under, *g_over, 0);
		ctl_over = new CTLSolver(id, *g_over, *g_under, 1);
		ctl_standalone_over = new CTLSolverStandalone(id, *g_over);
		ctl_standalone_under = new CTLSolverStandalone(id, *g_under);
		f = new CTLFormula(True, nullptr, nullptr);

		rnd_seed = opt_random_seed;

		// Compute AG EX True
		CTLFormula* fInfinitePaths1 = newCTLFormula();
		fInfinitePaths1->op = True;
		CTLFormula* fInfinitePaths2 = newCTLFormula();
		fInfinitePaths2->op = EX;
		fInfinitePaths2->operand1 = fInfinitePaths1;
		fAGEXTrue = newCTLFormula();
		fAGEXTrue->op = AG;
		fAGEXTrue->operand1 = fInfinitePaths2;
	}
	~CTLTheorySolver(){
		delete(g_under);
		delete(g_over);
		delete(ctl_under);
		delete(ctl_over);
	}
	
	void writeTheoryWitness(std::ostream& write_to) {

		for (FSMDetector * d : detectors) {
			write_to << "Graph " << this->getGraphID() << ", detector " << d->getID() << ":\n";
			d->printSolution(write_to);
		}
	}
	
	inline int getTheoryIndex() {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}
	inline int getGraphID() {
		return id;
	}

	inline bool isEdgeVar(Var v) const{
		assert(v < vars.size());
		return vars[v].type == EDGE;
	}

	inline bool isNodeAPVar(Var v) const{
		assert(v < vars.size());
		return vars[v].type == NODEAP;
	}

	inline bool isDetectorVar(Var v) const{
		assert(v < vars.size());
		return vars[v].type == DETECTOR;
	}

	inline int getEdgeID(Var v) const {
		assert(isEdgeVar(v));
		return vars[v].detector_node_edge;
	}

	inline int getNodeID(Var v) const {
		assert(isNodeAPVar(v));
		return vars[v].detector_node_edge;
	}

	inline int getAPID(Var v) const {
		assert(isNodeAPVar(v));
		return vars[v].ap;
	}

	inline Transition & getTransition(int edgeID){
			return edge_labels[edgeID];
		}

	inline int getDetector(Var v) const {
		assert(!isEdgeVar(v));
		return vars[v].detector_node_edge;
	}

	inline Var getTransitionVar(int edgeID) {
		Var v = getTransition(edgeID).v;
		assert(v < vars.size());
		return v;
	}

	inline Var getNodeAPVar(int node, int ap) {
		Var v = nodeAPVarLookup[node][ap];
		assert(v < vars.size());

		return v;
	}

	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1, o2);
		S->addClause(o1, ~o2);
	}
	void makeEqualInSolver(Lit l1, Lit l2) {
		S->addClause(~l1, l2);
		S->addClause(l1, ~l2);
	}
	void addClause(Lit l1) {
		Lit o1 = toSolver(l1);
		S->addClause(o1);
	}
	void addClause(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(o1, o2);
	}
	void addClause(Lit l1, Lit l2, Lit l3) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		Lit o3 = toSolver(l3);
		S->addClause(o1, o2, o3);
	}
	void addClause(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		S->addClause(tmp_clause);
	}
	void addClauseSafely(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		
		S->addClauseSafely(tmp_clause);
	}
	
	Var newAuxVar(int forDetector = -1, int ap = -1, VarType type = DETECTOR, bool connectToTheory = false) {
		Var s = S->newVar();
		return newVar(s, forDetector, ap, type, connectToTheory);
	}
	Var newVar(Var solverVar,  int detector_node_edge, int ap, VarType type, bool connectToTheory = true) {
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();

		vars.push();
		vars[v].type = type;
		vars[v].detector_node_edge = detector_node_edge;
		vars[v].ap = ap;
		vars[v].solverVar = solverVar;
		assigns.push(l_Undef);
		if (connectToTheory) {
			S->setTheoryVar(solverVar, getTheoryIndex(), v);
			assert(toSolver(v) == solverVar);
		}

		if(type == EDGE){
			assert(detector_node_edge>-1);
		}

		// FIXME is this right? Don't do anything in particular for DETECTOR
		//if ((type == DETECTOR) && detector_node_edge >= 0)
		//	detectors[detector_node_edge]->addVar(v);
		return v;
	}
	inline int level(Var v) {
		return S->level(toSolver(v));
	}
	inline int decisionLevel() {
		return trail_lim.size(); //S->decisionLevel();
	}
	inline int nVars() const {
		return vars.size(); //S->nVars();
	}
	inline Var toSolver(Var v) {
		//return v;
		assert(v < vars.size()); // FIXME we don't require that vars are assigned without gaps inbetween
		//assert(S->hasTheory(vars[v].solverVar));
		//assert(S->getTheoryVar(vars[v].solverVar)==v);
		return vars[v].solverVar;
	}
	
	inline Lit toSolver(Lit l) {
		//assert(S->hasTheory(vars[var(l)].solverVar));
		//assert(S->getTheoryVar(vars[var(l)].solverVar)==var(l));
		return mkLit(vars[var(l)].solverVar, sign(l));
	}
	
	void toSolver(vec<Lit> & c) {
		for (int i = 0; i < c.size(); i++) {
			c[i] = toSolver(c[i]);
		}
	}
	
	inline lbool value(Var v) {
		if (assigns[v] != l_Undef)
			assert(S->value(toSolver(v)) == assigns[v]);
		
		return assigns[v]; //S->value(toSolver(v));
	}
	inline lbool value(Lit l) {
		if (assigns[var(l)] != l_Undef) {
			assert(S->value(toSolver(l)) == (assigns[var(l)] ^ sign(l)));
		}
		return assigns[var(l)] ^ sign(l);
	}
	inline lbool dbg_value(Var v) {
		return S->value(toSolver(v));
	}
	inline lbool dbg_value(Lit l) {
		return S->value(toSolver(l));
	}
	inline bool enqueue(Lit l, CRef reason) {
		assert(assigns[var(l)]==l_Undef);
		
		Lit sl = toSolver(l);
		if (S->enqueue(sl, reason)) {
			enqueueTheory(l);
			return true;
		} else {
			return false;
		}
	}


	int newNode() {
		g_over->addEmptyState();

		seen.growTo(nNodes());

		return g_under->addEmptyState();
	}
	void newNodes(int n) {
		for (int i = 0; i < n; i++)
			newNode();
	}
	int nNodes() {
		return g_over->nodes();
	}
	bool isNode(int n) {
		return n >= 0 && n < nNodes();
	}

	void setCTL(CTLFormula& myf, int initial) {
		f = &myf;
		initialNode = initial;
	}

	void backtrackUntil(int level) {
		if(opt_verb>1)
			printf("Backtracking until level %d\n", level);

		static int it = 0;
		
		bool changed = false;
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if (e.type == EDGE) {
					//assert(dbg_value(e.var)==l_Undef);
					int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
					assert(edgeID==e.ID);
					if(opt_verb>1)
						printf("Backtracking edge assignment %d\n", edgeID);


					if (e.assign) {
						g_under->disableTransition(edgeID);
					} else {
						g_over->enableTransition(edgeID);
					}
				} else if (e.type == NODEAP) {
					//assert(dbg_value(e.var)==l_Undef);
					int nodeID = getNodeID(e.var);
					int apID = getAPID(e.var);
					assert(nodeID==e.ID);
					assert(apID==e.AP);
					if(opt_verb>1)
						printf("Backtracking nodeap assignment %d\n", nodeID);

					if (e.assign) {
						g_under->disableAPinStateLabel(nodeID, apID);
					} else {
						g_over->enableAPinStateLabel(nodeID, apID);
					}
				} else {
					//This is a detector literal
					//detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
					if(opt_verb>1)
						printf("Backtracking CTL assignment\n");

				}
				assigns[e.var] = l_Undef;
				changed = true;
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			
			if (changed) {
				requiresPropagation = true;
				/*				g.markChanged();
				 antig.markChanged();
				 cutGraph.markChanged();*/
			}
			
			for (FSMDetector * d : detectors) {
				d->backtrack(level);
			}
		}
		
	};

	Lit decideTheory() {
		if (!opt_decide_theories)
			return lit_Undef;
		double start = rtime(1);

		for (int i = 0; i < detectors.size(); i++) {
			FSMDetector * r = detectors[i];
			Lit l = r->decide(decisionLevel());
			if (l != lit_Undef) {
				stats_decisions++;
				r->stats_decisions++;
				stats_decision_time += rtime(1) - start;
				return toSolver(l);
			}
		}
		stats_decision_time += rtime(1) - start;
		return lit_Undef;
	}
	
	void backtrackUntil(Lit p) { printf("WARNING: Used dummy function backtrackUntil(Lit) in CTLTheory.h"); } // dummy function

	/*
	void backtrackUntil(Lit p) { // NOT actually important
		//need to remove and add edges in the two graphs accordingly.
		assert(value(p)==l_True);
		int i = trail.size() - 1;
		for (; i >= 0; i--) {
			Assignment e = trail[i];
			if (var(p) == e.var) {
				assert(sign(p) != e.assign);
				break;
			}
			if (e.type == EDGE) {
				int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
					g_under->disableTransition(edgeID);
				} else {
					g_over->enableTransition(edgeID);
				}
			} else {

				assigns[e.var] = l_Undef;
				detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
			}
		}
		
		trail.shrink(trail.size() - (i + 1));
		//if(i>0){
		requiresPropagation = true;
		//			g.markChanged();
		// antig.markChanged();
		// cutGraph.markChanged();
		//}

	/*
		for (FSMDetector * d : detectors) {
			d->backtrack(this->decisionLevel());
		}
	*/
		//while(trail_lim.size() && trail_lim.last()>=trail.size())
		//	trail_lim.pop();

		/*		for(int i = 0;i<reach_detectors.size();i++){
		 if(reach_detectors[i]->positive_reach_detector)
		 reach_detectors[i]->positive_reach_detector->update();
		 if(reach_detectors[i]->negative_reach_detector)
		 reach_detectors[i]->negative_reach_detector->update();
		 }*/
	/*
	}
	;
*/

	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason) { // IGNORED, not essential according to Sam
		/*
		CRef marker = S->reason(var(toSolver(p)));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef - marker;
		int d = marker_map[pos];
		//double initial_start = rtime(1);
		double start = rtime(1);
		backtrackUntil(p);
		
		assert(d < detectors.size());
		detectors[d]->buildReason(p, reason, marker);
		toSolver(reason);
		double finish = rtime(1);
		stats_reason_time += finish - start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;
		 */
		
	}
	


	
	bool dbg_graphsUpToDate() {

		return true;
	}

	void prepare_ids(CTLFormula *f){
		if(f->id<0){
			f->setID(next_id++);
			if(f->operand1 && f->operand1->id==-1){
				prepare_ids(f->operand1);
			}
			if(f->operand2 && f->operand2->id==-1){
				prepare_ids(f->operand2);
			}
		}

	}
	void cnfSymmetryBreaking(){
		//Every node must have an outgoing edge (which may be a self loop).
		Circuit<Solver> c(*S);
		Lit True = c.True();
		Lit False = c.False();
			if(opt_ctl_symmetry_cnf==1 || opt_ctl_symmetry_cnf==3){
				//Node 0 is reachable.
				//Every other node is reachable exactly if it has an incoming edge from a lower indexed node.
				assert(initialNode==0);
				//If a node is not reachable, then it has exactly one in/out edge, which is the self loop.
				//if node i>0 is reachable, then node i-1 must also be reachable
				vec<Lit> nodes_reachable;
				nodes_reachable.push(True);//node 0 is always reachable
				for(int i = 1;i<g_over->nodes();i++){
					Lit reachable = False;
					for(int j = 0;j<g_over->nIncoming(i);j++){
						int edgeID = g_over->incoming(i,j).id;
						int to = g_over->incoming(i,j).node;
						if(to<i){
							Lit edge_enabled = mkLit(this->toSolver(getTransitionVar(edgeID)));
							reachable = c.Or(edge_enabled,reachable);
							/*Lit reachable_or_edge_enabled = mkLit(S->newVar());
							S->addClause(reachable_or_edge_enabled,~edge_enabled);
							S->addClause(reachable_or_edge_enabled,~reachable);
							S->addClause(~reachable_or_edge_enabled,edge_enabled,reachable);*/
							//reachable=reachable_or_edge_enabled;
						}
					}
					nodes_reachable.push(reachable);
				}

				//if a node is not reachable, then all incoming and outgoing edges, except for the self loop, are disabled, and the self loop edge is enabled
				for(int i = 1;i<g_over->nodes();i++){
					Lit r = nodes_reachable[i];
					for(int j = 0;j<g_over->nIncoming(i);j++){
						int edgeID = g_over->incoming(i,j).id;
						int to = g_over->incoming(i,j).node;
						if(to>i){
							Lit edge_enabled = mkLit(this->toSolver(getTransitionVar(edgeID)));
							S->addClause(r,~edge_enabled);//if this node is _not_ reachable, then this incoming edge must be disabled.
							//note that the case for to<i doesn't need to be checked here, since it was checked above.
						}else if (to==i){
							//this is the self loop edge. If this node is not reachable, then force it to have a self loop
							Lit edge_enabled = mkLit(this->toSolver(getTransitionVar(edgeID)));
							S->addClause(r,edge_enabled);
						}
					}
					for(int j = 0;j<g_over->nIncident(i);j++){
						int edgeID = g_over->incident(i,j).id;
						int to = g_over->incident(i,j).node;
						if(to!=i){
							Lit edge_enabled = mkLit(this->toSolver(getTransitionVar(edgeID)));
							S->addClause(r,~edge_enabled);//if this node is _not_ reachable, then this outgoing edge must be disabled.
							//note that the case for to<i doesn't need to be checked here, since it was checked above.
						}
					}
				}

				//if node i is reachable, then node i-1 must also be reachable (can skip node 1 here, since node 0 is always reachable).
				for(int i = 2;i<g_over->nodes();i++){
					Lit r = nodes_reachable[i];
					Lit prev_r = nodes_reachable[i-1];
					S->addClause(~r, prev_r);//node i being reachable implies that the previous node must be reachable.
				}
			}
			if(opt_ctl_symmetry_cnf==2 || opt_ctl_symmetry_cnf==3){
				assert(initialNode==0);//haven't generalized this yet
				//interpret the property assignment of each state i as bitvector bv(i), and enforce that bv(i) <= bv(i+1).
				//but ignore the starting state in these constraints.
				vec<vec<Lit> > state_assignments;
				for(int i = 0;i<g_over->nodes();i++){
					state_assignments.push();
					vec<Lit> & state_assignment = state_assignments.last();
					for(int j =0;j<g_over->nProperties();j++){
						Var v= this->getNodeAPVar(i,j);
						Lit l = mkLit(toSolver(v));
						state_assignment.push(l);
					}
				}

				//enforce that each state assignment is less than its predecessor
				for(int i = 2;i<state_assignments.size();i++){
					vec<Lit> & p = state_assignments[i-1];
					vec<Lit> & n = state_assignments[i];
					assert(p.size()==n.size());
					//starting from the most significant bit, check if they are all the same
					Lit p_lt_n=False;
					Lit p_eq_n = True;
					for(int j = p.size()-1;j>=0;j--){
						//Lit pj_eq_n =  mkLit(S->newVar());
						Lit pj_eq_n = c.Xnor(n[j],p[j]);
						Lit pj_ltn = c.And(n[j],p[j]);
						Lit lt = c.And(p_eq_n,pj_ltn);
						p_lt_n = c.Or(lt,p_lt_n);
						p_eq_n=c.And(pj_eq_n,p_eq_n);
					}
					//The bv(i-1) must be less or equal to bv(i)
					S->addClause(p_eq_n,p_lt_n);
				}
				if(opt_ctl_symmetry_statelabelandedges>0){
					//if p_eq_n, then also enforce symmetry constraints on the edges

				}
			}
	}
	void preprocess() {
		original_f = f->copy();
		if(opt_ctl_symmetry_cnf>0){
			cnfSymmetryBreaking();
		}
		for (int i = 0; i < detectors.size(); i++) {
			detectors[i]->preprocess();
		}

		// AG EX True
		// This is somewhat slow, but whatever. Only done once in preprocessing
		vec<Lit> c;
		/*
		for (int i = 0; i < g_over->states(); i++) { // iterate over neighbours of current front of queue
			for (int j = 0; j < g_over->nIncident(i); j++) { // iterate over neighbours of current front of queue
				e = g_over->incident(i, j);

				Lit l = ~mkLit(getTransitionVar(e.id), true);
				c.push(l);
			}
			addClauseSafely(c);
			c.clear();
		}
		*/

		if(opt_force_all_states_reachable){
			if (initialNode!=0){
				throw std::runtime_error("Enforcing that all nodes are reachable currently requires the initial node to be 0.");
			}
			//Enforcing that every node is reachable by forcing all nodes to have at least one incoming transition from a lower
			//indexed node. Watchout - this may interact incorrectly with symmetry breaking!
			vec<Lit> tmp;
			for(int i = 1;i<g_over->states();i++){
				tmp.clear();
				for(int j = 0;j<g_over->nIncoming(i);j++){
					int edgeID = g_over->incoming(i,j).id;
					Lit l = mkLit(getTransitionVar(edgeID));
					tmp.push(toSolver(l));
				}
				S->addClause(tmp);
			}
		}

		if(opt_optimize_formula>=1){
			Circuit<Solver> c(*S);
			Lit hasSomeIncomingEdgeLit[g_over->statecount];
			for (int s = 0; s<g_over->statecount; s++) {
				if(s!=initialNode){
					vec<Lit> incomingEdges;
					//printf("hasSomeIncomingEdgeLit[%d] = ", s);
					for(int i = 0;i<g_over->nIncoming(s);i++){
						int edgeID = g_over->incoming(s,i).id;
						incomingEdges.push( toSolver( mkLit(this->getTransitionVar(edgeID))) );
							//printf("%d (%d) %d -> %d || ", toSolver( mkLit(this->getTransitionVar(edgeID))) , this->getTransitionVar(edgeID), g_over->getEdge(edgeID).from, g_over->getEdge(edgeID).to);
					}
					hasSomeIncomingEdgeLit[s] = c.Or(incomingEdges);
				}else{
					hasSomeIncomingEdgeLit[s]=c.True();//the initial node is always reachable.
				}
				//printf("  has var %d\n", hasSomeIncomingEdgeLit[s]);
			}
			if(opt_verb>1){
				printf("OPTIMIZING FORMULA\n");
				printFormula(f);printf("\n");
			}
			optimizeFormula(opt_force_all_states_reachable, c, hasSomeIncomingEdgeLit);
			if(opt_verb>1){
				printf("Done optimizing\n");
				printFormula(f);printf("\n");
			}
		}

		// Each process is in exactly one process-state
		int ap;
		if (opt_ctl_process_in_single_state && processes > 0) {
			for (int s = 0; s < g_over->states(); s++) {         // state
				for (int p = 0; p < processes; p++) {            // process
					// Each process has to be in one of statesperprocess process-states
					for (int j = 0; j < statesperprocess; j++) { // process-state
						ap = p*statesperprocess + j;
						Lit l = ~mkLit(getNodeAPVar(s, ap), true);
						c.push(l);
					}
					if (opt_verb > 1)
						printLearntClause(c);
					addClauseSafely(c);
					c.clear();
					// For every combination of two process-states such that they're not the same, add a clause with two literals stating that they can't be both true
					for (int i = 0; i < statesperprocess; i++) { // process-state
						for (int j = 0; j < statesperprocess; j++) { // process-state
							if (i != j) {
								ap = p*statesperprocess + i;
								Lit l1 = ~mkLit(getNodeAPVar(s, ap), false);
								c.push(l1);
								ap = p*statesperprocess + j;
								Lit l2 = ~mkLit(getNodeAPVar(s, ap), false);
								c.push(l2);
								if (opt_verb > 1)
									printLearntClause(c);
								addClauseSafely(c);
								c.clear();
							}
						}
					}
				}
			}
		}

		if(opt_ctl_learn_cache){
			prepare_ids(f);
		}

		if (opt_verb > 0) {
			printf("Formula:\n");
			printFormula(f);
			printf("\n");
		}
	}

	// Does optimizations on the formula
	// Right now, it considers all top-level AG phi constraints where phi is propositional. Top level means that the parent operators are only AND
	void optimizeFormula(bool all_nodes_reachable, Circuit<Solver>  & c, Lit hasSomeIncomingEdgeLit[]) {
		CTLFormula* position = f;
		optimizeFormulaRec(position,c,all_nodes_reachable,opt_optimize_formula>=2, hasSomeIncomingEdgeLit);
	}
	bool optimizeFormulaRec(CTLFormula* position,Circuit<Solver>  & c, bool all_nodes_reachable, bool clausify_non_nested_x, Lit hasSomeIncomingEdgeLit[]) {
/*		if (position->op == NEG) {
			printf("optimizeFormulaRec: has neg: "); printFormula(position); printf("\n");
		}*/
		bool removed=false;
		//This code should really be generalized to allow an arbitrary propositional CTL formula above AG/EF statements
		if (position->isPropositional(clausify_non_nested_x)){
			//purely propositional statements hold on the intitial node only
			Lit result = AGtoCNF(position, initialNode, &c);
			S->addClause(result);
			removed=true;
		}else if (position->op == AG && position->operand1->isPropositional(clausify_non_nested_x)) {
			if(opt_verb>1){
				printf("optimizeFormulaRec: optimizing \n"); printFormula(position);
			}
			for (int s = 0; s<g_over->statecount; s++) {
				Lit result = AGtoCNF(position->operand1, s, &c);
				if (opt_optimize_formula == 3 && s != initialNode)
					S->addClause( c.Or(result, ~hasSomeIncomingEdgeLit[s]) );
				else
					S->addClause(result);
			}
			removed=true;
		}else if (all_nodes_reachable && position->op == NEG && position->operand1->op == AG &&  position->operand1->operand1->isPropositional(clausify_non_nested_x) ) {

				if(opt_verb>1){
					printf("optimizeFormulaRec: optimizing \n"); printFormula(position);
				}
				vec<Lit> tmp;
				for (int s = 0; s<g_over->statecount; s++) {
					Lit result = ~AGtoCNF(position->operand1, s, &c);
					tmp.push(result);
				}
				//WARNING: this is only safe if we also enforce that all nodes are reachable
				S->addClause(tmp);
				removed=true;

		}else if (position->op == NEG && position->operand1->op == EF &&  position->operand1->operand1->isPropositional(clausify_non_nested_x)) {
			// AG phi = ~EF~ phi
			if(opt_verb>1){
				printf("optimizeFormulaRec: optimizing \n"); printFormula(position);
			}
			for (int s = 0; s<g_over->statecount; s++) {
				Lit result = ~AGtoCNF(position->operand1->operand1, s, &c);
				if (opt_optimize_formula == 3 && s != initialNode) {
					//printf("Added clause: %d or ~ %d\n ", result.x, (~hasSomeIncomingEdgeLit[s]).x);
					S->addClause( c.Or(result, ~hasSomeIncomingEdgeLit[s]) );
				} else
					S->addClause(result);
			}
			removed=true;
		}else if (all_nodes_reachable && position->op == EF &&  position->operand1->isPropositional(clausify_non_nested_x)) {
			// AG phi = ~EF~ phi
			if(opt_verb>1){
				printf("optimizeFormulaRec: optimizing \n"); printFormula(position);
			}
			vec<Lit> tmp;
			for (int s = 0; s<g_over->statecount; s++) {
				Lit result = AGtoCNF(position->operand1, s, &c);
				tmp.push(result);
			}
			//WARNING: this is only safe if we also enforce that all nodes are reachable
			S->addClause(tmp);
			removed=true;
		}else if (position->op == AND) {
			removed = optimizeFormulaRec(position->operand1,c,all_nodes_reachable,clausify_non_nested_x, hasSomeIncomingEdgeLit);
			removed &= optimizeFormulaRec(position->operand2,c,all_nodes_reachable,clausify_non_nested_x, hasSomeIncomingEdgeLit);
		}

		if(removed){
			position->op = True;
			if(position->operand1){
				delete(position->operand1);position->operand1=nullptr;
			}
			if(position->operand2){
				delete(position->operand2);position->operand2=nullptr;
			}
		}
		return removed;
	}
	vec<Lit> ag_to_cnf_tmp;
	// AG phi for phi propositional. Convert this into a literal using the given circuit.
	Lit AGtoCNF(CTLFormula* inner, int state, Circuit<Solver>* c) {
		if (opt_verb>1){
			printf("optimizeFormulaRec: at  \n"); printFormula(inner);
		}
		if (inner->op == AND) {
			return c->And(AGtoCNF(inner->operand1, state, c), AGtoCNF(inner->operand2, state, c));
		} else if (inner->op == OR) {
			return c->Or(AGtoCNF(inner->operand1, state, c), AGtoCNF(inner->operand2, state, c));
		} else if (inner->op == NEG) {
			return ~AGtoCNF(inner->operand1, state, c);
		} else if (inner->op == True) {
			return c->True();
		} else if (inner->op == ID) {
			return toSolver(mkLit(getNodeAPVar(state, inner->value)));//Sam: Changed this from mkLit(v, True), as that creates a negated literal. Also translating the variable from the theory's to the solver's namespace.
		}else if (inner->op == EX) {
			assert(inner->operand1->isPropositional());
			ag_to_cnf_tmp.clear();
			for(int i = 0;i<g_over->nIncident(state);i++){
				int edgeID = g_over->incident(state,i).id;
				int to = g_over->incident(state,i).node;
				Lit l =toSolver( mkLit(this->getTransitionVar(edgeID)));
				ag_to_cnf_tmp.push(c->And(l,AGtoCNF(inner->operand1,to,c)));
			}
			return c->Or(ag_to_cnf_tmp);
		}else if (inner->op == AX) {
			assert(inner->operand1->isPropositional());
			ag_to_cnf_tmp.clear();
			for(int i = 0;i<g_over->nIncident(state);i++){
				int edgeID = g_over->incident(state,i).id;
				int to = g_over->incident(state,i).node;
				Lit l =toSolver( mkLit(this->getTransitionVar(edgeID)));
				ag_to_cnf_tmp.push(c->And(l,AGtoCNF(inner->operand1,to,c)));
			}
			return c->And(ag_to_cnf_tmp);
		}  else {
			assert(false); // we only wanted propositional formulas!
		}
	}

	void setLiteralOccurs(Lit l, bool occurs) {
		if (isEdgeVar(var(l))) {
			//don't do anything
		} else {
			//this is a graph property detector var
			if (!sign(l) && vars[var(l)].occursPositive != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
			else if (sign(l) && vars[var(l)].occursNegative != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
		}

	}
	const char *getPropertySymbol(int property){
		if(property>=0 && property<this->property_symbols.size()){
			return this->property_symbols[property];
		}else{
			return "";
		}
	}
	int getPropertyFromSymbol(const char * symbol){
	/*	string s(symbol);
		if(property_symbol_map.count(s)){
			return property_symbol_map[s];
		}else{
			return -1;
		}*/
		return -1;
	}

	void enqueueTheory(Lit l) {
		Var v = var(l);
		
		int lev = level(v);
		
		assert(decisionLevel() <= lev);
		
		while (lev > trail_lim.size()) {
			newDecisionLevel();
		}
		
		if (assigns[var(l)] != l_Undef) {
			return;			//this is already enqueued.
		}
		if(opt_verb>1)
			printf("enqueueTheory(Lit %d), var %d, new sign: %d\n", (l.x+2), v+1, sign(l));

		assert(assigns[var(l)]==l_Undef);
		assigns[var(l)] = sign(l) ? l_False : l_True;
		requiresPropagation = true;
		//printf("enqueue %d\n", dimacs(l));
		
#ifndef NDEBUG
		{
			for (int i = 0; i < trail.size(); i++) {
				assert(trail[i].var != v);
			}
		}
#endif

		if (isEdgeVar(var(l))) {
			
			//this is an edge assignment
			int edgeID = getEdgeID(var(l)); //v-min_edge_var;
			assert(getTransition(edgeID).v == var(l));

			trail.push( { EDGE, !sign(l),edgeID, -1, v });

			if (!sign(l)) {
				g_under->enableTransition(edgeID);
			} else {
				g_over->disableTransition(edgeID);
			}
			
		} else if (isNodeAPVar(var(l))){
			int nodeID = getNodeID(var(l));
			int apID = getAPID(var(l));
			trail.push( { NODEAP, !sign(l), nodeID, apID, v });
			if (!sign(l)) {
				g_under->enableAPinStateLabel(nodeID, apID);
			} else {
				g_over->disableAPinStateLabel(nodeID, apID);
			}

		} else { // Detector AP
			trail.push( { DETECTOR, !sign(l),-1, -1, v });
			//this is an assignment to a non-edge atom. (eg, a reachability assertion)
			//detectors[getDetector(var(l))]->assign(l);
		}
		
	}
	;
	bool propagateTheory(vec<Lit> & conflict) override{
		return propagateTheory(conflict,false);
	}

	bool propagateTheory(vec<Lit> & conflict, bool forcePropagation) {

		stats_propagations++;

		if (!requiresPropagation) {
			stats_propagations_skipped++;
			assert(dbg_graphsUpToDate());
			return true;
		}
		stats_real_propagations++;

		bool skip_this_propagation = decisionLevel()>0 && ! forcePropagation  &&  (stats_real_propagations %((long)opt_ctl_skip_prop));
		if(skip_this_propagation){
			stats_propagations_skipped++;
			return true;//skip this theory propagation round
		}
		
		bool any_change = false;
		double startproptime = rtime(1);
		//static vec<int> detectors_to_check;
		
		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.
		
		//At level 0, need to propagate constant reaches/source nodes/edges...
		
		//stats_initial_propagation_time += rtime(1) - startproptime;

		assert(dbg_graphsUpToDate());

		if(opt_verb>1){
			printf("\n--------------------\npropagateTheory has been called\n");
			printf("\nOver:\n");
			g_over->draw();
			printf("Under:\n");
			g_under->draw();
		}
		/*
		//should leak exactly 2 bitsets ('bit_under' and 'bit_over') in the above calls
		long bitsets = Bitset::remainingBitsets();

		Bitset* bit_under = ctl_under->solveFormula(*f);
		Bitset* bit_over = ctl_over->solveFormula(*f);
		long remainingbitsets Bitset::remainingBitsets()-bitsets;
		if(remainingbitsets !=2){
			throw std::runtime_error("Leaked bitsets during solve call!");
		}
		*/

		if (value(ctl_lit)==l_False) {
			throw std::runtime_error("ctl_lit is assigned false by the SAT solver. So far this is not supported. Please negate the CTL formula and assign ctl_lit true");
		}
		if (value(ctl_lit)==l_Undef) {
			throw std::runtime_error("ctl_lit is unassigned by the SAT solver. This is not permitted?");
		}

		// Enforce that on every transition, exactly one process-state changes
		if (opt_ctl_only_one_process_moves > 0 && processes > 0) {
			processConflict.clear();
			checkExactlyOneProcessChangesOnTransition(processConflict);
			if (processConflict.size() > 0) {
				if (opt_verb>1) {
					printf("Process in Single State:\n");
					printLearntClause(processConflict);
				}
				if (opt_ctl_only_one_process_moves == 1) { // Always prefer single-process-moves clauses
					processConflict.copyTo(conflict);
					if(opt_verb>1)
						printf("propagateTheory: processConflict conflict is learned, since it is always preferred.\n");
					toSolver(conflict);
					theoryPropagationAppendix(startproptime);
					return false;
				}
			}
		}

		/*
		 * Try to find a clause to learn. With symmetry reduction disabled, this will be a CTL clause. Otherwise, it may also be a symmetry clause
		 * Which clause is preferred depends on the mode in opt_ctl_symmetry.
		 *
		 * There are different modes in opt_ctl_symmetry.
		 * 0: No symmetry reduction
		 * 1-3: NodeAP Symmetry (1: learn first symmetry clause, 2: learn smallest symmetry clause, 3: learn smallest symmetry/ctl clause)
		 * 4-6: Edge Symmetry (4: learn first edge symmetry clause, 5: learn smallest symmetry clause, 6: learn smallest symmetry/ctl clause)
		 */
		ctl_over->clearCache();
		ctl_under->clearCache();
		if (opt_ctl_symmetry > 0) { // >0 means we have symmetry reduction
			double start_time = rtime(2);
			symmetryConflict.clear();
			if (opt_ctl_symmetry <= 3) // reduce on nodeAPs
				checkStateLabelSymmetryConstraints(symmetryConflict, initialNode); // check if there is a conflict with the symmetry constraints under the current assignment, and if there is, build a clause to learn
			else // reduce on number of outgoing edges
				checkEdgeSymmetryConstraints(symmetryConflict, initialNode);
			minimizeClause(symmetryConflict);
			stats_theory_prop_symmetry_time+= rtime(2)-start_time;
			// Prefer the symmetry conflict without even of looking for a CTL conflict if:
			//  - Symmetry reduction is enabled
			//  - and, there is a symmetry conflict
			//  - and, we are are in modes 1, 2, 4, or 5

			if ((symmetryConflict.size() != 0) && (opt_ctl_symmetry < 3 || opt_ctl_symmetry == 4 || opt_ctl_symmetry == 5)) {
				symmetryConflict.copyTo(conflict); // Overwrite conflict

				// Overwrite with process conflict, if it is smaller.
				if (opt_ctl_only_one_process_moves > 1 && processConflict.size() > 0) {
					minimizeClause(processConflict);
					if (processConflict.size() < symmetryConflict.size())
						processConflict.copyTo(conflict);
				}

				toSolver(conflict);
				stats_symmetry_conflicts++;
				stats_symmetry_conflict_literals+=conflict.size();
				if(opt_verb>1){
					printFullClause();
			  		printf("Symmetry Clause:\n");
			  		printLearntClause(symmetryConflict);
					printf("propagateTheory returns false. Choosing to learn symmetry conflict rather than CTL conflict, since symmetry clauses are preferred due to opt_ctl_symmetry\n");
				}
				theoryPropagationAppendix(startproptime);
				return false;
			}

			// Otherwise, we have to look for a CTL conflict, too.
			double start_ctl_prop_time = rtime(2);
			bit_over = ctl_over->solveFormula(*f);
			stats_solve_time+=rtime(2)-start_ctl_prop_time;
			if (!ctl_over->bitsets[bit_over]->operator [](initialNode)) { // we know value(ctl_lit)==l_True
				double start_ctl_clause_learning_time = rtime(2);
				learnClausePos(conflict, *f, initialNode);
				minimizeClause(conflict);
				stats_theory_prop_clause_learning_time+=rtime(2)-start_ctl_clause_learning_time;
				if(opt_verb>1){
					printFullClause();
					printf("CTL Clause:\n");
					printLearntClause(conflict);
					printf("Symmetry Clause:\n");
					printLearntClause(symmetryConflict);
				}
				if (symmetryConflict.size()<conflict.size() && (symmetryConflict.size() != 0)) {
					if(opt_verb>1)
						printf("propagateTheory returns false. Choosing to learn symmetry conflict rather than CTL conflict, since it is smaller (%d literals vs %d literals)\n", symmetryConflict.size(), conflict.size());
					symmetryConflict.copyTo(conflict); // Overwrite conflict
				} else {
					if(opt_verb>1)
						printf("propagateTheory returns false, since formula is asserted true, but fails to hold in the overapproximation (and hence also fails to hold in the underapproximation) \n");
				}
				// Overwrite with process conflict, if it is smaller.
				if (opt_ctl_only_one_process_moves > 1 && processConflict.size() > 0) {
					minimizeClause(processConflict);
					if (processConflict.size() < conflict.size()) {
						processConflict.copyTo(conflict);
						if(opt_verb>1)
							printf("NOT, actually the processConflict conflict is learned.\n");
					}
				}
				toSolver(conflict);
				ctl_over->freeBitset(bit_over);
				stats_ctl_conflicts++;
				stats_ctl_conflict_literals+=conflict.size();
				theoryPropagationAppendix(startproptime);
				return false;
			}
			ctl_over->freeBitset(bit_over);
			// There is no CTL conflict, so just learn the symmetry conflict, if there is one
			if (symmetryConflict.size() != 0) {
				symmetryConflict.copyTo(conflict); // Overwrite conflict
				// Overwrite with process conflict, if it is smaller.
				if (opt_ctl_only_one_process_moves > 1 && processConflict.size() > 0) {
					minimizeClause(processConflict);
					if (processConflict.size() < symmetryConflict.size()) {
						processConflict.copyTo(conflict);
						if(opt_verb>1)
							printf("propagateTheory: processConflict conflict is learned.\n");
					}
				}
				toSolver(conflict);
				stats_symmetry_conflicts++;
				stats_symmetry_conflict_literals+=conflict.size();
				theoryPropagationAppendix(startproptime);
				return false;
			} else if (opt_ctl_only_one_process_moves > 1 && processConflict.size() > 0) {
				processConflict.copyTo(conflict);
				if(opt_verb>1)
					printf("propagateTheory: processConflict conflict is learned.\n");
				toSolver(conflict);
				theoryPropagationAppendix(startproptime);
				return false;
			}
		} else { // Ignore symmetry
			double start_time = rtime(2);
			bit_over = ctl_over->solveFormula(*f);
			stats_solve_time+=rtime(2)-start_time;

			if (!ctl_over->bitsets[bit_over]->operator [](initialNode)) { // we know value(ctl_lit)==l_True
				double start_ctl_clause_learning_time = rtime(2);
				learnClausePos(conflict, *f, initialNode);
				stats_theory_prop_clause_learning_time+=rtime(2)-start_ctl_clause_learning_time;
				minimizeClause(conflict);
				// Overwrite with process conflict, if it is smaller.
				if (opt_ctl_only_one_process_moves > 1 && processConflict.size() > 0) {
					minimizeClause(processConflict);
					if (processConflict.size() < symmetryConflict.size()) {
						processConflict.copyTo(conflict);
						if(opt_verb>1)
							printf("propagateTheory: processConflict conflict is learned.\n");
					}
				}
				toSolver(conflict);
				stats_ctl_conflicts++;
				stats_ctl_conflict_literals+=conflict.size();
				if(opt_verb>1){
					printFullClause();
			  		printf("CTL Clause:\n");
					printLearntClause(conflict);
			  		printf("Symmetry Clause:\n");
			  		printLearntClause(symmetryConflict);
					printf("propagateTheory returns false, since formula is asserted true, but fails to hold in the overapproximation (and hence also fails to hold in the underapproximation) \n");
				}
				ctl_over->freeBitset(bit_over);
				theoryPropagationAppendix(startproptime);
				return false; // It does not hold in the overapproximation
			} else if (opt_ctl_only_one_process_moves > 1 && processConflict.size() > 0) {
				processConflict.copyTo(conflict);
				if(opt_verb>1)
					printf("propagateTheory: processConflict conflict is learned.\n");
				toSolver(conflict);
				ctl_over->freeBitset(bit_over);
				theoryPropagationAppendix(startproptime);
				return false;
			}
			ctl_over->freeBitset(bit_over);
		}

		// ALLSAT
		// We can use the following code to simply print out all models to the CTL formula.
		// We do this by checking if the formula holds in the underapproximation. If so, then it must hold in every extension.
		// We learn the naive clause so that we don't get the exact same case again
		// FIXME this is broken
		if (opt_all_solutions > 0) {
			bit_under = ctl_under->solveFormula(*f);
			if (ctl_under->bitsets[bit_under]->operator [](initialNode)) { // It is true in the underapproximation, therefore it is a solution
				learnNaiveClause(conflict, initialNode);
				toSolver(conflict);
				printf("--------------------\nFound new Solution set (conflict size: %d):\n", conflict.size());
				drawCurrentAssignment();
				ctl_under->freeBitset(bit_under);
				theoryPropagationAppendix(startproptime);
				return false;
			} else {
				ctl_under->freeBitset(bit_under);
			}
		}
		theoryPropagationAppendix(startproptime);
		return true;
	}
	void theoryPropagationAppendix(int startproptime) {
		requiresPropagation = false;
		g_under->clearChanged();
		g_over->clearChanged();

		g_under->clearHistory();
		g_over->clearHistory();

		//detectors_to_check.clear();
		
		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;
	}
	;

	// Sort and remove duplicates in the clause
	// (This is important to call here, rather than in the SAT solver, so that multiple candidate clauses can be compared and the smallest returned)
	void minimizeClause(vec<Lit> & ps) {
    	sort(ps);
    	int i, j;
    	Lit p;
    	for (i = j = 0, p = lit_Undef; i < ps.size(); i++) {
    		if (ps[i] != p) {
    			ps[j++] = p = ps[i];
    		}
    	}
    	ps.shrink(i - j);
	}

	/*
	 * Make sure that exactly one process changes its local state on a transition, and all other processes stay the same.
	 */
	void checkExactlyOneProcessChangesOnTransition(vec<Lit> & conflict) {
		//printf("checkProcessInSingleState\n");
		int from, to;
		DynamicGraph<int>::FullEdge e;
		int fromP, toP;
		int fromP2, toP2;
		int processThatMoved;
		for (int edge = 0; edge < g_under->edges(); edge++) {         // all edges
			if (g_under->edgeEnabled(edge)) {
				processThatMoved = -1;
				e = g_under->getEdge(edge);
				assert (edge == e.id);
				from = e.from;
				to = e.to;

				// Check if both state labels are completely determined, and also completely equivalent. If so, then no process moves, which is a violation
				if (g_under->statelabel[from]->Equiv(*g_over->statelabel[to]) && g_under->statelabel[from]->Equiv(*g_over->statelabel[to])) {
					learnAtLeastOneProcessChangesOnTransition(conflict, e.id, from, to);
				}

				//printf("checkProcessInSingleState: edge %d->%d\n", from, to);

				// Find the first process that changed its process-state in to vs from
				for (int p = 0; p < processes && processThatMoved < 0; p++) {
					fromP = -1;
					toP = -1;
					//printf("checkProcessInSingleState: edge %d->%d, process %d\n", from, to, p);

					// Find out which process-state this process is in, both for from and for to
					for (int j = 0; j < statesperprocess && (fromP < 0 || toP < 0); j++) { // process-state
						if (g_under->isAPinStateLabel(from, p*statesperprocess + j)) {
							fromP = j;
							//printf("checkProcessInSingleState: edge %d->%d. process %d set fromP to %d\n", from, to, p, fromP);

						}
						if (g_under->isAPinStateLabel(to, p*statesperprocess + j)) {
							toP = j;
							//printf("checkProcessInSingleState: edge %d->%d. process %d set toP to %d\n", from, to, p, toP);
						}
					}
					// Register if process p has moved
					if (fromP >= 0 && toP >= 0 && fromP != toP) {
						processThatMoved = p;
						//printf("checkProcessInSingleState: edge %d->%d. first move in process %d that went from %d to %d\n", from, to, p, fromP, toP);
					}
				}

				// If we found a first process that changed, try to find a second one
				if (processThatMoved >=0) {
					for (int p = processThatMoved + 1; p < processes; p++) {
						fromP2 = -1;
						toP2 = -1;
						//printf("checkProcessInSingleState: edge %d->%d, second iteration, process %d\n", from, to, p);
						// Find out which process-state this process is in, both for from and for to
						for (int j = 0; j < statesperprocess && (fromP2 < 0 || toP2 < 0); j++) { // process-state
							if (g_under->isAPinStateLabel(from, p*statesperprocess + j)) {
								fromP2 = j;
							}
							if (g_under->isAPinStateLabel(to, p*statesperprocess + j)) {
								toP2 = j;
							}
						}
						// Register if process p has moved
						if (fromP2 >= 0 && toP2 >= 0 && fromP2 != toP2) {
							//printf("checkProcessInSingleState: edge %d->%d. SECOND move in process %d that went from %d to %d\n", from, to, p, fromP2, toP2);
							// we know there is a conflict. Learn it
							learnAtMostOneProcessChangesOnTransition(conflict, e.id, processThatMoved, fromP, toP, p, fromP2, toP2);
							return;
						}
					}
				}

			}
		}
	}

	/*
	 * Both states have exactly the same label. We learn that one of the enabled APs has to be disabled (which translates to one process being in a different local state)
	 */
	void learnAtLeastOneProcessChangesOnTransition(vec<Lit> & conflict, int eid, int from, int to) {
		for (int i = 0; i < g_under->statelabel[from]->size(); i++) {
			if (g_under->isAPinStateLabel(from, i)) {
				Lit l1 = ~mkLit(getNodeAPVar(from, i), false);
				conflict.push(l1);
				assert(value(l1)==l_False);

				Lit l2 = ~mkLit(getNodeAPVar(to, i), false);
				conflict.push(l2);
				assert(value(l2)==l_False);
			}
		}

		Lit l = ~mkLit(getTransitionVar(eid), false);
		conflict.push(l);
		assert(value(l)==l_False);
	}

	void learnAtMostOneProcessChangesOnTransition(vec<Lit> & conflict, int eid, int p, int fromP, int toP, int p2, int fromP2, int toP2) {
		Lit l1 = ~mkLit(getNodeAPVar(g_under->getEdge(eid).from, p*statesperprocess + fromP), false);
		conflict.push(l1);
		assert(value(l1)==l_False);
		Lit l2 = ~mkLit(getNodeAPVar(g_under->getEdge(eid).from, p2*statesperprocess + fromP2), false);
		conflict.push(l2);
		assert(value(l2)==l_False);

		Lit l3 = ~mkLit(getNodeAPVar(g_under->getEdge(eid).to, p*statesperprocess + toP), false);
		conflict.push(l3);
		assert(value(l3)==l_False);
		Lit l4 = ~mkLit(getNodeAPVar(g_under->getEdge(eid).to, p2*statesperprocess + toP2), false);
		conflict.push(l4);
		assert(value(l4)==l_False);

		Lit l = ~mkLit(getTransitionVar(eid), false);
		conflict.push(l);
		assert(value(l)==l_False);

		return;
	}


	/*
	 * Do Symmetry Reduction. Check if any of the symmetry constraints are violated, and if so, build a clause describing this conflict.
	 * Calling this function might yield a conflict, but it does not have to if there is none.
	 *
	 * 	opt_ctl_symmetry Modes 1-3
	 *
	 */
	void checkStateLabelSymmetryConstraints(vec<Lit> & conflict, int startNode) {
		for (int i = 0; i < g_over->states(); i++) {
			for (int j = i+1; j < g_over->states(); j++) {
				if (i != startNode && j != startNode) {
					if(opt_verb>1) {
						printf("SYMMETRY: checking %d > %d\n", i, j);
					}
					// Reduce on edges in case of state label equivalence
					// Mode 1: Require label to be exactly the same
					// Mode 2: Require only that it is possible that the labels are the same, and if not, that there is a label symmetry violation
					if (opt_ctl_symmetry_statelabelandedges > 0 && g_under->statelabel[i]->Equiv( *g_over->statelabel[j] )
							&& (opt_ctl_symmetry_statelabelandedges == 2 || g_under->statelabel[j]->Equiv( *g_over->statelabel[i]))){
						if (g_under->nIncidentEnabled(i) > g_over->nIncidentEnabled(j)) {
							if(opt_verb>1) {
								printf("SYMMETRY: %d and %d have possibly the same state label and %d has %d > %d edges\n", i, j, i, g_under->nIncidentEnabled(i), g_over->nIncidentEnabled(j));
							}
							tmpConflict.clear();
							learnClauseSymmetryConflictEquivStateLabels(tmpConflict, i, j);
							if (conflict.size() == 0 || conflict.size() > tmpConflict.size()) {
								if(opt_verb>1) {
									printf("SYMMETRY: The currently smallest symmetry conflict is with states %d and %d\n", i, j);
								}
								tmpConflict.copyTo(conflict);
								if (opt_ctl_symmetry.operator int() == 1) // Mode 1: return first conflict instead of looking for smallest
									return;
							}
						}
					}
					// Reduce if i has a larger state label than j (for all extensions) despite i<j
					if (g_under->statelabel[i]->GreaterThan( *g_over->statelabel[j] )) {
						if(opt_verb>1) {
							printf("SYMMETRY: %d has a higher state label than %d\n", i, j);
						}
						tmpConflict.clear();
						learnClauseSymmetryConflict(tmpConflict, i, j);
						if (conflict.size() == 0 || conflict.size() > tmpConflict.size()) {
							if(opt_verb>1) {
								printf("SYMMETRY: The currently smallest symmetry conflict is with states %d and %d\n", i, j);
							}
							tmpConflict.copyTo(conflict);
							if (opt_ctl_symmetry.operator int() == 1) // Mode 1: return first conflict instead of looking for smallest
								return;
						}
					}
				}
			}
		}
	}

	// opt_ctl_symmetry Modes 4-6
	void checkEdgeSymmetryConstraints(vec<Lit> & conflict, int startNode) {
		for (int i = 0; i < g_over->states(); i++) {
			for (int j = i+1; j < g_over->states(); j++) {
				if (i != startNode && j != startNode) {
					if(opt_verb>1)
						printf("SYMMETRY: checking if edges(%d) = %d > %d = edges(%d)\n", i, g_under->nIncidentEnabled(i), g_over->nIncidentEnabled(j), j);
					if (g_under->nIncidentEnabled(i) > g_over->nIncidentEnabled(j)) { // We have a conflict
						if(opt_verb>1) {
							printf("SYMMETRY: We have a conflict, with states %d (%d edges) and %d (%d edges)\n", i, g_under->nIncidentEnabled(i), j, g_over->nIncidentEnabled(j));
						}
						tmpConflict.clear();
						learnClauseSymmetryConflictEdges(tmpConflict, i, j);
						if (conflict.size() == 0 || conflict.size() > tmpConflict.size()) {
							if(opt_verb>1) {
								printf("SYMMETRY: The currently smallest symmetry conflict is with states %d (%d edges) and %d (%d edges)\n", i, g_under->nIncidentEnabled(i), j, g_over->nIncidentEnabled(j));
							}
							tmpConflict.copyTo(conflict);
							if (opt_ctl_symmetry.operator int() == 4) // Mode 4: return first conflict instead of looking for smallest
								return;
						}
					}
				}
			}
		}
	}

	// Build a clause that describes a symmetry conflict with respect to state labels. We assume that the state label for b is
	// under the current assignment smaller than the state label for a, despite b>a. Learn a conflict to change this.
	// We start at the most significant bit and work our way down from there. We basically match the longest matching prefix.
	// For everything in the prefix, we learn that the stateLabel has to become 1 in b's overapproximation (if it is 0 currently)
	// or 0 in a's underapproximation (if it is 1 currently).
	// The bits after the longest matching prefix HAVE to be 1 in a's underapprox and 0 in b's overapprox. We add to the clause,
	// that they both could change.
	// Example:                        2>1 but state(1)= 01101 > state(2)= 00101. Note that we have Least significant bit first!
	// We learn one of the following bits must change =   ^^ ^              ^ ^
	void learnClauseSymmetryConflict(vec<Lit> & conflict, int a, int b) { // b>a, but over_APsToInt(b) < under_APstoInt(a).
		int i = g_over->apcount;
		do {
			i--;
			if (g_under->isAPinStateLabel(a, i)) {
				if(opt_verb>1) {
					printf("SYMMETRY: label(%d)>label(%d): learning that AP %d in state %d could be false\n", b, a, i, a);
				}
				Lit l = ~mkLit(getNodeAPVar(a, i), false);
				conflict.push(l);
				assert(value(l)==l_False);
			}
			if (!g_over->isAPinStateLabel(b, i)) {
				if(opt_verb>1) {
					printf("SYMMETRY: label(%d)>label(%d): learning that AP %d in state %d could be true\n", b, a, i, b);
				}
				Lit l = ~mkLit(getNodeAPVar(b, i), true);
				conflict.push(l);
				assert(value(l)==l_False);
			}
		} while ((!g_under->isAPinStateLabel(a, i) || g_over->isAPinStateLabel(b, i)) && i != 0); // stop when under_AP(a, i)=1 and over_AP(b, i)=0
		// Assertion not valid anymore, because we use this code from learnClauseSymmetryConflictEquivStateLabels, where it is not the case
		//assert(g_under->isAPinStateLabel(a, i) && !g_over->isAPinStateLabel(b, i)); // there must be something that differentiates them, otherwise a and b have the same statelabel!
	}

	// i has more edges than j, despite i<j. Learn a conflict describing this
	void learnClauseSymmetryConflictEdges(vec<Lit> & conflict, int i, int j) {
		DynamicGraph<int>::Edge e;
		int from, to;
		for (int k = 0; k < g_under->nIncident(i); k++) {
			e = g_under->incident(i, k);
			to = g_under->getEdge(e.id).to;
			if (g_under->edgeEnabled(e.id)) {
				if(opt_verb>1) {
					printf("SYMMETRY-Edges: edges(%d)>edges(%d): learning that edge %d from state %d could be false\n", i, j, e.id, i);
				}
				Lit l = ~mkLit(getTransitionVar(e.id), false);
				conflict.push(l);
				assert(value(l)==l_False);
			}
		}
		// Enable j's edges
		for (int k = 0; k < g_over->nIncident(j); k++) {
			e = g_over->incident(j, k);
			to = g_over->getEdge(e.id).to;
			if (!g_over->edgeEnabled(e.id)) {
				if(opt_verb>1) {
					printf("SYMMETRY-Edges: states(%d)>states(%d): learning that edge %d from state %d could be true\n", i, j, e.id, j);
				}
				Lit l = ~mkLit(getTransitionVar(e.id), true);
				conflict.push(l);
				assert(value(l)==l_False);
			}
		}


	}
	/*
	 * Called when statelabels are equal and number of outgoing edges of i is greater than in j.
	 */
	void learnClauseSymmetryConflictEquivStateLabels(vec<Lit> & conflict, int i, int j) { // j>i, but over_numedges(j) < under_numedges(i) and statelabels are equivalent.
		// The following edge-related literals are for the extensions of the current assignment, where both state labels are equivalent.
		learnClauseSymmetryConflictEdges(conflict, i, j);
		// Besides the possibility of adding/removing edges, there is also the possibility that the state label can be non-equivalent.
		// But, because we have g_under->nIncidentEnabled(i) > g_over->nIncidentEnabled(j), the only way for it to be non-equivalent is
		// when the statelabel of i is greater than the state label of j. In this case, we have a state label symmetry violation.
		learnClauseSymmetryConflict(conflict, i, j);

	}

	// Learn the trivial clause. Used for debugging purposes or to check if our clause learning algorithms are correct
	void learnNaiveClause(vec<Lit> & conflict, int startNode) {
		int w;
		for (int v = 0; v < vars.size(); v++) {
			if(value(v)!=l_Undef){
				Lit l = ~mkLit(v,value(v)==l_False);
				assert(value(l)==l_False);
				conflict.push(l);
			}
		}
	}

	// Clause learning for when a CTL formula is supposed to be true, but is false in the overapproximation
    // Starting from some initial state startNode (on the most toplevel call this will be the KripkeStructure's initial state,
	// but not necessarily so on recursive calls).
	void learnClausePos(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		 // only for debugging, you may uncomment this:

		/*
		if(opt_verb>1)
			printf("Learning naive clause...\n");
		learnNaiveClause(conflict, initialNode);
		return;
		*/


		if(opt_verb>1) {
			printf("Clause learning subformula... ");
			printFormula(&subf);
			printf(", for startNode %d\n", startNode);
		}

		if (subf.fairnessConstraints.size() > 0 && subf.op != EG) {
			if(opt_verb>1)
				printf("Learning naive clause, since we don't support clause learning with fairness constraints on this operator yet...\n");
			learnNaiveClause(conflict, initialNode);
			return;
		}

		if (subf.op == AND) {
			if(opt_verb>1)
				printf("Clause learning case AND...\n");

			learnAND(conflict, *subf.operand1, *subf.operand2, startNode);
		} else if (subf.op == OR) {
			if(opt_verb>1)
				printf("Clause learning case OR...\n");

			learnOR(conflict, *subf.operand1, *subf.operand2, startNode);
		} else if (subf.op == NEG) { // We push negation all the way to the atomic propositions, such that we only have literals and otherwise negation free formulas
			if(opt_verb>1)
				printf("Clause learning case NEG...\n");

			if (subf.operand1->op == NEG) {
				learnClausePos(conflict, *subf.operand1->operand1, startNode); // Be clever about double negation
				return;
			}
			if (subf.operand1->op == True) {
				// don't learn anything...!
				return;
			}
			CTLFormula inner1 {NEG, subf.operand1->operand1, NULL, 0};
			// inner2 is only needed for some connectives and thus declared below

			if (subf.operand1->op == AND) {
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer {OR, &inner1, &inner2, 0};

				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == OR) {
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer {AND, &inner1, &inner2, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == AX) {
				CTLFormula outer {EX, &inner1, NULL, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == EX) {
				CTLFormula outer {AX, &inner1, NULL, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == AG) {
				CTLFormula outer {EF, &inner1, NULL, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == EG) {
				CTLFormula outer {AF, &inner1, NULL, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == EF) {
				CTLFormula outer {AG, &inner1, NULL, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == AF) {
				CTLFormula outer {EG, &inner1, NULL, 0};
				learnClausePos(conflict, outer, startNode);
			} else if (subf.operand1->op == ID) {
				int p = subf.operand1->value;

				assert( g_over->isAPinStateLabel(startNode, p) ); // Assert that the startnode we are looking at does satisfy p, otherwise we had no conflict to learn
				Lit l = ~mkLit(getNodeAPVar(startNode, p), false);
				conflict.push(l);
				assert(value(l)==l_False);
			} else if (subf.operand1->op == AU) { // ~(Ï†_1 AU Ï†_2) â‰¡ ~Ï†_2 EW ( ~Ï†_1 ^ ~Ï†_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {EW, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else if (subf.operand1->op == AW) { // ~(Ï†_1 AW Ï†_2) â‰¡ ~Ï†_2 EU ( ~Ï†_1 ^ ~Ï†_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {EU, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else if (subf.operand1->op == EU) { // ~(Ï†_1 EU Ï†_2) â‰¡ ~Ï†_2 AW ( ~Ï†_1 ^ ~Ï†_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {AW, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else if (subf.operand1->op == EW) { // ~(Ï†_1 EW Ï†_2) â‰¡ ~Ï†_2 AU ( ~Ï†_1 ^ ~Ï†_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {AU, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else {
				assert (false); // we should not have missed anything
			}

		}
		else if (subf.op == ID) { // Not recursive
			if(opt_verb>1)
				printf("Clause learning case ID...\n");

			int p = subf.value;

			assert( !g_over->isAPinStateLabel(startNode, p) ); // Assert that the startnode we are looking at does not satisfy p, otherwise we had no conflict to learn
			Lit l = ~mkLit(getNodeAPVar(startNode, p), true);
			conflict.push(l);
			assert(value(l)==l_False);
		}
		else if (subf.op == True) { // Not recursive
			if(opt_verb>1)
				printf("Clause learning case True...\n");

			assert( false ); // We should never have to go here, since True can never cause a conflict
		}
		else if (subf.op == EX) {
			if(opt_verb>1)
				printf("Clause learning case EX...\n");

			learnEX(conflict, *subf.operand1, startNode);
		}
		else if (subf.op == AX) {
			if(opt_verb>1)
				printf("Clause learning case AX...\n");

			learnAX(conflict, *subf.operand1, startNode);
		}
		else if (subf.op == EG) {
			if(opt_verb>1)
				printf("Clause learning case EG...\n");

			if (subf.fairnessConstraints.size() > 0)
				learnEGFair(conflict, subf, startNode); // Note that we are calling it with subf, not subf.operand1
			else
				learnEG(conflict, *subf.operand1, startNode);
		}
		else if (subf.op == AG) {
			if(opt_verb>1)
				printf("Clause learning case AG...\n");

			learnAG(conflict, subf, startNode);
		}
		else if (subf.op == EF) {
			if(opt_verb>1)
				printf("Clause learning case EF...\n");

			learnEF(conflict, *subf.operand1, startNode);
		}
		else if (subf.op == AF) {
			if(opt_verb>1)
				printf("Clause learning case AF...\n");

			learnAF(conflict, *subf.operand1, startNode);
		}
		else if (subf.op == EW) {
			if(opt_verb>1)
				printf("Clause learning case EW...\n");

			learnEW(conflict, *subf.operand1, *subf.operand2, startNode);
		}
		else if (subf.op == EU) {
			if(opt_verb>1)
				printf("Clause learning case EU...\n");

			learnEU(conflict, *subf.operand1, *subf.operand2, startNode);
		}
		else if (subf.op == AW) {
			if(opt_verb>1)
				printf("Clause learning case AW...\n");

			learnAW(conflict, *subf.operand1, *subf.operand2, startNode, subf); // we give the entire formula as well, even though it's redundant, for convenience
		}
		else if (subf.op == AU) {
			if(opt_verb>1)
				printf("Clause learning case AU...\n");

			learnAU(conflict, *subf.operand1, *subf.operand2, startNode, subf); // we give the entire formula as well, even though it's redundant, for convenience
		}
		else {
			if(opt_verb>1)
				printf("Clause learning case OTHER...\n");

			for (int v = 0; v < vars.size(); v++) {
				if(value(v)!=l_Undef){
					Lit l = ~mkLit(v,value(v)==l_False);
					assert(value(l)==l_False);
					conflict.push(l);
					assert(value(l)==l_False);
				}
			}
		}
	}

	// Fallback: learn full clause
	void learnFullClause(vec<Lit> & conflict) {
		for (int v = 0; v < vars.size(); v++) {
			if(value(v)!=l_Undef){
				Lit l = ~mkLit(v,value(v)==l_False);
				assert(value(l)==l_False);
				conflict.push(l);
			}
		}
	}

	// We describe the individual clause learning strategies referring to phi, which is the inner formula.
	// The to-learn clause representing the inner formula is discovered recursively

	// At least one of the two parts of AND must be false, and we use the one which is false to build a clause
	void learnAND(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) {
		int phi1_over = ctl_over->solveFormula(subf1);
		int phi2_over = ctl_over->solveFormula(subf2);

		if(opt_verb>1) {
			printf("learnAND: phi1_over: "); ctl_standalone_over->printStateSet(*ctl_over->bitsets[phi1_over]);
			printf("learnAND: phi2_over: "); ctl_standalone_over->printStateSet(*ctl_over->bitsets[phi2_over]);
		}


		/* STRATEGY 1: Pick one subformula that works, and work with it to build a clause
		if (!phi1_over->operator [](startNode) && !phi1_under->operator [](startNode)) {
			learnClausePos(conflict, subf1, startNode);
		} else {
			assert(!phi2_over->operator [](startNode) && !phi2_under->operator [](startNode)); // If this is violated, that means that both parts of the AND formula are satisfied -- then there should not be a conflict
			learnClausePos(conflict, subf2, startNode);
		}
		*/

		/* STRATEGY 2: Pick both subformulas if they both work, and compare which has the smallest clause */
		if (!ctl_over->bitsets[phi1_over]->operator [](startNode) && !ctl_over->bitsets[phi2_over]->operator [](startNode)) {
			// Both subformulas conflict and are suitable to learn a clause from. Build both clauses and learn the smaller one.
			vec<Lit> conflict2;
			conflict.copyTo(conflict2);
			learnClausePos(conflict,  subf1, startNode);
			learnClausePos(conflict2, subf2, startNode);
			if (conflict.size() > conflict2.size()) {
				conflict2.copyTo(conflict);
			}
		} else if (!ctl_over->bitsets[phi1_over]->operator [](startNode)) {
			learnClausePos(conflict, subf1, startNode);
		} else {
			assert(!ctl_over->bitsets[phi2_over]->operator [](startNode)); // If this is violated, that means that both parts of the AND formula are satisfied -- then there should not be a conflict
			learnClausePos(conflict, subf2, startNode);
		}
		ctl_over->freeBitset(phi1_over);
		ctl_over->freeBitset(phi2_over);
	}

	// Both parts of the OR must be false, and we OR together their recursively learned sub-clauses
	void learnOR(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) {
		learnClausePos(conflict, subf1, startNode);
		learnClausePos(conflict, subf2, startNode);
	}

	// At least one neighbour gets phi, or at least one disabled edge to a neighbour is enabled
	void learnEX(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		DynamicGraph<int>::Edge e;
		int from, to;
		for (int i = 0; i < g_over->nIncident(startNode); i++) {
			e = g_over->incident(startNode, i);
			from = g_over->getEdge(e.id).from;
			to = g_over->getEdge(e.id).to;
			if(opt_verb>1)
				printf("learnEX: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);

			// Two cases: Either the edge is enabled, then we need to add to the clause the possibility that the
			// successor (the "to" state) satisfies p. Or the edge is disabled and we add the possibility that the
			// edge is enabled.
			// The first case is recursively solved
			if (g_over->edgeEnabled(e.id)) {
				learnClausePos(conflict, subf, to);
			} else { // No recursive evaluation needed
				Lit l = ~mkLit(getTransitionVar(e.id), true);
				conflict.push(l);
				assert(value(l)==l_False);
			}
		}
	}

	// At least one neighbour enables phi, or at least one enabled edge to a neighbour is disabled
	// In order to make shorter clauses, we will only learn one clause: find an enabled edge to a state not satisfying
	//
	void learnAX(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		//long bitsets = Bitset::remainingBitsets();
		int phi_over = ctl_over->solveFormula(subf,opt_ctl_learn_cache);
		//long leaked = Bitset::remainingBitsets()-bitsets;
		//if(leaked!=2){
		//	throw std::runtime_error("Leaked bitsets during solve call in learnAX.1!");
		//}

		DynamicGraph<int>::Edge e;
		int from, to;
		for (int i = 0; i < g_over->nIncident(startNode); i++) {
			e = g_over->incident(startNode, i);
			from = g_over->getEdge(e.id).from;
			to = g_over->getEdge(e.id).to;
			if(opt_verb>1)
				printf("learnAX: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);

			// We are looking for one single enabled edge such that the destination does not satisfy phi, and require
			// that either it will satisfy phi or that the edge be disabled
			if (g_under->edgeEnabled(e.id) && ! ctl_over->bitsets[phi_over]->operator [](to)) {
				learnClausePos(conflict, subf, to);
				Lit l = ~mkLit(getTransitionVar(e.id), false);
				conflict.push(l);
				assert(value(l)==l_False);
				ctl_over->freeBitset(phi_over);
				return;
			}
		}
		assert(false); // No clause was learned, which means that there is no enabled edge such that phi does not hold in the destination -- a contradiction to the fact, that AX phi does not hold!

	}


	// Successor to "phi-reachable" state satisfies phi or enable transition from any "phi-reachable" state
	//
	void learnEG(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		int phi_over = ctl_over->solveFormula(subf,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		//ctl_over->printStateSet(*phi);
		//printFormula2(subf);

		DynamicGraph<int>::Edge e;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!ctl_over->bitsets[phi_over]->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnEG: initial state does not satisfy phi\n");
			learnClausePos(conflict, subf, startNode);
			ctl_over->freeBitset(phi_over);
			return;
		}

		std::queue <int> list; // this queue denotes all the states satisfying phi, whose neighbours have to be inspected
		list.push(startNode);
		int visited = ctl_over->getFreshBitset();//new Bitset(g_over->states()); // This bitset denotes all the visited nodes
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			if(opt_verb>1)
				printf("learnEG: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEG: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id)) {
					if (ctl_over->bitsets[phi_over]->operator [](to) && !ctl_over->bitsets[visited]->operator [](to)) {
						list.push(to);

					} else if (!ctl_over->bitsets[phi_over]->operator [](to)) {
						learnClausePos(conflict, subf, to);
					} // else we already have visited the phi-state "to" and recursively all its neighbours
				} else {
					Lit l = ~mkLit(getTransitionVar(e.id), true);
					conflict.push(l);
					assert(value(l)==l_False);
				}
			}
			ctl_over->bitsets[visited]->set(from);
		}
		ctl_over->freeBitset(visited);
		ctl_over->freeBitset(phi_over);
	}

	// Start out by learning EG. Then, for each fairness constraints, learn that it can switch to being true in every state
	// This is not ideal, ideally we'd restrict it to reachable states.
	void learnEGFair(vec<Lit> & conflict, CTLFormula &thisf, int startNode) {
		learnEG(conflict, *thisf.operand1, startNode);
		int cSet;
		for (int c = 0; c < thisf.fairnessConstraints.size(); c ++) {
			cSet = ctl_over->solveFormula(*thisf.fairnessConstraints[c],opt_ctl_learn_cache);
			for (int i = 0; i < g_over->states(); i++) {
				if (!ctl_over->bitsets[cSet]->operator [](i)) {
					learnClausePos(conflict, *thisf.fairnessConstraints[c], i);
				}
			}
			ctl_over->freeBitset(cSet);
		}
	}

	// We know there currently exists a fair path from startNode.
	// We learn that this path could become unfair, in which case startNode possibly is not a witness anymore for some outer formula (probably AG).
	// Currently we learn that any state satisfying a fairness constraint may stop to satisfy it. Ideally we'd limit this to reachable states
	void learnMakePathUnfair(vec<Lit> & conflict, CTLFormula &thisf, int startNode){
		int cSet;
		for (int c = 0; c < thisf.fairnessConstraints.size(); c ++) {
			CTLFormula* notc = newCTLFormula();
			notc->op = NEG;
			notc->operand1 = thisf.fairnessConstraints[c];
			cSet = ctl_under->solveFormula(*thisf.fairnessConstraints[c],opt_ctl_learn_cache);
			for (int i = 0; i < g_over->states(); i++) {
				if (ctl_under->bitsets[cSet]->operator [](i)) {
					learnClausePos(conflict, *notc, i);
				}
			}
			ctl_under->freeBitset(cSet);
		}
	}

	// Find a path to a not-phi state, either that state must satisfy phi, or one of the edges must be disabled
	void learnAG(vec<Lit> & conflict, CTLFormula &thisf, int startNode) {
		CTLFormula* subf = thisf.operand1;
		// Solve the inner subformula of the entire formula
		int phi_over = ctl_over->solveFormula(*subf,opt_ctl_learn_cache);
		if(opt_verb>1) {
			printf("learnAG: phi_over "); ctl_standalone_over->printStateSet(*ctl_over->bitsets[phi_over]);
		}

		// This set represents the set of states satisfying EG_C True, where C are the fairness constraint from thisf (AG_C phi).
		// If we don't have any fairness constraints ("else"), we simply use a bitset which is true for every state
		int fair_over;
		if (thisf.fairnessConstraints.size() > 0) {
			CTLFormula* fairFormula = newCTLFormula();
			CTLFormula* fairFormulaTrue = newCTLFormula();
			fairFormulaTrue->op = True;
			fairFormula->op = EG;
			fairFormula->operand1 = fairFormulaTrue;
			fairFormula->fairnessConstraints = thisf.fairnessConstraints;
			fair_over = ctl_over->solveFormula(*fairFormula,opt_ctl_learn_cache);
		} else {
			fair_over = ctl_over->getDirtyBitset();
			ctl_over->bitsets[fair_over]->memset(true);
		}

		DynamicGraph<int>::Edge e;
		int eid = e.id;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!ctl_over->bitsets[phi_over]->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnAG: initial state does not satisfy phi\n");
			assert(ctl_over->bitsets[fair_over]->operator [](startNode)); // We assert there exists a fair path. If there is not, AG would trivially hold
			// We learn that either the inner formula has to hold in the initial state, or we make sure no fair path exists
			learnClausePos(conflict, *subf, startNode);
			if (thisf.fairnessConstraints.size() > 0) {
				learnMakePathUnfair(conflict, thisf, startNode);
			}
			ctl_over->freeBitset(phi_over);
			ctl_over->freeBitset(fair_over);
			return;
		}

		std::queue <int> s;
		s.push(startNode);
		std::map <int,int> parent; // from this we retreive the path of edges from start state to some state that doesn't satisfy phi
		int visited = ctl_over->getFreshBitset(); // This bitset denotes all the visited nodes
		bool done = false;
		while (s.size() > 0 && !done) { // stack of visited elements
			from = s.front();
			s.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnAG: Considering state %d\n", from);

			for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAG: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid)) {
					if(opt_verb>1)
						printf("learnAG: Adding map parent(%d) = %d. phi_over(to): %d, visited: %d\n", to, from, ctl_over->bitsets[phi_over]->operator [](to), ctl_over->bitsets[visited]->operator [](to));
					// The state does not satisfy phi and has a path satisfying the fairness constraints.
					if (!ctl_over->bitsets[phi_over]->operator [](to) && ctl_over->bitsets[fair_over]->operator [](to)) {
						parent[to] = from;
						if(opt_verb>1)
							printf("learnAG: Found state no %d, which does not satisfy phi. edgeid: %d, from: %d, to: %d\n", to, e.id, from, to);
						// We learn that either the inner formula has to hold in this state, or we make sure no fair path exists from this state
						learnClausePos(conflict, *subf, to);
						if (thisf.fairnessConstraints.size() > 0) {
							learnMakePathUnfair(conflict, thisf, to);
						}
						done = true; // we have found a state that does not satisfy phi, exit loop and retreive path
					} else if (!ctl_over->bitsets[visited]->operator [](to)) {
						parent[to] = from;
						if(opt_verb>1)
							printf("learnAG: putting %d in queue\n", to);
						s.push(to);
					}
				}
			}
		}
		assert(done);

		if(opt_verb>1)
			printf("learnAG: reconstructing path from %d to %d\n", startNode, to);
		while (to != startNode) {
			int to1 = g_under->getEdge(eid).to;
			int from1 = g_under->getEdge(eid).from;
			if(opt_verb>1)
				printf("learnAG: %d -> %d\n", from1, to1);

			Lit l = ~mkLit(getTransitionVar(eid), false);
			conflict.push(l);
			assert(value(l)==l_False);

			to = from;
			from = parent[to];
			eid = g_under->getEdge(from, to);
		}
		ctl_over->freeBitset(phi_over);
		ctl_over->freeBitset(visited);
		ctl_over->freeBitset(fair_over);
	}

	// Find all reachable states, at least one of them should satisfy phi OR there should be a transition enabled
	// from a reachable state to an unreachable state
	void learnEF(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		DynamicGraph<int>::Edge e;
		int from, to;


		std::queue <int> list; // this queue denotes the current path of visited states, whose neighbours have to be inspected
		list.push(startNode);
		int visited = ctl_over->getFreshBitset();
		// We first employ BFS to find all reachable nodes. We add them to the clause
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnEF: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEF: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id) && !ctl_over->bitsets[visited]->operator [](to)) {
					list.push(to);
				}
			}
			learnClausePos(conflict, subf, from); // learn that this reachable node satisfy phi
		}
		// Now for part two of the clause:
		// iterate over edges of visited states to unvisited states; add literal to enable these edges
		for (int i = 0; i < ctl_over->bitsets[visited]->size(); i++) {
			if (ctl_over->bitsets[visited]->operator [](i)) {
				for (int j = 0; j < g_over->nIncident(i); j++) {
					e = g_over->incident(i, j);
					to = g_over->getEdge(e.id).to;

					if (!ctl_over->bitsets[visited]->operator [](to)) {
						if(opt_verb>1)
							printf("learnEF: Learning that the edge between %d and %d (edgeid: %d) could be enabled.\n", i, to, e.id);

						assert(!g_over->edgeEnabled(e.id)); // if edge was enabled, then j would be reachable, which means j would be visited
						Lit l = ~mkLit(getTransitionVar(e.id), true);
						conflict.push(l);
						assert(value(l)==l_False);
					}
				}
			}
		}
		ctl_over->freeBitset(visited);
	}

	// Find lasso of states that satisfy not phi, at least one of them should satisfy phi OR one of the edges of the lasso should become disabled
	void learnAF(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		int phi_over = ctl_over->solveFormula(subf); // Solve the inner subformula of the entire formula


		DynamicGraph<int>::Edge e;
		int from, to, eid, pred, predpred, to1, from1;

		std::stack <int> s;
		s.push(startNode);
		std::map <int,int> parent; // from this we retreive the path of edges from start state to some state that doesn't satisfy phi
		int visited = ctl_over->getFreshBitset();
		ctl_over->bitsets[visited]->set(startNode);
		bool done = false;
		while (s.size() > 0 && !done) { // stack of visited elements
			from = s.top();
			s.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnAG: Considering state %d\n", from);

			for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAF: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !ctl_over->bitsets[phi_over]->operator [](to)) {
					if (!ctl_over->bitsets[visited]->operator [](to)) { // explore new state in graph
					if(opt_verb>1)
						printf("learnAF: Adding map parent(%d) = %d. phi_over(to): %d, visited: %d. putting %d in queue\n", to, from, ctl_over->bitsets[phi_over]->operator [](to),  ctl_over->bitsets[visited]->operator [](to), to);

						s.push(to);
						parent[to] = from;
					} else { // check if this next state is part of the current path, in this case we have found a lasso
						if(opt_verb>1)
							printf("learnAF: Checking if %d makes a loop in current path\n", to);

						pred = to;
						predpred = from;
						if (predpred == to) {
							done = true; // We have found a lasso
							if(opt_verb>1)
								printf("learnAF: %d does make a loop in current path (via %d)\n", to, from);
						}
						while (predpred != startNode && !done) {
							to1 = g_under->getEdge(eid).to;
							from1 = g_under->getEdge(eid).from;
							if(opt_verb>1)
								printf("learnAF: %d -> %d\n", from1, to1);

							Lit l = ~mkLit(getTransitionVar(eid), false);
							assert(value(l)==l_False);

							pred = predpred;
							predpred = parent[predpred];
							eid = g_under->getEdge(predpred, pred);

							if (predpred == to) {
								done = true; // We have found a lasso
								if(opt_verb>1)
									printf("learnAF: %d does make a loop in current path (via %d)\n", to, from);
							}
						}
					}
				}
			}
		}

		assert(done);

		if(opt_verb>1)
			printf("learnAF: reconstructing path from %d to %d\n", startNode, to);

		// First, get the very last edge in the path accounted for (the one that completes the loop of the lasso)
		eid = g_under->getEdge(from, to);
		Lit l = ~mkLit(getTransitionVar(eid), false);
		conflict.push(l);
		assert(value(l)==l_False);
		learnClausePos(conflict, subf, from); // learn that this reachable node satisfy phi
		to = from;
		from = parent[to];
		eid = g_under->getEdge(from, to);

		while (to != startNode) {
			int to1 = g_under->getEdge(eid).to;
			int from1 = g_under->getEdge(eid).from;
			if(opt_verb>1)
				printf("learnAF: %d -> %d\n", from1, to1);

			// Learn that the edge might be turned off
			Lit l = ~mkLit(getTransitionVar(eid), false);
			conflict.push(l);
			assert(value(l)==l_False);

			learnClausePos(conflict, subf, from); // learn that this reachable node satisfy phi

			to = from;
			from = parent[to];
			eid = g_under->getEdge(from, to);
		}
		ctl_over->freeBitset(visited);
		ctl_over->freeBitset(phi_over);
	}



	// Find all phi-reachable states, at least one of them should satisfy psi OR there should be a transition enabled
	// from a phi-reachable state (to anywhere) OR a non-phi-reachable state that is a successor to a phi-reachable state satisfies phi, or psi

	// For EU this should be pretty much the same, except that we only add transitions from phi-reachable to non-phi-reachable
	// FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME Not done yet
	void learnEW(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) { // phi EW psi
		int phi = ctl_over->solveFormula(subf1,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int psi = ctl_over->solveFormula(subf2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula


		DynamicGraph<int>::Edge e;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi or psi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!ctl_over->bitsets[phi]->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnEW: initial state does not satisfy phi or psi\n");
			assert(!ctl_over->bitsets[psi]->operator [](startNode)); // If psi was satisfied in initial state, then we would not have a conflict
			learnClausePos(conflict, subf1, startNode); // Either phi must be true...
			learnClausePos(conflict, subf2, startNode); // ... or psi must be true
			ctl_over->freeBitset(phi);
			ctl_over->freeBitset(psi);
			return;
		}


		std::queue <int> list; // this queue denotes the current path of visited states, whose neighbours have to be inspected
		list.push(startNode);
		int visited = ctl_over->getFreshBitset();
		// We first employ BFS to find all phi-reachable nodes. We add them to the clause
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnEW: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEW: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id) && !ctl_over->bitsets[visited]->operator [](to) && ctl_over->bitsets[phi]->operator [](to)) {
					list.push(to);
				} else if (g_over->edgeEnabled(e.id) && !ctl_over->bitsets[visited]->operator [](to) && !ctl_over->bitsets[phi]->operator [](to)) {
					learnClausePos(conflict, subf1, to); // learn that this successor to a phi-reachable state satisfy phi
					learnClausePos(conflict, subf2, to); // learn that this successor to a phi-reachable state satisfy psi
				}
			}
			assert(!ctl_over->bitsets[psi]->operator [](startNode)); // If psi was satisfied in this phi-reachable state, then we would not have a conflict
			learnClausePos(conflict, subf2, from); // learn that this phi-reachable node satisfy psi
		}
		// Now for part two of the clause:
		// iterate over edges from visited states; add literal to enable these edges
		for (int i = 0; i < ctl_over->bitsets[visited]->size(); i++) {
			if (ctl_over->bitsets[visited]->operator [](i)) {
				for (int j = 0; j < g_over->nIncident(i); j++) {
					e = g_over->incident(i, j);
					to = g_over->getEdge(e.id).to;

					if (!g_over->edgeEnabled(e.id)) {
						if(opt_verb>1)
							printf("learnEW: Learning that the edge between %d and %d (edgeid: %d) could be enabled.\n", i, to, e.id);
						Lit l = ~mkLit(getTransitionVar(e.id), true);
						conflict.push(l);
						assert(value(l)==l_False);
					}
				}
			}
		}
		ctl_over->freeBitset(phi);
		ctl_over->freeBitset(psi);
		ctl_over->freeBitset(visited);
	}


	// Find all phi-reachable states, at least one of them should satisfy psi OR there should be a transition enabled
	// from a phi-reachable state (to anywhere) OR a non-phi-reachable state that is a successor to a phi-reachable state satisfies phi, or psi

	// For EU this should be pretty much the same, except that we only add transitions from phi-reachable to non-phi-reachable
	void learnEU(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) { // phi EW psi
		int phi = ctl_over->solveFormula(subf1,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int psi = ctl_over->solveFormula(subf2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula


		DynamicGraph<int>::Edge e;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi or psi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!ctl_over->bitsets[phi]->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnEU: initial state does not satisfy phi or psi\n");
			assert(!ctl_over->bitsets[psi]->operator [](startNode)); // If psi was satisfied in initial state, then we would not have a conflict
			learnClausePos(conflict, subf1, startNode); // Either phi must be true...
			learnClausePos(conflict, subf2, startNode); // ... or psi must be true
			ctl_over->freeBitset(phi);
			ctl_over->freeBitset(psi);
			return;
		}


		std::queue <int> list; // this queue denotes the current path of visited states, whose neighbours have to be inspected
		list.push(startNode);
		int visited = ctl_over->getFreshBitset();
		// We first employ BFS to find all phi-reachable nodes. We add them to the clause
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnEU: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEU: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id) && !ctl_over->bitsets[visited]->operator [](to) && ctl_over->bitsets[phi]->operator [](to)) {
					list.push(to);
				} else if (g_over->edgeEnabled(e.id) && !ctl_over->bitsets[visited]->operator [](to) && !ctl_over->bitsets[phi]->operator [](to)) {
					learnClausePos(conflict, subf1, to); // learn that this successor to a phi-reachable state satisfy phi
					learnClausePos(conflict, subf2, to); // learn that this successor to a phi-reachable state satisfy psi
				}
			}
			assert(!ctl_over->bitsets[psi]->operator [](startNode)); // If psi was satisfied in this phi-reachable state, then we would not have a conflict
			learnClausePos(conflict, subf2, from); // learn that this phi-reachable node satisfy psi
		}
		// Now for part two of the clause:
		// iterate over edges from visited states; add literal to enable these edges
		for (int i = 0; i < ctl_over->bitsets[visited]->size(); i++) {
			if (ctl_over->bitsets[visited]->operator [](i)) {
				for (int j = 0; j < g_over->nIncident(i); j++) {
					e = g_over->incident(i, j);
					to = g_over->getEdge(e.id).to;

					if (!g_over->edgeEnabled(e.id) && !ctl_over->bitsets[visited]->operator [](to)) { // This is where the difference to EW is... we don't learn the edge if it merely leads to a phi-reachable state
						if(opt_verb>1)
							printf("learnEU: Learning that the edge between %d and %d (edgeid: %d) could be enabled.\n", i, to, e.id);
						Lit l = ~mkLit(getTransitionVar(e.id), true);
						conflict.push(l);
						assert(value(l)==l_False);
					}
				}
			}
		}
		ctl_over->freeBitset(phi);
		ctl_over->freeBitset(psi);
		ctl_over->freeBitset(visited);
	}
/*
 * NOTE: This was a previous attempt, where we actually tried finding a lasso, instead of just a finite path. It *does* work, but it
 * is inefficient since it will return the full clause unless AG EX True holds in the underapproximation.
 *
	// We take an additional fAll parameter, which contains the entire formula. I.e. fAll = subf1 AW subf2. This is not checked.
	void learnAW(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode, CTLFormula &fAll) {

		// We will only proceed to do efficient clause learning if AG EX True holds. Otherwise learn default clause.
		// The strategy in learnAND will compare both this clause to the clause learned from AG EX True and will prefer the smaller one.
		// So it's not a problem that we learn a full clause here, it will not end up being learned (rather something on AG EX True will be)
		Bitset* infinitePaths = ctl_under->solveFormula(*fAGEXTrue);
		if(!infinitePaths->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnAW: learning full clause, since AG EX True does not hold\n");
			learnFullClause(conflict);
			delete infinitePaths;
			return;
		}

		Bitset* phi1_under = ctl_under->solveFormula(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi1_over = ctl_over->solveFormula(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi2_under = ctl_under->solveFormula(subf2); // Solve the inner subformula of the entire formula
		Bitset* phi2_over = ctl_over->solveFormula(subf2); // Solve the inner subformula of the entire formula
		Bitset* fAll_under = ctl_under->solveFormula(fAll); // Solve the inner subformula of the entire formula
		Bitset* fAll_over = ctl_over->solveFormula(fAll); // Solve the inner subformula of the entire formula

		DynamicGraph::Edge e;
		int from, to, eid, pred, predpred, to1, from1;



		std::stack <int> s;
		s.push(startNode);
		std::map <int,int> parent; // from this we retreive the path of edges from start state to some state that doesn't satisfy phi
		Bitset *visited = new Bitset(g_over->states()); // This bitset denotes all the visited nodes
		visited->set(startNode);
		bool done = false;
		while (s.size() > 0 && !done) { // stack of visited elements
			from = s.top();
			s.pop();
			visited->set(from);
			if(opt_verb>1)
				printf("learnAW: Considering state %d\n", from);

			for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAW: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !fAll_under->operator [](to) && !fAll_over->operator [](to)) {
					if (!visited->operator [](to)) { // explore new state in graph
						if(opt_verb>1)
							printf("learnAW: Adding map parent(%d) = %d. phi1_under(to): %d, phi1_over(to): %d, phi2_under(to): %d, phi2_over(to): %d, visited: %d. putting %d in queue\n", to, from, phi1_under->operator [](to), phi1_over->operator [](to), phi2_under->operator [](to), phi2_over->operator [](to), visited->operator [](to), to);

						s.push(to);
						parent[to] = from;
					} else { // check if this next state is part of the current path, in this case we have found a lasso
						if(opt_verb>1)
							printf("learnAW: Checking if %d makes a loop in current path\n", to);

						pred = to;
						predpred = from;
						if (predpred == to) {
							done = true; // We have found a lasso
							if(opt_verb>1)
								printf("learnAW: %d does make a loop in current path (via %d)\n", to, from);
						}
						while (predpred != startNode && !done) {
							to1 = g_under->getEdge(eid).to;
							from1 = g_under->getEdge(eid).from;
							if(opt_verb>1)
								printf("learnAW: %d -> %d\n", from1, to1);

							Lit l = ~mkLit(getTransitionVar(eid), false);
							assert(value(l)==l_False);

							pred = predpred;
							predpred = parent[predpred];
							eid = g_under->getEdge(predpred, pred);

							if (predpred == to) {
								done = true; // We have found a lasso
								if(opt_verb>1)
									printf("learnAW: %d does make a loop in current path (via %d)\n", to, from);
							}
						}
					}
				}
			}
		}

		assert(done);

		if(opt_verb>1)
			printf("learnAW: reconstructing path from %d to %d\n", startNode, to);

		// First, get the very last edge in the path accounted for (the one that completes the loop of the lasso)
		eid = g_under->getEdge(from, to);

		std::stack <int> reverse; // we have to invert the path, since currently it goes from target to "start", and we have to start from "start"
		reverse.push(eid);

		to = from;
		from = parent[to];
		eid = g_under->getEdge(from, to);

		while (to != startNode) {
			reverse.push(eid);
			to = from;
			from = parent[to];
			eid = g_under->getEdge(from, to);
		}

		while (reverse.size() > 0) {
			eid = reverse.top();
			reverse.pop();
			to = g_under->getEdge(eid).to;
			from = g_under->getEdge(eid).from;
			learnClausePos(conflict, subf2, from); // learn that this reachable node satisfy psi
			if (!phi1_under->operator [](from)) {
				// We have reached a state where ~phi1 ^ ~ phi2. We can conclude by adding to the clause that phi1 v phi2 should hold in this state
				assert(!phi1_over->operator [](from));
				learnClausePos(conflict, subf1, from); // learn that this reachable node satisfy phi
				delete phi1_under;
				delete phi1_over;
				delete phi2_under;
				delete phi2_over;
				delete fAll_under;
				delete fAll_over;
				delete infinitePaths;
				return;
			} else {
				// Learn that we can disable the transition
				Lit l = ~mkLit(getTransitionVar(eid), false);
				conflict.push(l);
				assert(value(l)==l_False);
			}
		}
		assert(false); // We should always be able to find a concluding state where ~phi1 ^ ~ phi2 along our path.
	}
*/

	// Find a finite path such that all except the last state satisfy phi, and no state satisfies psi. Learn a clause where you OR together:
	// - remove an edge along the path
	// - make the last state satisfy phi
	// - make any state on the path satisfy psi
	//
    // In theory we'd have to find a lasso. However, AG EX True guarantees that even if we have a finite path, we will learn
	// clauses that extend it to an infinite path. So we don't actually need to make sure that the path we find is infinite, which
	// makes everything a lot simpler.
	//
	// We take an additional fAll parameter, which contains the entire formula. I.e. fAll = subf1 AW subf2. This is not checked.
	void learnAW(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode, CTLFormula &fAll) {
		int phi1_under = ctl_under->solveFormula(subf1,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int phi1_over = ctl_over->solveFormula(subf1,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int phi2_under = ctl_under->solveFormula(subf2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int phi2_over = ctl_over->solveFormula(subf2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int fAll_under = ctl_under->solveFormula(fAll,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int fAll_over = ctl_over->solveFormula(fAll,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula

		DynamicGraph<int>::Edge e;
		int from, to, eid, pred, predpred, to1, from1;

		std::stack <int> s;
		s.push(startNode);
		std::map <int,int> parent; // from this we retreive the path of edges from start state to some state that doesn't satisfy phi
		int visited = ctl_over->getFreshBitset();
		ctl_over->bitsets[visited]->set(startNode);
		bool done = false;
		while (s.size() > 0 && !done) { // stack of visited elements
			from = s.top();
			s.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnAW/AU: Considering state %d\n", from);

			if (!ctl_over->bitsets[phi1_over]->operator [](from) && !ctl_over->bitsets[phi2_over]->operator [](from)) {
				done = true;
			}
			else for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAW/AU: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !ctl_under->bitsets[fAll_under]->operator [](to) && !ctl_over->bitsets[fAll_over]->operator [](to)) {
					if (!ctl_over->bitsets[visited]->operator [](to)) { // explore new state in graph
						if(opt_verb>1)
							printf("learnAW/AU: Adding map parent(%d) = %d. phi1_under(to): %d, phi1_over(to): %d, phi2_under(to): %d, phi2_over(to): %d, visited: %d. putting %d in queue\n", to, from, ctl_under->bitsets[phi1_under]->operator [](to), ctl_over->bitsets[phi1_over]->operator [](to), ctl_under->bitsets[phi2_under]->operator [](to), ctl_over->bitsets[phi2_over]->operator [](to), ctl_over->bitsets[visited]->operator [](to), to);

						s.push(to);
						parent[to] = from;
					}
				}
			}
		}
		assert(done);

		if(opt_verb>1)
			printf("learnAW/AU: reconstructing path from %d to %d\n", startNode, from);

		if (from == startNode) {
			if(opt_verb>1)
				printf("learnAW/AU: phi and psi fail to hold in startNode\n");
			learnClausePos(conflict, subf1, from); // learn that this initial node satisfy phi
			learnClausePos(conflict, subf2, from); // learn that this initial node satisfy psi
			ctl_over->freeBitset(phi1_over);
			ctl_over->freeBitset(phi2_over);
			ctl_over->freeBitset(fAll_over);
			ctl_under->freeBitset(phi1_under);
			ctl_under->freeBitset(phi2_under);
			ctl_under->freeBitset(fAll_under);
			ctl_over->freeBitset(visited);
			return;
		}
		learnClausePos(conflict, subf1, from); // learn that this reachable node satisfy phi
		learnClausePos(conflict, subf2, from); // learn that this reachable node satisfy psi

		while (from != startNode) {
			to = from;
			from = parent[from];
			eid = g_under->getEdge(from, to);

			learnClausePos(conflict, subf2, from); // learn that this reachable node satisfy psi

			// Learn that we can disable the transition
			Lit l = ~mkLit(getTransitionVar(eid), false);
			conflict.push(l);
			assert(value(l)==l_False);
		}
		ctl_over->freeBitset(phi1_over);
		ctl_over->freeBitset(phi2_over);
		ctl_over->freeBitset(fAll_over);
		ctl_under->freeBitset(phi1_under);
		ctl_under->freeBitset(phi2_under);
		ctl_under->freeBitset(fAll_under);
		ctl_over->freeBitset(visited);
	}


	// Our approach is very similar to learnAW:
	// Find a finite path such that all except the last state satisfy phi, and no state satisfies psi. Learn a clause where you OR together:
	// - remove an edge along the path
	// - make the last state satisfy phi
	// - make any state on the path satisfy psi
	//
	// However, such a path may not exist. If no path like that exists, then we learnAF on psi.
	//
	// In theory we'd have to find a lasso. However, AG EX True guarantees that even if we have a finite path, we will learn
	// clauses that extend it to an infinite path. So we don't actually need to make sure that the path we find is infinite, which
	// makes everything a lot simpler.
	//
	// We take an additional fAll parameter, which contains the entire formula. I.e. fAll = subf1 AW subf2. This is not checked.
	void learnAU(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode, CTLFormula &fAll) {
		int phi1_under = ctl_under->solveFormula(subf1,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int phi1_over = ctl_over->solveFormula(subf1,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int phi2_under = ctl_under->solveFormula(subf2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int phi2_over = ctl_over->solveFormula(subf2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int fAll_under = ctl_under->solveFormula(fAll,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
		int fAll_over = ctl_over->solveFormula(fAll,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula

		DynamicGraph<int>::Edge e;
		int from, to, eid, pred, predpred, to1, from1;

		std::stack <int> s;
		s.push(startNode);
		std::map <int,int> parent; // from this we retreive the path of edges from start state to some state that doesn't satisfy phi
		int visited = ctl_over->getFreshBitset();
		ctl_over->bitsets[visited]->set(startNode);
		bool done = false;
		while (s.size() > 0 && !done) { // stack of visited elements
			from = s.top();
			s.pop();
			ctl_over->bitsets[visited]->set(from);
			if(opt_verb>1)
				printf("learnAU: Considering state %d\n", from);

			if (!ctl_over->bitsets[phi1_over]->operator [](from) && !ctl_over->bitsets[phi2_over]->operator [](from)) {
				done = true;
			}
			else for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAU: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !ctl_under->bitsets[fAll_under]->operator [](to) && !ctl_over->bitsets[fAll_over]->operator [](to)) {
					if (!ctl_over->bitsets[visited]->operator [](to)) { // explore new state in graph
						if(opt_verb>1)
							printf("learnAU: Adding map parent(%d) = %d. phi1_under(to): %d, phi1_over(to): %d, phi2_under(to): %d, phi2_over(to): %d, visited: %d. putting %d in queue\n", to, from, ctl_under->bitsets[phi1_under]->operator [](to), ctl_over->bitsets[phi1_over]->operator [](to), ctl_under->bitsets[phi2_under]->operator [](to), ctl_over->bitsets[phi2_over]->operator [](to), ctl_over->bitsets[visited]->operator [](to), to);

						s.push(to);
						parent[to] = from;
					}
				}
			}
		}
		if(!done) { // we know that p AW q holds (since we did not find any paths to ~p^~q), but p AU q does not (by assumption)
			if(opt_verb>1)
				printf("learnAU: No path exists, where eventually ~p^~q holds. This implies p AW q holds. Learn AF q\n");

			CTLFormula* fAFphi2 = newCTLFormula();
			fAFphi2->op = AF;
			fAFphi2->operand1 = &subf2;
			int AFphi2_under = ctl_under->solveFormula(*fAFphi2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula
			int AFphi2_over = ctl_over->solveFormula(*fAFphi2,opt_ctl_learn_cache); // Solve the inner subformula of the entire formula

			if(opt_verb>1) {
				printf("learnAU: These are the under and overapproximations for AF q: ");
				ctl_over->printStateSet(*ctl_under->bitsets[AFphi2_under]);
				ctl_over->printStateSet(*ctl_over->bitsets[AFphi2_over]);
				printf("\n");
			}

			assert(!ctl_over->bitsets[AFphi2_over]->operator [](startNode));
			assert(!ctl_under->bitsets[AFphi2_under]->operator [](startNode));
			learnAF(conflict, subf2, startNode);
			ctl_over->freeBitset(phi1_over);
			ctl_over->freeBitset(phi2_over);
			ctl_over->freeBitset(fAll_over);
			ctl_under->freeBitset(phi1_under);
			ctl_under->freeBitset(phi2_under);
			ctl_under->freeBitset(fAll_under);
			ctl_under->freeBitset(AFphi2_under);
			ctl_over->freeBitset(AFphi2_over);
			ctl_over->freeBitset(visited);
			return;
		}

		if(opt_verb>1)
			printf("learnAU: reconstructing path from %d to %d\n", startNode, from);

		if (from == startNode) {
			if(opt_verb>1)
				printf("learnAU: phi and psi fail to hold in startNode\n");
			learnClausePos(conflict, subf1, from); // learn that this initial node satisfy phi
			learnClausePos(conflict, subf2, from); // learn that this initial node satisfy psi
			ctl_over->freeBitset(phi1_over);
			ctl_over->freeBitset(phi2_over);
			ctl_over->freeBitset(fAll_over);
			ctl_under->freeBitset(phi1_under);
			ctl_under->freeBitset(phi2_under);
			ctl_under->freeBitset(fAll_under);
			ctl_over->freeBitset(visited);
			return;
		}
		learnClausePos(conflict, subf1, from); // learn that this reachable node satisfy phi
		learnClausePos(conflict, subf2, from); // learn that this reachable node satisfy psi

		while (from != startNode) {
			to = from;
			from = parent[from];
			eid = g_under->getEdge(from, to);

			learnClausePos(conflict, subf2, from); // learn that this reachable node satisfy psi

			// Learn that we can disable the transition
			Lit l = ~mkLit(getTransitionVar(eid), false);
			conflict.push(l);
			assert(value(l)==l_False);
		}
		ctl_over->freeBitset(phi1_over);
		ctl_over->freeBitset(phi2_over);
		ctl_over->freeBitset(fAll_over);
		ctl_under->freeBitset(phi1_under);
		ctl_under->freeBitset(phi2_under);
		ctl_under->freeBitset(fAll_under);
		ctl_over->freeBitset(visited);
	}



	void printLearntClause(vec<Lit> & c) {
      	for (int v = 0; v < c.size(); v++) {
      		if (sign(c[v])) {
      			printf("-%d ", var(c[v])+1);
      		} else {
          		printf("%d ", var(c[v])+1);
      		}
		}
		printf("     ");
		for (int v = 0; v < c.size(); v++) {
			if (sign(c[v])) {
				printf("-%d ", var(c[v])+1);
				if (vars[var(c[v])].type == EDGE)
					printf("-(Edge %d -> %d), ", g_over->getEdge(vars[var(c[v])].detector_node_edge).from, g_over->getEdge(vars[var(c[v])].detector_node_edge).to);
				if (vars[var(c[v])].type == NODEAP)
					printf("-(Node %d AP %d), ", vars[var(c[v])].detector_node_edge, vars[var(c[v])].ap);
			} else {
				printf("%d ", var(c[v])+1);
				if (vars[var(c[v])].type == EDGE)
					printf("(Edge %d -> %d), ", g_over->getEdge(vars[var(c[v])].detector_node_edge).from, g_over->getEdge(vars[var(c[v])].detector_node_edge).to);
				if (vars[var(c[v])].type == NODEAP)
					printf("(Node %d AP %d), ", vars[var(c[v])].detector_node_edge, vars[var(c[v])].ap);
			}
		}
		printf("\n");
	}

	void printFullClause() {
		printf("Full Clause:\n");

		vec<Lit> c;
		for (int v = 0; v < vars.size(); v++) {
			if(value(v)!=l_Undef){
				Lit l = ~mkLit(v,value(v)==l_False);
				assert(value(l)==l_False);
				c.push(l);
			}
		}

		printLearntClause(c);
	}

	bool solveTheory(vec<Lit> & conflict) {
		requiresPropagation = true;		//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict,true);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}
	;

	void drawFull(int from, int to) {
	}

	void printStats(int detailLevel) {
		printf("CTL Theory:\n");
		printf("\nUsed %d / %d bitset allocations for %d / %d bitsets\n", ctl_over->bitsetsmax, ctl_under->bitsetsmax, ctl_over->bitsetcounter, ctl_under->bitsetcounter);

		printf("Propagations: %ld (%f s, avg: %f s, %ld skipped)\n", stats_propagations, propagationtime,
				(propagationtime) / ((double) stats_propagations + 1), stats_propagations_skipped);
		printf("Decisions: %ld (%f s, avg: %f s)\n", stats_decisions, stats_decision_time,
				(stats_decision_time) / ((double) stats_decisions + 1));
		printf("Conflicts: %ld \n", stats_num_conflicts);

		printf("CTL propagation solve time: %f s\n", stats_solve_time);
		printf("CTL Conflicts: %d (%f lits/clause), time: %f s\n", stats_ctl_conflicts,stats_ctl_conflict_literals/((double) stats_symmetry_conflicts+1),stats_theory_prop_clause_learning_time);
		if(opt_ctl_symmetry>0){
			printf("Symmetry conflicts: %d (%f lits/clause), time: %f s\n", stats_symmetry_conflicts,stats_symmetry_conflict_literals/((double) stats_symmetry_conflicts+1), stats_theory_prop_symmetry_time);
		}

	}

	bool check_solved() {
		int infinitePaths = ctl_over->getReachableStates(initialNode);

		// hacked in something to print out the solution
		printf("\n--------------------\nSolution\n");
		g_under->draw(initialNode, -1, true);
		printf("\n--------------------\nSolution restricted to reachable states\n");
		g_under->drawRestricted(ctl_over->bitsets[infinitePaths], initialNode, -1, true);

		// Sanity check: Make sure that the solution is fully determined (no Undef) and the Kripke structure corresponds to the variables
		for (int edgeID = 0; edgeID < edge_labels.size(); edgeID++) {
			if (getTransition(edgeID).v < 0)
				continue;
			Var v = getTransition(edgeID).v;
			if(v==var_Undef)
				continue;
			lbool val = value(v);
			if (val == l_Undef) {
				return false;
			}

			if (val == l_True) {
				if (!g_under->transitionEnabled(edgeID)) {
					return false;
				}
				if (!g_over->transitionEnabled(edgeID)) {
					return false;
				}
			} else {
				if (g_under->transitionEnabled(edgeID)) {
					return false;
				}
				if (g_over->transitionEnabled(edgeID)) {
					return false;
				}
			}
		}
		for (int ap = 0; ap < g_over->apcount; ap++) {
			for (int state = 0; state < g_over->states(); state++) {
				Var v = getNodeAPVar(state, ap);
				lbool val = value(v);
				if (val == l_Undef) {
					return false;
				}
				if (val == l_True) {
					if (!g_under->isAPinStateLabel(state, ap)) {
						return false;
					}
					if (!g_over->isAPinStateLabel(state, ap)) {
						return false;
					}
				} else {
					if (g_under->isAPinStateLabel(state, ap)) {
						return false;
					}
					if (g_over->isAPinStateLabel(state, ap)) {
						return false;
					}
				}
			}
		}

		// Check solution against our own standalone CTL solver.
		//important change here - checking against the original formula, before any preprocessing or transformations!
		Bitset* bit_standalone_over = ctl_standalone_over->solveFormula(*original_f);
		Bitset* bit_standalone_under = ctl_standalone_under->solveFormula(*original_f);
		if(opt_verb>1)
			printf("check_solved making sure that solution agrees with standalone CTL solver...\n");

		if (value(ctl_lit)==l_True &&
				(!bit_standalone_over->operator [](initialNode) || !bit_standalone_under->operator [](initialNode))) {
			throw std::runtime_error("Solution found does not agree with standalone CTL model checker!");
		}
		if (value(ctl_lit)==l_False &&
				(bit_standalone_over->operator [](initialNode) || bit_standalone_under->operator [](initialNode))) {
			throw std::runtime_error("No solution despite the standalone CTL model checker stating that the current Kripke structure is a solution!");
		}

		printPctlOutput();
		printCTLSATOutput();

		if(opt_verb>1)
			printf("check_solved making sure that solution agrees with NuSMV solver...\n");

		// Even more paranoidly, check solution against nuSMV's CTL solver to make sure it is sound
		std::ofstream inputConvertedToNuSMVInput;
		//Sam: moved this to /tmp so that it will work if monosat is run from other directories.
		inputConvertedToNuSMVInput.open("/tmp/inputConvertedToNuSMVInput.txt", std::ios_base::out);
		inputConvertedToNuSMVInput << getNuSMVInput(ctl_over->bitsets[infinitePaths]);
		inputConvertedToNuSMVInput.close();
		std::system("NuSMV /tmp/inputConvertedToNuSMVInput.txt | grep Counterexample");

		// Solution passed all checks
		return true;
	}

	/*
	 * Generate NuSMV Input so that check_solved() can make sure that the solution found is indeed a solution
	 * The translation to the NuSMV input format is somewhat involved in that no states are allowed that don't have a successor.
	 * We thus compute the set of states satisfying AG EX True and only output those.
	 * */
	std::string getNuSMVInput(Bitset* infinitePaths) {
		/*
		 * Format for NuSMV:
MODULE main
VAR
  v0 : boolean;
  state : {0,1,2};
ASSIGN
  init(state) := 0;
  next(state) :=
    case
      state = 0 : {1,2};
      state = 1 : {1,2};
      state = 2 : {0,1};
    esac;
  v0 :=
    case
      state = 0 : TRUE;
      state = 1 : FALSE;
      state = 2 : FALSE;
    esac;

SPEC
  AG (v0 -> EX AF !v0)
		 *
		 * */

		std::string nuSMVInput = "MODULE main\nVAR\n  state: {"; // prints Output sentence on screen

		// Print set of states, but only those that belong on an infinite path
		for (int i = 0; i < g_under->states(); i++) {
			if (infinitePaths->operator [](i))
				nuSMVInput = nuSMVInput + std::to_string(i) + ", ";
		}
		nuSMVInput = nuSMVInput.substr(0, nuSMVInput.size()-2);
		nuSMVInput = nuSMVInput + "};\n";

		// print out all state variables
		for (int i = 0; i < g_under->statelabel[0]->size(); i++) {
			nuSMVInput = nuSMVInput + "  v" + std::to_string(i) + " : boolean;\n";
		}
		nuSMVInput += "ASSIGN\n  init(state) := 0;\n  next(state) :=\n    case\n";


		// Transitions
		DynamicGraph<int>::Edge e;
		int from, to, eid;
		for (int i = 0; i < g_under->states(); i++) {
			std::string nuSMVInputEdge = "       state = " +std::to_string(i) + " : {";
			for (int j = 0; j < g_under->nIncident(i); j++) { // iterate over neighbours of current front of queue
				e = g_under->incident(i, j);
				from = g_under->getEdge(e.id).from;
				to = g_under->getEdge(e.id).to;
				if (g_under->edgeEnabled(e.id)) {
					nuSMVInputEdge += std::to_string(to) + ", ";
				}
			}
			if (infinitePaths->operator [](i)) { // only do this for states that belong on an infinite paths
				nuSMVInputEdge = nuSMVInputEdge.substr(0, nuSMVInputEdge.size()-2);
				nuSMVInputEdge += "};\n";
				nuSMVInput += nuSMVInputEdge;
			}
		}
		nuSMVInput += "    esac;\n\n";

		// State literal assignments

		/*
  v0 :=
    case
      state = 0 : TRUE;
      state = 1 : FALSE;
      state = 2 : FALSE;
    esac;
    */

		for (int i = 0; i < g_under->statelabel[0]->size(); i++) {
			nuSMVInput += "  v"+std::to_string(i)+" :=\n    case\n";
			for (int j = 0; j < g_under->states(); j++) {
				if (infinitePaths->operator [](j)) { // only do this for states that belong on an infinite paths
					if (g_under->isAPinStateLabel(j, i))
						nuSMVInput += "       state = "+std::to_string(j)+" : TRUE;\n";
					else
						nuSMVInput += "       state = "+std::to_string(j)+" : FALSE;\n";
				}
			}
			nuSMVInput += "    esac;\n";
		}

		nuSMVInput += "\nSPEC\n  ";
		nuSMVInput += getFormulaNuSMVFormat(*original_f);
		nuSMVInput += "\n";
		return nuSMVInput;
	}

	// currently not called FIXME
	void printPctlOutput() {
		// This should probably not be here, but I can't figure out a better place:
		// Generate input file for PCTL-SMT
		std::string pctlInput = getFormulaPctlFormat(*original_f);

		// dirty hack to remove double negations, since PCTL-SMT doesn't like it in its input
		std::string notnot = "not not";
	    size_t start_pos = 0;
	    while((start_pos = pctlInput.find(notnot, start_pos)) != std::string::npos) {
	    	pctlInput.replace(start_pos, notnot.length(), "");
	    }
		//printf("pctlInput:\n\n");
		//std::cout << pctlInput;
		//printf("\n\n");

		std::ofstream inputConvertedToPctl;

		inputConvertedToPctl.open("regression-testing/inputConvertedToPctl.txt", std::ios_base::out);
		inputConvertedToPctl << pctlInput;
		inputConvertedToPctl.close();

/*		   char cwd[1024];
		   if (getcwd(cwd, sizeof(cwd)) != NULL)
		       fprintf(stdout, "Current working dir: %s\n", cwd);*/
		if (opt_verb > 0) {
			std::string pctlsyscall = "java -jar regression-testing/pctl.jar regression-testing/inputConvertedToPctl.txt "+std::to_string(nNodes())+" | regression-testing/yices-1.0.40/bin/yices | grep sat";
			printf("Run this command to verify result with PCTL, in case phi is in the intersection of PCTL and CTL:\n");
			std::cout << pctlsyscall;
			char pctlsyscallchar[1024];
			strncpy(pctlsyscallchar, pctlsyscall.c_str(), sizeof(pctlsyscallchar));
			pctlsyscallchar[sizeof(pctlsyscallchar) - 1] = 0;
			printf("\n\n");
			//std::system(pctlsyscallchar);
		}
	}

	void printCTLSATOutput() {
		std::string CTLSATInput = getFormulaCTLSATFormat(original_f);
		std::ofstream inputConvertedToCTLSAT;

		inputConvertedToCTLSAT.open("regression-testing/inputConvertedToCTLSAT.txt", std::ios_base::out);
		inputConvertedToCTLSAT << CTLSATInput;
		inputConvertedToCTLSAT.close();
		if (opt_verb > 0) {
			std::cout << CTLSATInput;
			printf("\n");
		}

	}


	void drawCurrentAssignment() {
		printf("digraph{\n");

		for (int i = 0; i < g_over->nodes(); i++) {
			std::string s = "node [label=\"" + std::to_string(i) + std::string(": {");
			for (int j = 0; j < g_over->apcount; j++) {
				if (g_over->statelabel[i]->operator [](j) == g_under->statelabel[i]->operator [](j))
					s += std::to_string(g_over->statelabel[i]->operator [](j));
				else
					s += "x";
			}
			s += "}\"] ";
			std::cout << s << "[shape=circle] " << i << ";\n";
		}

		if(initialNode>=0){
			printf("node [label=\"start\",shape=plaintext] start;\n");
			printf("start->%d\n",initialNode);
		}
		// determined enabled transitions
		for(int i = 0;i<g_over->transitions.size();i++){
			if (g_over->transitions[i] && g_under->transitions[i]){
				printf("%d->%d\n", g_over->getEdge(i).from,g_over->getEdge(i).to);
			}
		}
		// undetermined transitions
		for(int i = 0;i<g_over->transitions.size();i++){
			if (g_over->transitions[i] && !g_under->transitions[i]){
				printf("%d->%d [style=dashed]\n", g_over->getEdge(i).from,g_over->getEdge(i).to);
			}
		}
		printf("}\n");
	}


	bool dbg_solved() {

		return true;
	}

	void drawCurrent() {

	}
	int nEdges() {
		return edge_labels.size();
	}
	CRef newReasonMarker(int detectorID) {
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef - reasonMarker;
		marker_map.growTo(mnum + 1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}

	Lit newTransition(int from, int to, Var outerVar = var_Undef) {
		assert(outerVar!=var_Undef);
		int edgeID=-1;
		if(g_under->states()>from && g_under->states()>to && (edgeID=g_under->getEdge(from,to))>-1){
			if(g_over->transitionEnabled(edgeID)){
				//we already have this transition implemented
				Var ov =getTransition(edgeID).outerVar;
				assert(ov!=var_Undef);
				assert(getTransition(edgeID).from ==from);
				assert(getTransition(edgeID).to ==to);
				makeEqualInSolver(mkLit(outerVar),mkLit(ov));
				return lit_Undef;
			}
		}

		g_under->addTransition(from,to,edgeID,false);
		edgeID =g_over->addTransition(from,to,edgeID,true);
		Var v = newVar(outerVar,edgeID, -1, EDGE, true);

		edge_labels.growTo(edgeID+1);

		getTransition(edgeID).v = v;
		getTransition(edgeID).outerVar = outerVar;
		getTransition(edgeID).from = from;
		getTransition(edgeID).to = to;

		return mkLit(v, false);
	}

	void initNodeAPVarLookup(int nNodes, int nAp){


		nodeAPVarLookup.growTo(nNodes);
		for (int i = 0; i<nNodes; i++) {
			nodeAPVarLookup[i].growTo(nAp);
		}
	}

	Lit newNodeAP(int node, int ap, Var outerVar = var_Undef) {
		while (node >= g_under->statelabel.size()) {
			g_under->addEmptyState();
			g_over->addEmptyState();
		}
		g_under->disableAPinStateLabel(node,ap);
		g_over->enableAPinStateLabel(node,ap);

		Var nodeap = newVar(outerVar, node, ap, NODEAP, true);

		nodeAPVarLookup[node][ap] = nodeap;

		return mkLit(nodeap, false);
	}

	Lit newCTLVar(Var outerVar = var_Undef) {
		Var ctlVar = newVar(outerVar, 0, -1, DETECTOR, true);
		ctl_lit = mkLit(ctlVar, false);
		return ctl_lit;
	}

	void printSolution() {

		for (auto * d : detectors) {
			assert(d);
			d->printSolution();
		}
	}

	void setStrings(vec<vec<int>>* strings){
		assert(!this->strings);
		this->strings=strings;
	}
};

}
;

#endif
