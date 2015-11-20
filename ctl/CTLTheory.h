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
	int initialNode; // Start state
	Lit ctl_lit; // the Lit corresponding to the ctl formula

	//// Vector of CTL Formulas on the specific Kripke structure this solver deals with
	//vec<CTLFormula> formulas;
	
	vec<Assignment> trail;
	vec<int> trail_lim;


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


	// Compute AG EX True
	CTLFormula* fAGEXTrue;


	CTLTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id){
		g_under = new DynamicKripke(id);
		g_over = new DynamicKripke(id);
		ctl_under = new CTLSolver(id, *g_under, *g_over);
		ctl_over = new CTLSolver(id, *g_over, *g_under);
		ctl_standalone_over = new CTLSolverStandalone(id, *g_over);
		ctl_standalone_under = new CTLSolverStandalone(id, *g_under);
		f = newCTLFormula();

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
	// FIXME do we need this?
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
		assert(v < vars.size());
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

	void backtrackUntil(int level) { // FIXME important
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

	void preprocess() {
		for (int i = 0; i < detectors.size(); i++) {
			detectors[i]->preprocess();
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

	void enqueueTheory(Lit l) { // FIXME important
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
	bool propagateTheory(vec<Lit> & conflict) { // TODO this is the most important function
		static int itp = 0;
		if (++itp == 62279) {
			int a = 1;
		}
		stats_propagations++;

		if (!requiresPropagation) {
			stats_propagations_skipped++;
			assert(dbg_graphsUpToDate());
			return true;
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
		long bitsets = Bitset::remainingBitsets();
		Bitset* bit_under = ctl_under->solve(*f);
		Bitset* bit_over = ctl_over->solve(*f);
		long leaked = Bitset::remainingBitsets()-bitsets;
		//should leak exactly 2 bitsets ('bit_under' and 'bit_over') in the above calls
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call!");
		}
		ctl_under->resetSwap();
		ctl_over->resetSwap();

		// If we have a conflict, populate conflict set
		/*if ((value(ctl_lit)==l_True &&  !bit_over->operator [](initialNode)) ||
				(value(ctl_lit)==l_False &&  bit_under->operator [](initialNode))) {
			for (int i = 0; i < vars.size(); i++) {
				if (vars[i].occursNegative) {
					assert(value(mkLit(vars[i].solverVar, true))==l_False);
					conflict.push(mkLit(vars[i].solverVar, true));
				} else {
					assert(value(mkLit(vars[i].solverVar, false))==l_False);
					conflict.push(mkLit(vars[i].solverVar, false));
				}
			}
		}*/

		vec<Lit> symmetryConflict;//this allocates a new symmetryConflict vector at each propagateTheory call, which is needlessly expensive.
    	checkSymmetryConstraints(symmetryConflict, initialNode); // check if there is a conflict with the symmetry constraints under the current assignment, and if there is, build a clause to learn

        if (value(ctl_lit)==l_True &&  !bit_over->operator [](initialNode)) {

        	// THIS IS WHERE I DO CLAUSE LEARNING 1
        	learnClausePos(conflict, *f, initialNode);
        	/*  // old method, rather naive clause
        	for (int v = 0; v < vars.size(); v++) {
        		if(value(v)!=l_Undef){
        			Lit l = ~mkLit(v,value(v)==l_False);
        			assert(value(l)==l_False);
					conflict.push(l);
        		}
        	}
			*/

			if(opt_verb>1){
				printFullClause();
		  		printf("Learnt CTL Clause:\n");
				printLearntClause(conflict);
		  		printf("Learnt Symmetry Clause:\n");
		  		printLearntClause(symmetryConflict);
			}

			if (symmetryConflict.size() != 0 && symmetryConflict.size()<conflict.size()) {
				if(opt_verb>1)
					printf("Choosing to learn symmetry conflict rather than CTL conflict, since it is smaller (%d literals vs %d literals)", symmetryConflict.size(), conflict.size());
				symmetryConflict.copyTo(conflict);
			}

			toSolver(conflict);
			if(opt_verb>1){
				printf("ctl_lit: %d, bit_over: %d", value(ctl_lit) == l_True, bit_over->operator [](initialNode));
				printf("\npropagateTheory returns false, since formula is asserted true, but fails to hold in the overapproximation (and hence also fails to hold in the underapproximation) \n");
			}
			delete bit_under;
			delete bit_over;
        	return false; // It does not hold in the overapproximation
        } else if (value(ctl_lit)==l_False &&  bit_under->operator [](initialNode)) {

        	// THIS IS WHERE I DO CLAUSE LEARNING 2
        	// This is not really well implemented yet (i.e. no good clause learning, no symmetry reduction).

          	for (int v = 0; v < vars.size(); v++) {
				if(value(v)!=l_Undef){
					Lit l = ~mkLit(v,value(v)==l_False);
					assert(value(l)==l_False);
					conflict.push(l);
				}
			}
			toSolver(conflict);
			if(opt_verb>1){
				printf("ctl_lit: %d, bit_over: %d", value(ctl_lit) == l_True, bit_over->operator [](initialNode));
				printf("\npropagateTheory returns false, since formula is asserted false, but it holds in the underapproximation (and hence also holds in the overapproximation)\n");
			}
			delete bit_under;
			delete bit_over;
            return false; // It does not hold in the underapproximation
        } else if (symmetryConflict.size() != 0 && value(ctl_lit)!=l_Undef) { // CTL formula does not seem to be unsatisfiable, but there is a conflict with the symmetry constraints
        	if(opt_verb>1) {
        		printf("Learning symmetry conflict, since there is nothing else to learn\n");

        		printFullClause();
        		printf("Learnt Symmetry Clause:\n");
        		printLearntClause(symmetryConflict);

        	}
        	symmetryConflict.copyTo(conflict);
        	toSolver(conflict);

			delete bit_under;
			delete bit_over;
            return false; // Symmetry condition violated
        }

		requiresPropagation = false;
		g_under->clearChanged();
		g_over->clearChanged();

		g_under->clearHistory();
		g_over->clearHistory();

		//detectors_to_check.clear();
		
		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;
		delete bit_under;
		delete bit_over;
		return true;
	}
	;

	/*
	 * Do Symmetry Reduction. Check if any of the symmetry constraints are violated, and if so, build a clause describing this conflict.
	 * Calling this function might yield a conflict, but it does not have to if there is none.
	 */
	void checkSymmetryConstraints(vec<Lit> & conflict, int startNode) {
		vec<Lit> tmpConflict;
		for (int i = 0; i < g_over->states(); i++) {
			for (int j = i+1; j < g_over->states(); j++) {
				if (i != startNode && j != startNode) {
					if(opt_verb>1) {
						printf("SYMMETRY: checking %d > %d\n", i, j);
					}
					if (g_under->statelabel[i]->Equiv( *g_over->statelabel[j] )) {
						if(opt_verb>1) {
							printf("SYMMETRY: %d and %d have the same state label\n", i, j);
						}
					}

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
		assert(g_under->isAPinStateLabel(a, i) && !g_over->isAPinStateLabel(b, i)); // there must be something that differentiates them, otherwise a and b have the same statelabel!
	}


	// Clause learning for when a CTL formula is supposed to be true, but is false in the overapproximation
    // Starting from some initial state startNode (on the most toplevel call this will be the KripkeStructure's initial state,
	// but not necessarily so on recursive calls).
	void learnClausePos(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		if(opt_verb>1) {
			printf("Clause learning subformula... ");
			printFormula2(subf);
			printf(", for startNode %d\n", startNode);
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
			} else if (subf.operand1->op == AU) { // ~(φ_1 AU φ_2) ≡ ~φ_2 EW ( ~φ_1 ^ ~φ_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {EW, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else if (subf.operand1->op == AW) { // ~(φ_1 AW φ_2) ≡ ~φ_2 EU ( ~φ_1 ^ ~φ_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {EU, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else if (subf.operand1->op == EU) { // ~(φ_1 EU φ_2) ≡ ~φ_2 AW ( ~φ_1 ^ ~φ_2 )
				CTLFormula inner2 {NEG, subf.operand1->operand2, NULL, 0};
				CTLFormula outer1 {AND, &inner1, &inner2, 0};
				CTLFormula outer2 {AW, &inner2, &outer1, 0};
				learnClausePos(conflict, outer2, startNode);
			} else if (subf.operand1->op == EW) { // ~(φ_1 EW φ_2) ≡ ~φ_2 AU ( ~φ_1 ^ ~φ_2 )
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

			learnEG(conflict, *subf.operand1, startNode);
		}
		else if (subf.op == AG) {
			if(opt_verb>1)
				printf("Clause learning case AG...\n");

			learnAG(conflict, *subf.operand1, startNode);
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

		long bitsets = Bitset::remainingBitsets();
		Bitset* phi1_under = ctl_under->solve(subf1);
		Bitset* phi1_over = ctl_over->solve(subf1);
		long leaked = Bitset::remainingBitsets()-bitsets;
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call in learnAND.1!");
		}

		ctl_under->resetSwap();
		ctl_over->resetSwap();

		bitsets = Bitset::remainingBitsets();
		Bitset* phi2_under = ctl_under->solve(subf2);
		Bitset* phi2_over = ctl_over->solve(subf2);
		leaked = Bitset::remainingBitsets()-bitsets;
		//should leak exactly 2 bitsets ('bit_under' and 'bit_over') in the above calls
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call in learnAND.2!");
		}
		ctl_under->resetSwap();
		ctl_over->resetSwap();

		if(opt_verb>1) {
			printf("learnAND: phi1_over: "); ctl_standalone_over->printStateSet(*phi1_over);
			printf("learnAND: phi1_under: "); ctl_standalone_under->printStateSet(*phi1_under);
			printf("learnAND: phi2_over: "); ctl_standalone_over->printStateSet(*phi2_over);
			printf("learnAND: phi2_under: "); ctl_standalone_under->printStateSet(*phi2_under);
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
		if (!phi1_over->operator [](startNode) && !phi1_under->operator [](startNode) && !phi2_over->operator [](startNode) && !phi2_under->operator [](startNode)) {
			// Both subformulas conflict and are suitable to learn a clause from. Build both clauses and learn the smaller one.
			vec<Lit> conflict2;
			conflict.copyTo(conflict2);
			learnClausePos(conflict,  subf1, startNode);
			learnClausePos(conflict2, subf2, startNode);
			if (conflict.size() > conflict2.size()) {
				conflict2.copyTo(conflict);
			}
			//delete conflict2; // FIXME doesn't work...?
		} else if (!phi1_over->operator [](startNode) && !phi1_under->operator [](startNode)) {
			learnClausePos(conflict, subf1, startNode);
		} else {
			assert(!phi2_over->operator [](startNode) && !phi2_under->operator [](startNode)); // If this is violated, that means that both parts of the AND formula are satisfied -- then there should not be a conflict
			learnClausePos(conflict, subf2, startNode);
		}
		delete phi1_under;
		delete phi1_over;
		delete phi2_under;
		delete phi2_over;

	}

	// Both parts of the OR must be false, and we OR together their recursively learned sub-clauses
	void learnOR(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) {

		long bitsets = Bitset::remainingBitsets();
		Bitset* phi1_under = ctl_under->solve(subf1);
		Bitset* phi1_over = ctl_over->solve(subf1);
		long leaked = Bitset::remainingBitsets()-bitsets;
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call in learnOR.1!");
		}

		ctl_under->resetSwap();
		ctl_over->resetSwap();

		bitsets = Bitset::remainingBitsets();
		Bitset* phi2_under = ctl_under->solve(subf2);
		Bitset* phi2_over = ctl_over->solve(subf2);
		leaked = Bitset::remainingBitsets()-bitsets;
		//should leak exactly 2 bitsets ('bit_under' and 'bit_over') in the above calls
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call in learnOR.2!");
		}
		ctl_under->resetSwap();
		ctl_over->resetSwap();


		if(opt_verb>1) {
			printf("phi1_over: "); ctl_standalone_over->printStateSet(*phi1_over);
			printf("phi1_under: "); ctl_standalone_under->printStateSet(*phi1_under);
			printf("phi2_over: "); ctl_standalone_over->printStateSet(*phi2_over);
			printf("phi2_under: "); ctl_standalone_under->printStateSet(*phi2_under);
		}

		learnClausePos(conflict, subf1, startNode);
		learnClausePos(conflict, subf2, startNode);
		delete phi1_under;
		delete phi1_over;
		delete phi2_under;
		delete phi2_over;
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
				Lit l = ~mkLit(e.id, true);
				conflict.push(l);
				assert(value(l)==l_False);
			}
		}
	}

	// At least one neighbour enables phi, or at least one enabled edge to a neighbour is disabled
	// In order to make shorter clauses, we will only learn one clause: find an enabled edge to a state not satisfying
	//
	void learnAX(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		long bitsets = Bitset::remainingBitsets();
		Bitset* phi_under = ctl_under->solve(subf);
		Bitset* phi_over = ctl_over->solve(subf);
		long leaked = Bitset::remainingBitsets()-bitsets;
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call in learnAX.1!");
		}

		ctl_under->resetSwap();
		ctl_over->resetSwap();

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
			if (g_under->edgeEnabled(e.id) && ! phi_over->operator [](to)) {
				learnClausePos(conflict, subf, to);
				Lit l = ~mkLit(e.id, false);
				conflict.push(l);
				assert(value(l)==l_False);
				delete phi_under;
				delete phi_over;
				return;
			}
		}
		assert(false); // No clause was learned, which means that there is no enabled edge such that phi does not hold in the destination -- a contradiction to the fact, that AX phi does not hold!

	}


	// Successor to "phi-reachable" state satisfies phi or enable transition from any "phi-reachable" state
	//
	void learnEG(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		Bitset* phi = ctl_over->solve(subf); // Solve the inner subformula of the entire formula
		//ctl_over->printStateSet(*phi);
		//printFormula2(subf);

		DynamicGraph<int>::Edge e;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!phi->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnEG: initial state does not satisfy phi\n");
			learnClausePos(conflict, subf, startNode);
			delete phi;
			return;
		}

		std::queue <int> list; // this queue denotes all the states satisfying phi, whose neighbours have to be inspected
		list.push(startNode);
		Bitset *visited = new Bitset(g_over->states()); // This bitset denotes all the visited nodes
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
					if (phi->operator [](to) && !visited->operator [](to)) {
						list.push(to);

					} else if (!phi->operator [](to)) {
						learnClausePos(conflict, subf, to);
					} // else we already have visited the phi-state "to" and recursively all its neighbours
				} else {
					Lit l = ~mkLit(e.id, true);
					conflict.push(l);
					assert(value(l)==l_False);
				}
			}
			visited->set(from);
		}
		delete visited;
		delete phi;
	}


	// Find a path to a not-phi state, either that state must satisfy phi, or one of the edges must be disabled
	void learnAG(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		// Solve the inner subformula of the entire formula
		long bitsets = Bitset::remainingBitsets();
		Bitset* phi_under = ctl_under->solve(subf);
		Bitset* phi_over = ctl_over->solve(subf);
		long leaked = Bitset::remainingBitsets()-bitsets;
		if(leaked!=2){
			throw std::runtime_error("Leaked bitsets during solve call in learnAG.1!");
		}

		ctl_under->resetSwap();
		ctl_over->resetSwap();

		Bitset* phi_under_so = ctl_standalone_under->solve(subf);
		Bitset* phi_over_so = ctl_standalone_over->solve(subf);

		if(opt_verb>1) {
			printf("learnAG: phi_under "); ctl_standalone_over->printStateSet(*phi_under);
			printf("learnAG: phi_over "); ctl_standalone_over->printStateSet(*phi_over);
			printf("learnAG: phi_under_so "); ctl_standalone_over->printStateSet(*phi_under_so);
			printf("learnAG: phi_over_so "); ctl_standalone_over->printStateSet(*phi_over_so);
		}

		DynamicGraph<int>::Edge e;
		int eid = e.id;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!phi_under->operator [](startNode) && !phi_over->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnAG: initial state does not satisfy phi\n");
			learnClausePos(conflict, subf, startNode);
			delete phi_under_so;
			delete phi_over_so;
			delete phi_under;
			delete phi_over;
			return;
		}

		std::queue <int> s;
		s.push(startNode);
		std::map <int,int> parent; // from this we retreive the path of edges from start state to some state that doesn't satisfy phi
		Bitset *visited = new Bitset(g_over->states()); // This bitset denotes all the visited nodes
		bool done = false;
		while (s.size() > 0 && !done) { // stack of visited elements
			from = s.front();
			s.pop();
			visited->set(from);
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
						printf("learnAG: Adding map parent(%d) = %d. phi_under(to): %d, phi_over(to): %d, visited: %d\n", to, from, phi_under->operator [](to), phi_over->operator [](to), visited->operator [](to));
					if (!phi_under->operator [](to) && !phi_over->operator [](to)) {
						parent[to] = from;
						if(opt_verb>1)
							printf("learnAG: Found state no %d, which does not satisfy phi. edgeid: %d, from: %d, to: %d\n", to, e.id, from, to);
						learnClausePos(conflict, subf, to);
						done = true; // we have found a state that does not satisfy phi, exit loop and retreive path
					} else if (!visited->operator [](to)) {
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

			Lit l = ~mkLit(eid, false);
			conflict.push(l);
			assert(value(l)==l_False);

			to = from;
			from = parent[to];
			eid = g_under->getEdge(from, to);
		}
		delete phi_under_so;
		delete phi_over_so;
		delete phi_under;
		delete phi_over;
		delete visited;
	}

	// Find all reachable states, at least one of them should satisfy phi OR there should be a transition enabled
	// from a reachable state to an unreachable state
	void learnEF(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		Bitset* phi = ctl_over->solve(subf); // Solve the inner subformula of the entire formula

		DynamicGraph<int>::Edge e;
		int from, to;


		std::queue <int> list; // this queue denotes the current path of visited states, whose neighbours have to be inspected
		list.push(startNode);
		Bitset *visited = new Bitset(g_over->states()); // This bitset denotes all the visited nodes
		// We first employ BFS to find all reachable nodes. We add them to the clause
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			visited->set(from);
			if(opt_verb>1)
				printf("learnEF: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEF: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id) && !visited->operator [](to)) {
					list.push(to);
				}
			}
			learnClausePos(conflict, subf, from); // learn that this reachable node satisfy phi
		}
		// Now for part two of the clause:
		// iterate over edges of visited states to unvisited states; add literal to enable these edges
		for (int i = 0; i < visited->size(); i++) {
			if (visited->operator [](i)) {
				for (int j = 0; j < g_over->nIncident(i); j++) {
					e = g_over->incident(i, j);
					to = g_over->getEdge(e.id).to;

					if (!visited->operator [](to)) {
						if(opt_verb>1)
							printf("learnEF: Learning that the edge between %d and %d (edgeid: %d) could be enabled.\n", i, to, e.id);

						assert(!g_over->edgeEnabled(e.id)); // if edge was enabled, then j would be reachable, which means j would be visited
						Lit l = ~mkLit(e.id, true);
						conflict.push(l);
						assert(value(l)==l_False);
					}
				}
			}
		}
	}

	// Find lasso of states that satisfy not phi, at least one of them should satisfy phi OR one of the edges of the lasso should become disabled
	void learnAF(vec<Lit> & conflict, CTLFormula &subf, int startNode) {
		Bitset* phi_under = ctl_under->solve(subf); // Solve the inner subformula of the entire formula
		Bitset* phi_over = ctl_over->solve(subf); // Solve the inner subformula of the entire formula

		DynamicGraph<int>::Edge e;
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
				printf("learnAG: Considering state %d\n", from);

			for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAF: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !phi_under->operator [](to) && !phi_over->operator [](to)) {
					if (!visited->operator [](to)) { // explore new state in graph
					if(opt_verb>1)
						printf("learnAF: Adding map parent(%d) = %d. phi_under(to): %d, phi_over(to): %d, visited: %d. putting %d in queue\n", to, from, phi_under->operator [](to), phi_over->operator [](to), visited->operator [](to), to);

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

							Lit l = ~mkLit(eid, false);
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
		Lit l = ~mkLit(eid, false);
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
			Lit l = ~mkLit(eid, false);
			conflict.push(l);
			assert(value(l)==l_False);

			learnClausePos(conflict, subf, from); // learn that this reachable node satisfy phi

			to = from;
			from = parent[to];
			eid = g_under->getEdge(from, to);
		}
	}



	// Find all phi-reachable states, at least one of them should satisfy psi OR there should be a transition enabled
	// from a phi-reachable state (to anywhere) OR a non-phi-reachable state that is a successor to a phi-reachable state satisfies phi, or psi

	// For EU this should be pretty much the same, except that we only add transitions from phi-reachable to non-phi-reachable
	// FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME Not done yet
	void learnEW(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) { // phi EW psi
		Bitset* phi = ctl_over->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* psi = ctl_over->solve(subf2); // Solve the inner subformula of the entire formula


		DynamicGraph<int>::Edge e;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi or psi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!phi->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnEW: initial state does not satisfy phi or psi\n");
			assert(!psi->operator [](startNode)); // If psi was satisfied in initial state, then we would not have a conflict
			learnClausePos(conflict, subf1, startNode); // Either phi must be true...
			learnClausePos(conflict, subf2, startNode); // ... or psi must be true
			delete phi;
			delete psi;
			return;
		}


		std::queue <int> list; // this queue denotes the current path of visited states, whose neighbours have to be inspected
		list.push(startNode);
		Bitset *visited = new Bitset(g_over->states()); // This bitset denotes all the visited nodes
		// We first employ BFS to find all phi-reachable nodes. We add them to the clause
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			visited->set(from);
			if(opt_verb>1)
				printf("learnEW: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEW: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id) && !visited->operator [](to) && phi->operator [](to)) {
					list.push(to);
				} else if (g_over->edgeEnabled(e.id) && !visited->operator [](to) && !phi->operator [](to)) {
					learnClausePos(conflict, subf1, to); // learn that this successor to a phi-reachable state satisfy phi
					learnClausePos(conflict, subf2, to); // learn that this successor to a phi-reachable state satisfy psi
				}
			}
			assert(!psi->operator [](startNode)); // If psi was satisfied in this phi-reachable state, then we would not have a conflict
			learnClausePos(conflict, subf2, from); // learn that this phi-reachable node satisfy psi
		}
		// Now for part two of the clause:
		// iterate over edges from visited states; add literal to enable these edges
		for (int i = 0; i < visited->size(); i++) {
			if (visited->operator [](i)) {
				for (int j = 0; j < g_over->nIncident(i); j++) {
					e = g_over->incident(i, j);
					to = g_over->getEdge(e.id).to;

					if (!g_over->edgeEnabled(e.id)) {
						if(opt_verb>1)
							printf("learnEW: Learning that the edge between %d and %d (edgeid: %d) could be enabled.\n", i, to, e.id);
						Lit l = ~mkLit(e.id, true);
						conflict.push(l);
						assert(value(l)==l_False);
					}
				}
			}
		}
		delete phi;
		delete psi;
	}


	// Find all phi-reachable states, at least one of them should satisfy psi OR there should be a transition enabled
	// from a phi-reachable state (to anywhere) OR a non-phi-reachable state that is a successor to a phi-reachable state satisfies phi, or psi

	// For EU this should be pretty much the same, except that we only add transitions from phi-reachable to non-phi-reachable
	// FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME Not done yet
	void learnEU(vec<Lit> & conflict, CTLFormula &subf1, CTLFormula &subf2, int startNode) { // phi EW psi
		Bitset* phi = ctl_over->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* psi = ctl_over->solve(subf2); // Solve the inner subformula of the entire formula


		DynamicGraph<int>::Edge e;
		int from, to;

		// if the start node does not satisfy phi, then we just ask phi or psi to be satisfied here. This is treated separately, so
		// that in the queue later we can assume every item in the queue to satisfy phi.
		if (!phi->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnEU: initial state does not satisfy phi or psi\n");
			assert(!psi->operator [](startNode)); // If psi was satisfied in initial state, then we would not have a conflict
			learnClausePos(conflict, subf1, startNode); // Either phi must be true...
			learnClausePos(conflict, subf2, startNode); // ... or psi must be true
			delete phi;
			delete psi;
			return;
		}


		std::queue <int> list; // this queue denotes the current path of visited states, whose neighbours have to be inspected
		list.push(startNode);
		Bitset *visited = new Bitset(g_over->states()); // This bitset denotes all the visited nodes
		// We first employ BFS to find all phi-reachable nodes. We add them to the clause
		while (list.size() > 0) { // FIFO queue of visited elements
			from = list.front();
			list.pop();
			visited->set(from);
			if(opt_verb>1)
				printf("learnEU: Considering state %d\n", from);

			for (int i = 0; i < g_over->nIncident(from); i++) { // iterate over neighbours of current front of queue
				e = g_over->incident(from, i);
				to = g_over->getEdge(e.id).to;
				if(opt_verb>1)
					printf("learnEU: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_over->edgeEnabled(e.id) && !visited->operator [](to) && phi->operator [](to)) {
					list.push(to);
				} else if (g_over->edgeEnabled(e.id) && !visited->operator [](to) && !phi->operator [](to)) {
					learnClausePos(conflict, subf1, to); // learn that this successor to a phi-reachable state satisfy phi
					learnClausePos(conflict, subf2, to); // learn that this successor to a phi-reachable state satisfy psi
				}
			}
			assert(!psi->operator [](startNode)); // If psi was satisfied in this phi-reachable state, then we would not have a conflict
			learnClausePos(conflict, subf2, from); // learn that this phi-reachable node satisfy psi
		}
		// Now for part two of the clause:
		// iterate over edges from visited states; add literal to enable these edges
		for (int i = 0; i < visited->size(); i++) {
			if (visited->operator [](i)) {
				for (int j = 0; j < g_over->nIncident(i); j++) {
					e = g_over->incident(i, j);
					to = g_over->getEdge(e.id).to;

					if (!g_over->edgeEnabled(e.id) && !visited->operator [](to)) { // This is where the difference to EW is... we don't learn the edge if it merely leads to a phi-reachable state
						if(opt_verb>1)
							printf("learnEU: Learning that the edge between %d and %d (edgeid: %d) could be enabled.\n", i, to, e.id);
						Lit l = ~mkLit(e.id, true);
						conflict.push(l);
						assert(value(l)==l_False);
					}
				}
			}
		}
		delete phi;
		delete psi;
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
		Bitset* infinitePaths = ctl_under->solve(*fAGEXTrue);
		if(!infinitePaths->operator [](startNode)) {
			if(opt_verb>1)
				printf("learnAW: learning full clause, since AG EX True does not hold\n");
			learnFullClause(conflict);
			delete infinitePaths;
			return;
		}

		Bitset* phi1_under = ctl_under->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi1_over = ctl_over->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi2_under = ctl_under->solve(subf2); // Solve the inner subformula of the entire formula
		Bitset* phi2_over = ctl_over->solve(subf2); // Solve the inner subformula of the entire formula
		Bitset* fAll_under = ctl_under->solve(fAll); // Solve the inner subformula of the entire formula
		Bitset* fAll_over = ctl_over->solve(fAll); // Solve the inner subformula of the entire formula

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

							Lit l = ~mkLit(eid, false);
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
				Lit l = ~mkLit(eid, false);
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
		Bitset* phi1_under = ctl_under->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi1_over = ctl_over->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi2_under = ctl_under->solve(subf2); // Solve the inner subformula of the entire formula
		Bitset* phi2_over = ctl_over->solve(subf2); // Solve the inner subformula of the entire formula
		Bitset* fAll_under = ctl_under->solve(fAll); // Solve the inner subformula of the entire formula
		Bitset* fAll_over = ctl_over->solve(fAll); // Solve the inner subformula of the entire formula

		DynamicGraph<int>::Edge e;
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
				printf("learnAW/AU: Considering state %d\n", from);

			if (!phi1_over->operator [](from) && !phi2_over->operator [](from)) {
				done = true;
			}
			else for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAW/AU: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !fAll_under->operator [](to) && !fAll_over->operator [](to)) {
					if (!visited->operator [](to)) { // explore new state in graph
						if(opt_verb>1)
							printf("learnAW/AU: Adding map parent(%d) = %d. phi1_under(to): %d, phi1_over(to): %d, phi2_under(to): %d, phi2_over(to): %d, visited: %d. putting %d in queue\n", to, from, phi1_under->operator [](to), phi1_over->operator [](to), phi2_under->operator [](to), phi2_over->operator [](to), visited->operator [](to), to);

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
			delete phi1_under;
			delete phi1_over;
			delete phi2_under;
			delete phi2_over;
			delete fAll_under;
			delete fAll_over;
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
			Lit l = ~mkLit(eid, false);
			conflict.push(l);
			assert(value(l)==l_False);
		}
		delete phi1_under;
		delete phi1_over;
		delete phi2_under;
		delete phi2_over;
		delete fAll_under;
		delete fAll_over;
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
		Bitset* phi1_under = ctl_under->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi1_over = ctl_over->solve(subf1); // Solve the inner subformula of the entire formula
		Bitset* phi2_under = ctl_under->solve(subf2); // Solve the inner subformula of the entire formula
		Bitset* phi2_over = ctl_over->solve(subf2); // Solve the inner subformula of the entire formula
		Bitset* fAll_under = ctl_under->solve(fAll); // Solve the inner subformula of the entire formula
		Bitset* fAll_over = ctl_over->solve(fAll); // Solve the inner subformula of the entire formula

		DynamicGraph<int>::Edge e;
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
				printf("learnAU: Considering state %d\n", from);

			if (!phi1_over->operator [](from) && !phi2_over->operator [](from)) {
				done = true;
			}
			else for (int i = 0; i < g_under->nIncident(from) && !done; i++) { // iterate over neighbours of current front of queue
				e = g_under->incident(from, i);
				eid = e.id;
				to = g_under->getEdge(eid).to;
				if(opt_verb>1)
					printf("learnAU: Neighbour no %d, edgeid: %d, from: %d, to: %d\n", i, e.id, from, to);


				if (g_under->edgeEnabled(eid) && !fAll_under->operator [](to) && !fAll_over->operator [](to)) {
					if (!visited->operator [](to)) { // explore new state in graph
						if(opt_verb>1)
							printf("learnAU: Adding map parent(%d) = %d. phi1_under(to): %d, phi1_over(to): %d, phi2_under(to): %d, phi2_over(to): %d, visited: %d. putting %d in queue\n", to, from, phi1_under->operator [](to), phi1_over->operator [](to), phi2_under->operator [](to), phi2_over->operator [](to), visited->operator [](to), to);

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
			Bitset* AFphi2_under = ctl_under->solve(*fAFphi2); // Solve the inner subformula of the entire formula
			Bitset* AFphi2_over = ctl_over->solve(*fAFphi2); // Solve the inner subformula of the entire formula

			if(opt_verb>1) {
				printf("learnAU: These are the under and overapproximations for AF q: ");
				ctl_over->printStateSet(*AFphi2_under);
				ctl_over->printStateSet(*AFphi2_over);
				printf("\n");
			}

			assert(!AFphi2_over->operator [](startNode));
			assert(!AFphi2_under->operator [](startNode));
			learnAF(conflict, subf2, startNode);
			delete phi1_under;
			delete phi1_over;
			delete phi2_under;
			delete phi2_over;
			delete fAll_under;
			delete fAll_over;
			delete AFphi2_under;
			delete AFphi2_over;
			return;
		}

		if(opt_verb>1)
			printf("learnAU: reconstructing path from %d to %d\n", startNode, from);

		if (from == startNode) {
			if(opt_verb>1)
				printf("learnAU: phi and psi fail to hold in startNode\n");
			learnClausePos(conflict, subf1, from); // learn that this initial node satisfy phi
			learnClausePos(conflict, subf2, from); // learn that this initial node satisfy psi
			delete phi1_under;
			delete phi1_over;
			delete phi2_under;
			delete phi2_over;
			delete fAll_under;
			delete fAll_over;
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
			Lit l = ~mkLit(eid, false);
			conflict.push(l);
			assert(value(l)==l_False);
		}
		delete phi1_under;
		delete phi1_over;
		delete phi2_under;
		delete phi2_over;
		delete fAll_under;
		delete fAll_over;
	}



	void printLearntClause(vec<Lit> & c) {
		if(opt_verb>1){
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
	}

	void printFullClause() {
		if(opt_verb>1){
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
	}

	bool solveTheory(vec<Lit> & conflict) {
		requiresPropagation = true;		//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}
	;

	void drawFull(int from, int to) {
	}

	bool check_solved() {
		// TODO sanity check for NodeAPs. Thus far we only check for edges and whether the CTL formula is satisfied


		// hacked in something to print out the solution
		printf("\n--------------------\nSolution\n");
		printf("\nOver:\n");
		g_over->draw(0, -1, true);
		printf("Under:\n");
		g_under->draw(0, -1, true);

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
				/*	if(!g.hasEdge(e.from,e.to)){
						 return false;
						 }
						 if(!antig.hasEdge(e.from,e.to)){
						 return false;
						 }*/
				if (!g_under->transitionEnabled(edgeID)) {
					return false;
				}
				if (!g_over->transitionEnabled(edgeID)) {
					return false;
				}
			} else {
				/*if(g.hasEdge(e.from,e.to)){
						 return false;
						 }*/
				if (g_under->transitionEnabled(edgeID)) {
					return false;
				}
				if (g_over->transitionEnabled(edgeID)) {
					return false;
				}
				/*if(antig.hasEdge(e.from,e.to)){
						 return false;
						 }*/
			}
		}
		Bitset* bit_standalone_over = ctl_standalone_over->solve(*f);
		Bitset* bit_standalone_under = ctl_standalone_under->solve(*f);
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

		// Compute set of states satisfying AG EX True
		Bitset* infinitePaths = ctl_standalone_over->solve(*fAGEXTrue);

		if(opt_verb>1) {
			printf("Infinite paths AG EX True: ");
			ctl_standalone_over->printStateSet(*infinitePaths);
		}

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
		nuSMVInput += getFormulaNuSMVFormat(*f);
		nuSMVInput += "\n";

		std::ofstream inputConvertedToNuSMVInput;

		inputConvertedToNuSMVInput.open("regression-testing/inputConvertedToNuSMVInput.txt", std::ios_base::out);
		inputConvertedToNuSMVInput << nuSMVInput;
		inputConvertedToNuSMVInput.close();

		std::system("NuSMV regression-testing/inputConvertedToNuSMVInput.txt | grep Counterexample");


		return true;
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
		g_under->disableAPinStateLabel(node,ap);
		g_over->enableAPinStateLabel(node,ap);

		Var nodeap = newVar(outerVar, node, ap, NODEAP, true);

		nodeAPVarLookup[node][ap] = outerVar;

		return mkLit(nodeap, false);
	}

	// FIXME WTF am I doing here anyway?
	Lit newCTLVar(Var outerVar = var_Undef) {

		Var ctlVar = newVar(outerVar, 0, -1, DETECTOR, true);

		ctl_lit = mkLit(ctlVar, false);
		if(opt_verb>1)
			printf("Initialization: ctl_lit: true? %d, false? %d, undef? %d\n", value(ctl_lit) == l_True, value(ctl_lit) == l_False, value(ctl_lit) == l_Undef);

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

/*
	void addAcceptLit(int source ,int reach, int strID, Var outer_var){
		assert(g_under);
		accepts.growTo(source+1);
		if(!accepts[source]){
			accepts[source] = new FSMAcceptDetector(detectors.size(), this, g_under,g_over, source,*strings,drand(rnd_seed));
			detectors.push(accepts[source]);
		}
		accepts[source]->addAcceptLit(reach,strID,outer_var);
	}

	void addGenerateLit(int fsmID,int source, int strID, Var outer_var){
		assert(g_under[fsmID]);
		DynamicKripke & g_under = *g_under[fsmID];
		DynamicKripke & g_over = *g_over[fsmID];
		generates.growTo(source+1);
		if(!generates[source]){
			generates[source] = new FSMGeneratesDetector(detectors.size(), this, g_under,g_over, source,*strings,drand(rnd_seed));
			detectors.push(generates[source]);
		}
		generates[source]->addGeneratesLit(strID,outer_var);
	}
	void addTransduceLit(int fsmID, int source,int dest, int strID, int strID2, Var outer_var){
		assert(g_under[fsmID]);
		DynamicKripke & g_under = *g_under[fsmID];
		DynamicKripke & g_over = *g_over[fsmID];
		transduces.growTo(source+1);
		if(!transduces[source]){
			transduces[source] = new FSMTransducesDetector(detectors.size(), this, g_under,g_over, source,*strings,drand(rnd_seed));
			detectors.push(transduces[source]);
		}
		transduces[source]->addTransducesLit(dest,strID,strID2,outer_var);
	}
	*/
};

}
;

#endif
