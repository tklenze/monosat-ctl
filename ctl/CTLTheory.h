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


	CTLTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id){
		g_under = new DynamicKripke(id);
		g_over = new DynamicKripke(id);
		ctl_under = new CTLSolver(id, *g_under, *g_over);
		ctl_over = new CTLSolver(id, *g_over, *g_under);
		f = newCTLFormula();

		rnd_seed = opt_random_seed;
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
		//assert(vars[v].isEdge);
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
	
	Var newAuxVar(int forDetector = -1, VarType type = DETECTOR, bool connectToTheory = false) {
		Var s = S->newVar();
		return newVar(s, forDetector,type, connectToTheory);
	}
	Var newVar(Var solverVar,  int detector_node_edge, VarType type, bool connectToTheory = true) {
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();

		vars.push();
		vars[v].type = type;
		vars[v].detector_node_edge = detector_node_edge;
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
					assert(dbg_value(e.var)==l_Undef);
					int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
					assert(edgeID==e.ID);

					printf("Backtracking edge assignment %d\n", edgeID);


					if (e.assign) {
						g_under->disableTransition(edgeID);
					} else {
						g_over->enableTransition(edgeID);
					}
				} else if (e.type == NODEAP) {
					assert(dbg_value(e.var)==l_Undef);
					int nodeID = getNodeID(e.var);
					int apID = getAPID(e.var);
					assert(nodeID==e.ID);
					assert(apID==e.AP);
					printf("Backtracking nodeap assignment %d\n", nodeID);

					if (e.assign) {
						g_under->disableAPinStateLabel(nodeID, apID);
					} else {
						g_over->enableAPinStateLabel(nodeID, apID);
					}
				} else {
					//This is a detector literal
					//detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
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

		assert(dbg_graphsUpToDate());
		
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
		/*			g.markChanged();
		 antig.markChanged();
		 cutGraph.markChanged();*/
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

		printf("EnqueueTheory(Lit %d), var %d, new sign: %d\n", l.x, v, sign(l));

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
		
		/*
		for (int d = 0; d < detectors.size(); d++) {
			assert(conflict.size() == 0);
			bool r = detectors[d]->propagate(conflict);
			if (!r) {
				stats_num_conflicts++;
				toSolver(conflict);
				propagationtime += rtime(1) - startproptime;
				return false;
			}
		}
		*/

		printf("\n--------------------\npropagateTheory has been called\n");
		printf("Under:\n");
		g_under->draw();
		printf("\nOver:\n");
		g_over->draw();

		
		Bitset* bit_under = ctl_under->solve(*f);
		Bitset* bit_over = ctl_over->solve(*f);

		ctl_under->resetSwap();
		ctl_over->resetSwap();

		printf("\nSolution:\nUnder:\n");
		ctl_under->printStateSet(*bit_under);
		printf("Over:\n");
		ctl_over->printStateSet(*bit_over);

        if (value(ctl_lit)==l_True &&  !bit_over->operator [](initialNode)) {
        	printf("ctl_lit: %d, bit_over: %d", value(ctl_lit) == l_True, bit_over->operator [](initialNode));
        	printf("\npropagateTheory returns false, since formula is asserted true, but fails to hold in the overapproximation (and hence also fails to hold in the underapproximation) \n");
            return false; // It does not hold in the overapproximation
        }

        if (value(ctl_lit)==l_False &&  bit_under->operator [](initialNode)) {
        	printf("ctl_lit: %d, bit_over: %d", value(ctl_lit) == l_True, bit_over->operator [](initialNode));
            printf("\npropagateTheory returns false, since formula is asserted false, but it holds in the underapproximation (and hence also holds in the overapproximation)\n");
            return false; // It does not hold in the underapproximation
        }

		requiresPropagation = false;
		g_under->clearChanged();
		g_over->clearChanged();

		g_under->clearHistory();
		g_over->clearHistory();

		//detectors_to_check.clear();
		
		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;

		return true;
	}
	;

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
		return true; // Not implemented!


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
		for (int i = 0; i < detectors.size(); i++) {
			if (!detectors[i]->checkSatisfied()) {
				return false;
			}
		}
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
		Var v = newVar(outerVar,edgeID, EDGE, true);

		edge_labels.growTo(edgeID+1);

		getTransition(edgeID).v = v;
		getTransition(edgeID).outerVar = outerVar;
		getTransition(edgeID).from = from;
		getTransition(edgeID).to = to;

		return mkLit(v, false);
	}

	Lit newNodeAP(int node, int ap, Var outerVar = var_Undef) {
		g_under->disableAPinStateLabel(node,ap);
		g_over->enableAPinStateLabel(node,ap);

		Var nodeap = newVar(outerVar, node, NODEAP, true);

		return mkLit(nodeap, false);
	}

	// FIXME WTF am I doing here anyway?
	Lit newCTLVar(Var outerVar = var_Undef) {

		Var ctlVar = newVar(outerVar, 0, DETECTOR, true);

		ctl_lit = mkLit(ctlVar, false);

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
