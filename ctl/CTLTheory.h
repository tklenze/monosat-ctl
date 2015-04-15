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
#include "dgl/DynamicGraph.h"
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
	struct Transition{
		Var v=var_Undef;
		Var outerVar=var_Undef;
		int from=-1;
		int to=-1;
	};
	struct Assignment {
		bool isEdge :1;
		bool assign :1;
		int edgeID:30;
		Var var;
		Assignment(bool isEdge, bool assign, int edgeID, Var v):isEdge(isEdge),assign(assign),edgeID(edgeID),var(v){

		}
	};

	double rnd_seed;
	//vec<vec<int>> * strings=nullptr;
private:
	Solver * S;
	int local_q = 0;
public:
	int id;

	vec<lbool> assigns;


	//vec<vec<vec<Transition>>> edge_labels;
	DynamicKripke* g_unders;
	DynamicKripke* g_overs;

	CTLSolver ctl_under, ctl_over; //one for each over and under

	// Vector of CTL Formulas on the specific Kripke structure this solver deals with
	vec<CTLFormula> formulas;
	
	vec<Assignment> trail;
	vec<int> trail_lim;


public:
	//vec<FSMDetector*> detectors;
	//vec<FSMAcceptDetector*> reach_detectors;


	vec<int> marker_map;



	bool requiresPropagation = true;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {
		int isEdge :1;
		int occursPositive :1;
		int occursNegative :1;
		int detector_edge :29;	//the detector this variable belongs to, or its edge number, if it is an edge variable
		int input;
		int output;
		int fsmID;
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


		rnd_seed = opt_random_seed;
	}
	~CTLTheorySolver(){
		delete(g_unders);
		delete(g_overs);
		delete(CTLSolver);
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
		return vars[v].isEdge;
	}
	inline int getEdgeID(Var v) const {
		assert(isEdgeVar(v));
		return vars[v].detector_edge;
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
	
	Var newAuxVar(int forDetector = -1, bool connectToTheory = false) {
		Var s = S->newVar();
		return newVar(s, forDetector,-1, -1,-1,false, connectToTheory);
	}
	Var newVar(Var solverVar,  int detector_edge,int fsmID=-1, int label=-1,int output=-1, bool isEdge = false, bool connectToTheory = true) {
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();

		vars.push();
		vars[v].isEdge = isEdge;
		vars[v].fsmID=fsmID;
		vars[v].detector_edge = detector_edge;
		vars[v].solverVar = solverVar;
		vars[v].input=label;
		vars[v].output = output;
		assigns.push(l_Undef);
		if (connectToTheory) {
			S->setTheoryVar(solverVar, getTheoryIndex(), v);
			assert(toSolver(v) == solverVar);
		}
		if(isEdge){
			assert(label>-1);
			assert(detector_edge>-1);
		}
		if (!isEdge && detector_edge >= 0)
			detectors[detector_edge]->addVar(v);
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

	void backtrackUntil(int level) { // FIXME important
		static int it = 0;
		
		bool changed = false;
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if (e.isEdge) {
					assert(dbg_value(e.var)==l_Undef);
					int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
					assert(edgeID==e.edgeID);

					if (e.assign) {
						g_unders->disableTransition(edgeID);

					} else {
						g_overs->enableTransition(edgeID);

					}
				} else {
					//This is a CTL or a State AP assignment literal
					//FIXME do nothing for now
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
	
	void backtrackUntil(Lit p) { // FIXME important
		//need to remove and add edges in the two graphs accordingly.
		assert(value(p)==l_True);
		int i = trail.size() - 1;
		for (; i >= 0; i--) {
			Assignment e = trail[i];
			if (var(p) == e.var) {
				assert(sign(p) != e.assign);
				break;
			}
			if (e.isEdge) {
				int edgeID = getEdgeID(e.var); //e.var-min_edge_var;
				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
					g_unders->disableTransition(edgeID);
				} else {
					g_overs->enableTransition(edgeID);
				}
			} else {

				assigns[e.var] = l_Undef;
			}
		}
		
		trail.shrink(trail.size() - (i + 1));
		//if(i>0){
		requiresPropagation = true;
		/*			g.markChanged();
		 antig.markChanged();
		 cutGraph.markChanged();*/
		//}
		//while(trail_lim.size() && trail_lim.last()>=trail.size())
		//	trail_lim.pop();

		/*		for(int i = 0;i<reach_detectors.size();i++){
		 if(reach_detectors[i]->positive_reach_detector)
		 reach_detectors[i]->positive_reach_detector->update();
		 if(reach_detectors[i]->negative_reach_detector)
		 reach_detectors[i]->negative_reach_detector->update();
		 }*/
	}
	;

	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason) {
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
			assert(getTransition(fsmID,edgeID,input,output).v == var(l));

			trail.push( { true, !sign(l),edgeID, v });

			if (!sign(l)) {
				g_unders->enableTransition(edgeID);
			} else {
				g_overs->disableTransition(edgeID);
			}
			
		} else {
			
			trail.push( { false, !sign(l),-1, v });
		}
		
	}
	;
	bool propagateTheory(vec<Lit> & conflict) { // FIXME important
		static int itp = 0;
		if (++itp == 62279) {
			int a = 1;
		}
		stats_propagations++;

		if (!requiresPropagation) {
			stats_propagations_skipped++;
			//assert(dbg_graphsUpToDate());
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

		//assert(dbg_graphsUpToDate());
		
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
		

		
		requiresPropagation = false;
		g_unders->clearChanged();
		g_overs->clearChanged();

		g_unders->clearHistory();
		g_overs->clearHistory();

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

};

}
;

#endif
