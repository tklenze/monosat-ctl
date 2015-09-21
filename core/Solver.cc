/***************************************************************************************[Solver.cc]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless
 Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 Copyright (c) 2007-2010, Niklas Sorensson

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

#include <math.h>

#include "mtl/Sort.h"
#include "core/Solver.h"
#include "core/Config.h"
#include "graph/GraphTheory.h"
#include <unistd.h>
using namespace Monosat;

//=================================================================================================
// Options:
// Collected in Config.h

//=================================================================================================
// Constructor/Destructor:

Solver::Solver() :
		
		// Parameters (user settable):
		//
		verbosity(opt_verb), var_decay(opt_var_decay), clause_decay(opt_clause_decay), random_var_freq(
				opt_random_var_freq), random_seed(opt_random_seed), luby_restart(opt_luby_restart), ccmin_mode(
				opt_ccmin_mode), phase_saving(opt_phase_saving), rnd_pol(false), rnd_init_act(opt_rnd_init_act), garbage_frac(
				opt_garbage_frac), restart_first(opt_restart_first), restart_inc(opt_restart_inc)

		// Parameters (the rest):
		//
				, learntsize_factor((double) 1 / (double) 3), learntsize_inc(1.1)

		// Parameters (experimental):
		//
				, learntsize_adjust_start_confl(100), learntsize_adjust_inc(1.5)

		// Statistics: (formerly in 'SolverStats')
		//
				, solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0), stats_pure_lits(
				0), stats_pure_theory_lits(0), pure_literal_detections(0), stats_removed_clauses(0), dec_vars(0), clauses_literals(
				0), learnts_literals(0), max_literals(0), tot_literals(0), stats_pure_lit_time(0),  ok(
				true), cla_inc(1), var_inc(1), watches(WatcherDeleted(ca)), qhead(0), simpDB_assigns(-1), simpDB_props(
				0), order_heap(VarOrderLt(activity, priority)), progress_estimate(0), remove_satisfied(true)

		// Resource constraints:
		//
				, conflict_budget(-1), propagation_budget(-1) {
#ifdef DEBUG_SOLVER
	dbg_solver=NULL;
#endif
}

Solver::~Solver() {
	for (Theory * t : theories) {
		delete (t);
	}
}

//=================================================================================================
// Minor methods:

// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar) {
	
#ifdef DEBUG_SOLVER
	if( dbg_solver) {
		dbg_solver->newVar();
	}
#endif
	int v = nVars();
	watches.init(mkLit(v, false));
	watches.init(mkLit(v, true));
	assigns.push(l_Undef);
	vardata.push(mkVarData(CRef_Undef, 0));
	int p = 0;
	if (max_decision_var > 0 && v > max_decision_var)
		p = 1;
	if (v < min_decision_var)
		p = 1;
	priority.push(p);
	theory_vars.push();
	activity.push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
	seen.push(0);
	polarity.push(opt_init_rnd_phase ? irand(random_seed, 1) : sign);
	decision.push();
	trail.capacity(v + 1);
	if (max_decision_var > 0 && v > max_decision_var)
		dvar = false;
	if (v < min_decision_var)
		dvar = false;
	setDecisionVar(v, dvar);
	return v;
}

bool Solver::addClause_(vec<Lit>& ps) {
	
#ifdef DEBUG_SOLVER
	if( dbg_solver) {
		static vec<Lit> c;
		ps.copyTo(c);
		dbg_solver->addClause(c);
	}
#endif
	assert(decisionLevel() == 0);
	if (!ok)
		return false;
	resetInitialPropagation();    //Ensure that super solver call propagate on this solver at least once.
	// Check if clause is satisfied and remove false/duplicate literals:
	sort(ps);
	Lit p;
	int i, j;
	for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
		if (value(ps[i]) == l_True || ps[i] == ~p)
			return true;
		else if (value(ps[i]) != l_False && ps[i] != p)
			ps[j++] = p = ps[i];
	ps.shrink(i - j);
	
	if (ps.size() == 0)
		return ok = false;
	else if (ps.size() == 1) {
		uncheckedEnqueue(ps[0]);
		ok = (propagate(false) == CRef_Undef); //do NOT propagate theory solvers here, or else adding unit clauses can become very expensive in some circumstances (such as when constructing the initial CNF for example)
				
		return ok;
	} else {
		CRef cr = ca.alloc(ps, false);
		clauses.push(cr);
		attachClause(cr);
	}
	
	return true;
}

CRef Solver::attachClauseSafe(vec<Lit> & ps) {
	
	//sort(ps);
	Lit p;
	int i, j;
	for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
		if (((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
			return CRef_Undef;
		else if ((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p)
			ps[j++] = p = ps[i];
	ps.shrink(i - j);
	
	CRef confl_out = CRef_Undef;
	if (ps.size() == 0) {
		ok = false;
		cancelUntil(0);
		return CRef_Undef;
	} else if (ps.size() == 1) {
		cancelUntil(0);
		assert(var(ps[0]) < nVars());
		dbg_check_propagation(ps[0]);
		if (!enqueue(ps[0])) {
			ok = false;
			
		}
		return CRef_Undef;
	} else {
		//find the highest level in the conflict (should be the current decision level, but we won't require that)
		if (ps.size() > opt_temporary_theory_reasons) {
			if (tmp_clause == CRef_Undef) {
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else if (tmp_clause_sz < ps.size()) {
				ca[tmp_clause].mark(1); //is this needed?
				ca.free(tmp_clause);
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else {
				Clause & c = ca[tmp_clause];
				assert(tmp_clause_sz >= ps.size());
				assert(tmp_clause_sz >= c.size());
				c.grow(tmp_clause_sz - c.size());
				for (int i = 0; i < ps.size(); i++) {
					c[i] = ps[i];
				}
				c.shrink(c.size() - ps.size());
			}
			
			confl_out = tmp_clause;
		} else {
			CRef cr = ca.alloc(ps, !opt_permanent_theory_conflicts);
			ca[cr].setFromTheory(true);
			if (opt_permanent_theory_conflicts)
				clauses.push(cr);
			else {
				learnts.push(cr);
				if (--learntsize_adjust_cnt <= 0) {
					learntsize_adjust_confl *= learntsize_adjust_inc;
					learntsize_adjust_cnt = (int) learntsize_adjust_confl;
					max_learnts *= learntsize_inc;
					
					if (verbosity >= 1)
						printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %ld removed |\n", (int) conflicts,
								(int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
								(int) clauses_literals, (int) max_learnts, nLearnts(),
								(double) learnts_literals / nLearnts(), stats_removed_clauses);
				}
				
			}
			attachClause(cr);
			
			confl_out = cr;
		}
		return confl_out;
	}
}

void Solver::attachClause(CRef cr) {
	const Clause& c = ca[cr];
	assert(c.size() > 1);
#ifndef NDEBUG
	for (int p = 0; p <= 1; p++)
		if (value(c[p]) == l_False && level(var(c[p])) == 0) {
			//then c[0] must not have been propagated yet - it must be after qhead. Otherwise this clause will never be enforced
			if (decisionLevel() > 0) {
				assert(false); //because we have already propagated all the level 0 elements
			} else {
				bool found = false;
				for (int i = qhead; i < trail.size(); i++) {
					if (var(trail[i]) == var(c[p])) {
						found = true;
						break;
					}
				}
				assert(found);
			}
		}
#endif
	watches[~c[0]].push(Watcher(cr, c[1]));
	watches[~c[1]].push(Watcher(cr, c[0]));
	if (c.learnt())
		learnts_literals += c.size();
	else
		clauses_literals += c.size();
}

void Solver::detachClause(CRef cr, bool strict) {
	const Clause& c = ca[cr];
	assert(c.size() > 1);
	
	if (strict) {
		remove(watches[~c[0]], Watcher(cr, c[1]));
		remove(watches[~c[1]], Watcher(cr, c[0]));
	} else {
		// Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
		watches.smudge(~c[0]);
		watches.smudge(~c[1]);
	}
	
	if (c.learnt())
		learnts_literals -= c.size();
	else
		clauses_literals -= c.size();
}

void Solver::removeClause(CRef cr) {
	Clause& c = ca[cr];
	detachClause(cr);
	// Don't leave pointers to free'd memory!
	if (locked(c))
		vardata[var(c[0])].reason = CRef_Undef;
	c.mark(1);
	ca.free(cr);
}

bool Solver::satisfied(const Clause& c) const {
	for (int i = 0; i < c.size(); i++)
		if (value(c[i]) == l_True)
			return true;
	return false;
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
	if (decisionLevel() > level) {
		for (int c = trail.size() - 1; c >= trail_lim[level]; c--) {
			Var x = var(trail[c]);
			assigns[x] = l_Undef;
			if (phase_saving > 1 || ((phase_saving == 1) && c > trail_lim.last()))
				polarity[x] = sign(trail[c]);
			insertVarOrder(x);
		}
		qhead = trail_lim[level];
		if (local_qhead > qhead) {
			local_qhead = qhead;
		}
		if (S && super_qhead > S->qhead) {
			super_qhead = S->qhead;
		}
		
		trail.shrink(trail.size() - trail_lim[level]);
		trail_lim.shrink(trail_lim.size() - level);
		
		if (decisionLevel() < track_min_level) {
			track_min_level = decisionLevel();
		}
		for (int i : theory_queue) {
			in_theory_queue[i] = false;
		}
		theory_queue.clear();
		for (int i = 0; i < theories.size(); i++) {
			theories[i]->backtrackUntil(level);
		}
	}
}

void Solver::backtrackUntil(int level) {
	if (S->trail.size() < super_qhead)
		super_qhead = S->trail.size();
}

//=================================================================================================
// Major methods:

Lit Solver::pickBranchLit() {
	Var next = var_Undef;
	
	// Random decision:
	if (drand(random_seed) < random_var_freq && !order_heap.empty()) {
		next = order_heap[irand(random_seed, order_heap.size())];
		if (value(next) == l_Undef && decision[next])
			rnd_decisions++;
	}
	
	// Activity based decision:
	while (next == var_Undef || value(next) != l_Undef || !decision[next])
		if (order_heap.empty()) {
			next = var_Undef;
			break;
		} else
			next = order_heap.removeMin();
	
	return next == var_Undef ? lit_Undef : mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]);
}

/*_________________________________________________________________________________________________
 |
 |  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
 |  
 |  Description:
 |    Analyze conflict and produce a reason clause.
 |  
 |    Pre-conditions:
 |      * 'out_learnt' is assumed to be cleared.
 |      * Current decision level must be greater than root level.
 |  
 |    Post-conditions:
 |      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
 |      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
 |        rest of literals. There may be others from the same level though.
 |  
 |________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel) {
	int pathC = 0;
	Lit p = lit_Undef;
	assert(confl != CRef_Undef);
	// Generate conflict clause:
	//
	out_learnt.push();      // (leave room for the asserting literal)
	int index = trail.size() - 1;
	
	do {
		if (confl != CRef_Undef) {
			assert(!isTheoryCause(confl));
			Clause& c = ca[confl];
			
			if (c.learnt())
				claBumpActivity(c);
			
			for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
				Lit q = c[j];
				
				if (!seen[var(q)] && level(var(q)) > 0) {
					varBumpActivity(var(q));
					seen[var(q)] = 1;
					if (level(var(q)) >= decisionLevel())
						pathC++;
					else
						out_learnt.push(q);
				}
			}
		} else {
			out_learnt.push(~p);
		}
		// Select next clause to look at:
		while (!seen[var(trail[index--])])
			;
		assert(index >= -1);
		p = trail[index + 1];
		confl = reason(var(p));
		if (isTheoryCause(confl)) {
			//lazily construct the reason for this theory propagation now that we need it
			confl = constructReason(p);
		}
		seen[var(p)] = 0;
		pathC--;
		
	} while (pathC > 0);
	out_learnt[0] = ~p;
	dbg_check(out_learnt);
#ifndef NDEBUG
	for (Lit p : out_learnt)
		assert(value(p)==l_False);
#endif
	// Simplify conflict clause:
	//
	int i, j;
	out_learnt.copyTo(analyze_toclear);
	if (ccmin_mode == 2) {
		uint32_t abstract_level = 0;
		for (i = 1; i < out_learnt.size(); i++)
			abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)
					
		for (i = j = 1; i < out_learnt.size(); i++)
			if (!ca.isClause(reason(var(out_learnt[i]))) || !litRedundant(out_learnt[i], abstract_level))
				out_learnt[j++] = out_learnt[i];
		
	} else if (ccmin_mode == 1) {
		for (i = j = 1; i < out_learnt.size(); i++) {
			Var x = var(out_learnt[i]);
			
			if (reason(x) == CRef_Undef)
				out_learnt[j++] = out_learnt[i];
			else {
				Clause& c = ca[reason(var(out_learnt[i]))];
				for (int k = 1; k < c.size(); k++)
					if (!seen[var(c[k])] && level(var(c[k])) > 0) {
						out_learnt[j++] = out_learnt[i];
						break;
					}
			}
		}
	} else
		i = j = out_learnt.size();
	
	max_literals += out_learnt.size();
	out_learnt.shrink(i - j);
	tot_literals += out_learnt.size();
	
	// Find correct backtrack level:
	//
	if (out_learnt.size() == 1)
		out_btlevel = 0;
	else {
		int max_i = 1;
		// Find the first literal assigned at the next-highest level:
		for (int i = 2; i < out_learnt.size(); i++)
			if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
				max_i = i;
		// Swap-in this literal at index 1:
		Lit p = out_learnt[max_i];
		out_learnt[max_i] = out_learnt[1];
		out_learnt[1] = p;
		out_btlevel = level(var(p));
	}
#ifndef NDEBUG
	for (Lit p : out_learnt)
		assert(value(p)==l_False);
#endif
	dbg_check(out_learnt);
	for (int j = 0; j < analyze_toclear.size(); j++)
		seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels) {
	analyze_stack.clear();
	analyze_stack.push(p);
	int top = analyze_toclear.size();
	while (analyze_stack.size() > 0) {
		assert(reason(var(analyze_stack.last())) != CRef_Undef);
		if (!ca.isClause(reason(var(analyze_stack.last())))) {
			//Don't pursue this if it is a subsolver reason
			for (int j = top; j < analyze_toclear.size(); j++)
				seen[var(analyze_toclear[j])] = 0;
			analyze_toclear.shrink(analyze_toclear.size() - top);
			analyze_stack.clear();
			return false;
		}
		
		Clause& c = ca[reason(var(analyze_stack.last()))];
		analyze_stack.pop();
		
		for (int i = 1; i < c.size(); i++) {
			Lit p = c[i];
			if (!seen[var(p)] && level(var(p)) > 0) {
				if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0) {
					seen[var(p)] = 1;
					analyze_stack.push(p);
					analyze_toclear.push(p);
				} else {
					for (int j = top; j < analyze_toclear.size(); j++)
						seen[var(analyze_toclear[j])] = 0;
					analyze_toclear.shrink(analyze_toclear.size() - top);
					return false;
				}
			}
		}
	}
	
	return true;
}

/*_________________________________________________________________________________________________
 |
 |  analyzeFinal : (p : Lit)  ->  [void]
 |  
 |  Description:
 |    Specialized analysis procedure to express the final conflict in terms of assumptions.
 |    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
 |    stores the result in 'out_conflict'.
 |________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict) {
	out_conflict.clear();
	out_conflict.push(p);
	
	if (decisionLevel() == 0)
		return;
	
	seen[var(p)] = 1;
	
	for (int i = trail.size() - 1; i >= trail_lim[0]; i--) {
		Var x = var(trail[i]);
		if (seen[x]) {
			if (reason(x) == CRef_Undef) {
				assert(level(x) > 0);
				out_conflict.push(~trail[i]);
			} else {
				if (isTheoryCause(reason(x))) {
					constructReason(trail[i]);
				}
				Clause& c = ca[reason(x)];
				assert(var(c[0]) == x);
				for (int j = 1; j < c.size(); j++)
					if (level(var(c[j])) > 0)
						seen[var(c[j])] = 1;
			}
			seen[x] = 0;
		}
	}
	
	seen[var(p)] = 0;
}

void Solver::uncheckedEnqueue(Lit p, CRef from) {
	assert(value(p) == l_Undef);
	assigns[var(p)] = lbool(!sign(p));
	vardata[var(p)] = mkVarData(from, decisionLevel());
	trail.push_(p);
	//printf("%d\n",dimacs(p));
	if (hasTheory(p)) {
		int theoryID = getTheoryID(p);
		theories[theoryID]->enqueueTheory(getTheoryLit(p));
		if (!in_theory_queue[theoryID]) {
			in_theory_queue[theoryID] = true;
			theory_queue.push(theoryID);
			assert(theory_queue.size() <= theories.size());
		}
	}
}

void Solver::analyzeFinal(CRef confl, Lit skip_lit, vec<Lit>& out_conflict) {
	out_conflict.clear();
	if (decisionLevel() == 0)
		return;
	
	Clause & c = ca[confl];
	
	for (int i = 0; i < c.size(); i++) {
		Var x = var(c[i]);
		if (x == var(skip_lit))
			continue;
		
		assert(x >= 0);
		assert(x < nVars());
		if (level(x) > 0)
			seen[x] = 1;
		assert(value(x)!=l_Undef);
	}
	
	int start = trail.size() - 1;
	int i;
	for (i = start; i >= trail_lim[0]; i--) {
		Var x = var(trail[i]);
		assert(value(x)!=l_Undef);
		assert(x >= 0);
		assert(x < nVars());
		
		if (seen[x]) {
			assert(x != var(skip_lit));
			CRef r = reason(x);
			if (r == CRef_Undef) {
				Var v = var(trail[i]);
				int lev = level(v);
				out_conflict.push(~trail[i]);
			} else {
				if (isTheoryCause(r)) {
					r = constructReason(trail[i]);
				}
				Clause& c = ca[r];
				for (int j = 0; j < c.size(); j++) {
					
					if (level(var(c[j])) > 0) {
						seen[var(c[j])] = 1;
					}
				}
			}
			seen[x] = 0;
		}
	}
}

void Solver::buildReason(Lit p, vec<Lit> & reason) {
	Lit local_l = fromSuper(p);    // mkLit(var(p)-super_offset, sign(p));
	analyzeFinal(local_l, reason);
	toSuper(reason, reason);
	interpolant.push();
	reason.copyTo(interpolant.last());    //Add this clause into the interpolant vector
}

void Solver::enqueueTheory(Lit l) {
	
}

//Propagate assignments from the super solver's interface variables to this solver (and, if this solver makes further assignments to the interface, pass those back to the super solver)
bool Solver::propagateTheory(vec<Lit> & conflict_out) {
	//backtrack as needed
	if (S->trail.size() && super_qhead > 0) {
		Lit last_super = S->trail[super_qhead - 1];
		cancelUntil(S->level(var(last_super)));
	} else {
		cancelUntil(0);
	}
	assert(decisionLevel() <= S->decisionLevel());
	
	CRef confl = CRef_Undef;
	int curLev = decisionLevel();
	while (confl == CRef_Undef && super_qhead < S->qhead) {
		Lit out_l = S->trail[super_qhead++];
		
		if (var(out_l) < min_super || var(out_l) > max_super)
			continue;    //this lit is not on the interface
			
		int lev = S->level(var(out_l));
		assert(decisionLevel() <= lev);
		while (decisionLevel() < lev) {
			
			//We are going to start a new decision level; make sure the current one is fully propagated first
			initial_level = decisionLevel();
			track_min_level = initial_level;
			
			confl = propagate();
			if (confl != CRef_Undef || track_min_level < initial_level) {
				cancelUntil(track_min_level);
				S->cancelUntil(track_min_level);
				goto conflict;
			}
			
			newDecisionLevel();
		}
		
		Lit local_l = fromSuper(out_l);
		
		if (!enqueue(local_l, CRef_Undef)) {
			//this is a conflict
			conflict_out.clear();
			analyzeFinal(~local_l, conflict_out);
			toSuper(conflict_out, conflict_out);
			interpolant.push();
			conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
			return false;
		}
		
	}
	initial_level = decisionLevel();
	track_min_level = initial_level;
	confl = propagate();
	if (confl != CRef_Undef || track_min_level < initial_level) {
		cancelUntil(track_min_level);
		S->cancelUntil(track_min_level);
	}
	conflict:

	if (confl != CRef_Undef) {
		//then we have a conflict which we need to instantiate in S
		conflict_out.clear();
		analyzeFinal(confl, lit_Undef, conflict_out);
		toSuper(conflict_out, conflict_out);
		interpolant.push();
		conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
		return false;
	} else {
		//find lits to prop
		while (local_qhead < qhead) {
			Lit local_l = trail[local_qhead++];
			if (var(local_l) < min_local || var(local_l) > max_local)
				continue;    //this lit is not on the interface
			Lit out_l = toSuper(local_l);
			lbool v = S->value(out_l);
			
			if (S->value(out_l) != l_True || S->level(var(out_l)) > level(var(local_l))) {
				
				if (S->decisionLevel() > level(var(local_l))) {
					cancelUntil(level(var(local_l)));
					S->cancelUntil(level(var(local_l)));
					
				}
				if (!S->enqueue(out_l, cause_marker)) {    //can probably increment super_qhead here...
					//this is a conflict
					conflict_out.clear();
					analyzeFinal(local_l, conflict_out);
					toSuper(conflict_out, conflict_out);
					interpolant.push();
					conflict_out.copyTo(interpolant.last());			//Add this clause into the interpolant vector
					return false;
				}
			}
		}
		
		while (decisionLevel() < S->decisionLevel())
			newDecisionLevel();
		
		assert(decisionLevel() == S->decisionLevel());
		return true;
	}
}
/*_________________________________________________________________________________________________
 |
 |  propagate : [void]  ->  [Clause*]
 |  
 |  Description:
 |    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
 |    otherwise CRef_Undef.
 |  
 |    Post-conditions:
 |      * the propagation queue is empty, even if there was a conflict.
 |________________________________________________________________________________________________@*/
CRef Solver::propagate(bool propagate_theories) {
	if (qhead == trail.size() && !initialPropagate) {
		assert(!theory_queue.size());
		return CRef_Undef;
	}
	CRef confl = CRef_Undef;
	int num_props = 0;
	watches.cleanAll();
	if (decisionLevel() == 0 && !propagate_theories) {
		initialPropagate = true;//we will need to propagate this assignment to the theories at some point in the future.
	}
	do {
		
		while (qhead < trail.size()) {
			if (opt_early_theory_prop) {
				//propagate theories;
				while (propagate_theories && theory_queue.size() && confl == CRef_Undef) {
					theory_conflict.clear();
					int theoryID = theory_queue.last();
					theory_queue.pop();
					in_theory_queue[theoryID] = false;
					if (!theories[theoryID]->propagateTheory(theory_conflict)) {
						if (!addConflictClause(theory_conflict, confl)) {
							in_theory_queue[theoryID] = false;
							qhead = trail.size();
							return confl;
						}
					}
				}
			}
			
			Lit p = trail[qhead++];     // 'p' is enqueued fact to propagate.
			vec<Watcher>& ws = watches[p];
			Watcher *i, *j, *end;
			num_props++;
			
			for (i = j = (Watcher*) ws, end = i + ws.size(); i != end;) {
				// Try to avoid inspecting the clause:
				Lit blocker = i->blocker;
				if (value(blocker) == l_True) {
					*j++ = *i++;
					continue;
				}
				
				// Make sure the false literal is data[1]:
				CRef cr = i->cref;
				Clause& c = ca[cr];
				Lit false_lit = ~p;
				if (c[0] == false_lit)
					c[0] = c[1], c[1] = false_lit;
				assert(c[1] == false_lit);
				i++;
				
				// If 0th watch is true, then clause is already satisfied.
				Lit first = c[0];
				Watcher w = Watcher(cr, first);
				if (first != blocker && value(first) == l_True) {
					*j++ = w;
					continue;
				}
				
				// Look for new watch:
				for (int k = 2; k < c.size(); k++)
					if (value(c[k]) != l_False) {
						c[1] = c[k];
						c[k] = false_lit;
						watches[~c[1]].push(w);
						goto NextClause;
					}
				
				// Did not find watch -- clause is unit under assignment:
				*j++ = w;
				if (value(first) == l_False) {
					confl = cr;
					qhead = trail.size();
					// Copy the remaining watches:
					while (i < end)
						*j++ = *i++;
				} else
					uncheckedEnqueue(first, cr);
				
				NextClause: ;
			}
			ws.shrink(i - j);
		}
		
		if (initialPropagate && propagate_theories) {
			assert(decisionLevel() == 0);
			//propagate any as yet unpropagated literals to each theory
			for (int i = 0; i < qhead; i++) {
				Lit p = trail[i];
				if (hasTheory(p)) {
					int theoryID = getTheoryID(p);
					theories[theoryID]->enqueueTheory(getTheoryLit(p));
					
				}
			}
			//Have to check _all_ the theories, even if we haven't eneueued of their literals, in case they have constants to propagate up.
			for (int theoryID = 0; theoryID < theories.size(); theoryID++) {
				if (!in_theory_queue[theoryID]) {
					in_theory_queue[theoryID] = true;
					theory_queue.push(theoryID);
				}
			}
			initialPropagate = false;
		}
		static int iter = 0;
		if (++iter == 1452024) {
			int a = 1;
		}
		//propagate theories;
		while (propagate_theories && theory_queue.size() && (opt_early_theory_prop || qhead == trail.size())
				&& confl == CRef_Undef) {
			theory_conflict.clear();
			
			int theoryID = theory_queue.last();
			theory_queue.pop();
			in_theory_queue[theoryID] = false;
			if (!theories[theoryID]->propagateTheory(theory_conflict)) {
				if (!addConflictClause(theory_conflict, confl)) {
					in_theory_queue[theoryID] = false;
					qhead = trail.size();
					return confl;
				}else{
					throw std::runtime_error("Learnt clause is satisfiable!");
				}
			}
		}
		
		//solve theories if this solver is completely assigned
		/*	for(int i = 0;i<theories.size() && qhead == trail.size() && confl==CRef_Undef && nAssigns()==nVars();i++){
		 if(opt_subsearch==3 && track_min_level<initial_level)
		 continue;//Disable attempting to solve sub-solvers if we've backtracked past the super solver

		 if(!theories[i]->solve(theory_conflict)){
		 if(!addConflictClause(theory_conflict,confl))
		 return confl;

		 }
		 }*/

	} while (qhead < trail.size());
	propagations += num_props;
	simpDB_props -= num_props;
	
	return confl;
}

/*_________________________________________________________________________________________________
 |
 |  reduceDB : ()  ->  [void]
 |  
 |  Description:
 |    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
 |    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
 |________________________________________________________________________________________________@*/
struct reduceDB_lt {
	ClauseAllocator& ca;
	reduceDB_lt(ClauseAllocator& ca_) :
			ca(ca_) {
	}
	bool operator ()(CRef x, CRef y) {
		return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity());
		/*if(ca[x].size()==2)
		 return false;
		 if(ca[y].size()==2)
		 return true;
		 if(ca[x].fromTheory() && !ca[y].fromTheory())
		 return true;
		 if(ca[y].fromTheory() && ! ca[x].fromTheory())
		 return false;
		 return ca[x].activity() < ca[y].activity();*/
	}
	
};
void Solver::reduceDB() {
	int i, j;
	double extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity
			
	sort(learnts, reduceDB_lt(ca));
	// Don't delete binary or locked clauses. From the rest, delete clauses from the first half
	// and clauses with activity smaller than 'extra_lim':
	for (i = j = 0; i < learnts.size(); i++) {
		Clause& c = ca[learnts[i]];
		if (c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim)) {
			stats_removed_clauses++;
			removeClause(learnts[i]);
		} else
			learnts[j++] = learnts[i];
	}
	learnts.shrink(i - j);
	checkGarbage();
}

void Solver::removeSatisfied(vec<CRef>& cs) {
	int i, j;
	for (i = j = 0; i < cs.size(); i++) {
		Clause& c = ca[cs[i]];
		if (satisfied(c))
			removeClause(cs[i]);
		else
			cs[j++] = cs[i];
	}
	cs.shrink(i - j);
}

void Solver::rebuildOrderHeap() {
	vec<Var> vs;
	for (Var v = 0; v < nVars(); v++)
		if (decision[v] && value(v) == l_Undef)
			vs.push(v);
	order_heap.build(vs);
}

/*_________________________________________________________________________________________________
 |
 |  simplify : [void]  ->  [bool]
 |  
 |  Description:
 |    Simplify the clause database according to the current top-level assigment. Currently, the only
 |    thing done here is the removal of satisfied clauses, but more things can be put here.
 |________________________________________________________________________________________________@*/
bool Solver::simplify() {
	assert(decisionLevel() == 0);
	
	if (!ok || propagate() != CRef_Undef || !ok) //yes, the second ok check is now required, because propagation of a theory can make the solver unsat at this point...
		return ok = false;
	
	if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
		return true;
	
	// Remove satisfied clauses:
	removeSatisfied(learnts);
	if (remove_satisfied)        // Can be turned off.
		removeSatisfied(clauses);
	checkGarbage();
	rebuildOrderHeap();
	
	simpDB_assigns = nAssigns();
	simpDB_props = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)
			
	//Detect any theory literals that are occur only in positive or negative (or neither) polarity
	if (opt_detect_pure_lits) {
		//if(pure_literal_detections==0){
		
		pure_literal_detections++;
		double startTime = rtime(1);
		lit_counts.growTo(nVars() * 2);
		for (int i = 0; i < lit_counts.size(); i++) {
			lit_counts[i].seen = false;
		}
		assert(decisionLevel() == 0);
		
		//instead of counting the number of occurrences, just check if there are any occurences
		for (Lit l : trail) {
			lit_counts[toInt(l)].seen = true;
		}
		
		for (CRef cr : clauses) {
			Clause & c = ca[cr];
			for (Lit l : c) {
				lit_counts[toInt(l)].seen = true;
			}
		}
		for (CRef cr : learnts) {
			Clause & c = ca[cr];
			for (Lit l : c) {
				lit_counts[toInt(l)].seen = true;
			}
		}
		for (Var v = 0; v < nVars(); v++) {
			
			Lit l = mkLit(v, false);
			if (lit_counts[toInt(l)].seen && !lit_counts[toInt(l)].occurs) {
				stats_pure_lits--;
				lit_counts[toInt(l)].occurs = true;
				if (hasTheory(v)) {
					stats_pure_theory_lits--;
					theories[getTheoryID(v)]->setLiteralOccurs(getTheoryLit(l), true);
				}
			}
			if (lit_counts[toInt(~l)].seen && !lit_counts[toInt(~l)].occurs) {
				stats_pure_lits--;
				lit_counts[toInt(~l)].occurs = true;
				if (hasTheory(v)) {
					stats_pure_theory_lits--;
					theories[getTheoryID(v)]->setLiteralOccurs(getTheoryLit(~l), true);
				}
			}
			
			if (lit_counts[toInt(l)].occurs && !lit_counts[toInt(l)].seen) {
				lit_counts[toInt(l)].occurs = false;
				stats_pure_lits++;
				if (hasTheory(v)) {
					stats_pure_theory_lits++;
					theories[getTheoryID(v)]->setLiteralOccurs(getTheoryLit(l), false);
					//setPolarity(v,false);
				} else {
					//we can safely assign this now...
					if (!lit_counts[toInt(~l)].seen) {
						/*	if(decision[v])
						 setDecisionVar(v,false);*/
					} else {
						if (value(l) == l_Undef) {
							//uncheckedEnqueue(~l);
						}
					}
				}
			}
			if (lit_counts[toInt(~l)].occurs && !lit_counts[toInt(~l)].seen) {
				lit_counts[toInt(~l)].occurs = false;
				stats_pure_lits++;
				if (hasTheory(v)) {
					stats_pure_theory_lits++;
					theories[getTheoryID(v)]->setLiteralOccurs(getTheoryLit(~l), false);
					//setPolarity(v,true);
					if (!lit_counts[toInt(l)].occurs) {
						//setDecisionVar(v,false);//If v is pure in both polarities, and is a theory var, then just don't assign it at all - it is unconstrained.
						//This _should_ be safe to do, if the theory semantics make sense... although there are some clause learning schemes that introduce previosuly unused literals that might break with this...
					}
				} else {
					//we can safely assign this now...
					if (!lit_counts[toInt(l)].seen) {
						/*		if(decision[v])
						 setDecisionVar(v,false);*/
					} else if (value(l) == l_Undef) {
						//uncheckedEnqueue(l);
					}
				}
			}
			
		}
		assert(stats_pure_lits <= nVars() * 2);
		stats_pure_lit_time += rtime(1) - startTime;
		//}
	}
	
	return true;
}

bool Solver::addConflictClause(vec<Lit> & ps, CRef & confl_out, bool permanent) {
	dbg_check(theory_conflict);
	sort(ps);
	Lit p;
	int i, j;
	for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
		if (((value(ps[i]) == l_True && level(var(ps[i])) == 0)) || ps[i] == ~p)
			return true;
		else if ((value(ps[i]) != l_False || level(var(ps[i])) != 0) && ps[i] != p)
			ps[j++] = p = ps[i];
	ps.shrink(i - j);
	
	confl_out = CRef_Undef;
	if (ps.size() == 0) {
		ok = false;
		cancelUntil(0);
		return false;
	} else if (ps.size() == 1) {
		cancelUntil(0);
		assert(var(ps[0]) < nVars());
		dbg_check_propagation(ps[0]);
		if (!enqueue(ps[0])) {
			ok = false;
			return false;
		}
		
	} else {
		//find the highest level in the conflict (should be the current decision level, but we won't require that)
		bool conflicting = true;
		int nfalse = 0;
		int max_lev = 0;
		bool satisfied = false;
		int notFalsePos1 = -1;
		int notFalsePos2 = -1;
		for (int j = 0; j < ps.size(); j++) {
			assert(var(ps[j]) < nVars());
			//assert(value(ps[j])==l_False);
			if (value(ps[j]) == l_False) {
				nfalse++;
			} else {
				conflicting = false;
				if (value(ps[j]) == l_True)
					satisfied = true;
				if (notFalsePos1 < 0)
					notFalsePos1 = j;
				else if (notFalsePos2 < 0) {
					notFalsePos2 = j;
				}
			}
			if (value(ps[j]) != l_Undef) {
				int l = level(var(ps[j]));
				if (l > max_lev) {
					max_lev = l;
				}
			}
		}
		if (!conflicting) {
			assert(notFalsePos1 >= 0);
			if (notFalsePos1 >= 0 && notFalsePos2 >= 0) {
				assert(notFalsePos1 != notFalsePos2);
				if (notFalsePos1 == 1) {
					std::swap(ps[0], ps[notFalsePos2]);
				} else {
					std::swap(ps[0], ps[notFalsePos1]);
					std::swap(ps[1], ps[notFalsePos2]);
				}
			} else {
				std::swap(ps[0], ps[notFalsePos1]);
			}
			assert(value(ps[0])!=l_False);
			if (notFalsePos2 >= 0) {
				assert(value(ps[1])!=l_False);
			}
		}
		
#ifndef NDEBUG
		if (conflicting) {
			for (int j = 0; j < ps.size(); j++)
				assert(value(ps[j])==l_False);
		}
#endif
		if (!permanent && ps.size() > opt_temporary_theory_conflicts) {
			if (tmp_clause == CRef_Undef) {
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else if (tmp_clause_sz < ps.size()) {
				ca[tmp_clause].mark(1);							//is this needed?
				ca.free(tmp_clause);
				tmp_clause = ca.alloc(ps, false);
				tmp_clause_sz = ps.size();
			} else {
				Clause & c = ca[tmp_clause];
				assert(tmp_clause_sz >= ps.size());
				assert(tmp_clause_sz >= c.size());
				c.grow(tmp_clause_sz - c.size());
				for (int i = 0; i < ps.size(); i++) {
					c[i] = ps[i];
				}
				c.shrink(c.size() - ps.size());
			}
			
			confl_out = tmp_clause;
		} else {
			
			CRef cr = ca.alloc(ps, !permanent && !opt_permanent_theory_conflicts);
			ca[cr].setFromTheory(true);
			if (permanent || opt_permanent_theory_conflicts)
				clauses.push(cr);
			else {
				learnts.push(cr);
				if (--learntsize_adjust_cnt <= 0) {
					learntsize_adjust_confl *= learntsize_adjust_inc;
					learntsize_adjust_cnt = (int) learntsize_adjust_confl;
					max_learnts *= learntsize_inc;
					
					if (verbosity >= 1)
						printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %ld removed |\n", (int) conflicts,
								(int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
								(int) clauses_literals, (int) max_learnts, nLearnts(),
								(double) learnts_literals / nLearnts(), stats_removed_clauses);
				}
				
			}
			
			attachClause(cr);
			confl_out = cr;
		}
		if (!satisfied)
			cancelUntil(max_lev);
		
		if (!satisfied && nfalse == ps.size() - 1) {
			assert(value(ps[0])!=l_False);
			assert(value(ps[1])==l_False);
			if (value(ps[0]) == l_Undef) {
				uncheckedEnqueue(ps[0], confl_out);
			}
			confl_out = CRef_Undef;
		}
		return !conflicting;
	}
	return true;
}

bool Solver::addDelayedClauses(CRef & conflict_out) {
	conflict_out = CRef_Undef;
	while (clauses_to_add.size() && ok) {
		if (!addConflictClause(clauses_to_add.last(), conflict_out, true)) {
			clauses_to_add.pop();
			return false;
		}
		clauses_to_add.pop();
	}
	return true;
}

/*_________________________________________________________________________________________________
 |
 |  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
 |  
 |  Description:
 |    Search for a model the specified number of conflicts. 
 |    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
 |  
 |  Output:
 |    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
 |    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
 |    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
 |________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts) {
	assert(ok);
	int backtrack_level;
	int conflictC = 0;
	vec<Lit> learnt_clause;
	
	starts++;
	//last_dec = var_Undef;
	for (;;) {
		static int iter = 0;
		if (++iter == 5) {
			int a = 1;
		}
		propagate: CRef confl = propagate();
		conflict: if (!okay() || (confl != CRef_Undef)) {
			// CONFLICT
			conflicts++;
			conflictC++;
			if (decisionLevel() == 0)
				return l_False;
			learnt_clause.clear();
			analyze(confl, learnt_clause, backtrack_level);
			
			cancelUntil(backtrack_level);
			
			//this is now slightly more complicated, if there are multiple lits implied by the super solver in the current decision level:
			//The learnt clause may not be asserting.
			
			if (learnt_clause.size() == 1) {
				uncheckedEnqueue(learnt_clause[0]);
			} else {
				CRef cr = ca.alloc(learnt_clause, true);
				learnts.push(cr);
				attachClause(cr);
				claBumpActivity(ca[cr]);
				
				if (value(learnt_clause[0]) == l_Undef) {
					uncheckedEnqueue(learnt_clause[0], cr);
				} else {
					assert(S);
					//this is _not_ an asserting clause, its a conflict that must be passed up to the super solver.
					analyzeFinal(cr, lit_Undef, conflict);
					
					varDecayActivity();
					claDecayActivity();
					return l_False;
				}
				
			}
			
			varDecayActivity();
			claDecayActivity();
			
			if (--learntsize_adjust_cnt <= 0) {
				learntsize_adjust_confl *= learntsize_adjust_inc;
				learntsize_adjust_cnt = (int) learntsize_adjust_confl;
				max_learnts *= learntsize_inc;
				
				if (verbosity >= 1)
					printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %ld removed |\n", (int) conflicts,
							(int) dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(),
							(int) clauses_literals, (int) max_learnts, nLearnts(),
							(double) learnts_literals / nLearnts(), stats_removed_clauses);
			}
			
		} else {
			assert(theory_queue.size() == 0);
			
			if (!addDelayedClauses(confl))
				goto conflict;
			
			if (opt_subsearch == 0 && decisionLevel() < initial_level) {
				return l_Undef;            //give up if we have backtracked past the super solvers decisions
			} else if (opt_subsearch == 2 && S && confl == CRef_Undef) {
				if (super_qhead < S->qhead && !propagateTheory(theory_conflict)) //re-enqueue any super solver decisions that we have backtracked past, and keep going
						{
					if (!addConflictClause(theory_conflict, confl))
						return l_False;
					goto conflict;
				}
			}
			
			// NO CONFLICT
			if ((nof_conflicts >= 0 && conflictC >= nof_conflicts) || !withinBudget()) {
				// Reached bound on number of conflicts:
				progress_estimate = progressEstimate();
				cancelUntil(initial_level);
				return l_Undef;
			}
			
			// Simplify the set of problem clauses:
			if (decisionLevel() == 0 && !simplify())
				return l_False;
			
			if (learnts.size() - nAssigns() >= max_learnts)
				// Reduce the set of learnt clauses:
				reduceDB();
			
			Lit next = lit_Undef;
			while (decisionLevel() < assumptions.size()) {
				// Perform user provided assumption:
				Lit p = assumptions[decisionLevel()];
				if (value(p) == l_True) {
					// Dummy decision level:
					newDecisionLevel();
				} else if (value(p) == l_False) {
					analyzeFinal(~p, conflict);
					dbg_check(conflict);
					dbg_check_propagation(~p);
					return l_False;
				} else {
					next = p;
					break;
				}
			}
			
			if (opt_decide_theories && next == lit_Undef && drand(random_seed) < opt_random_theory_freq) {
				/**
				 * Give the theory solvers a chance to make decisions
				 */
				for (int i = 0; i < decidable_theories.size() && next == lit_Undef; i++) {
					
					next = decidable_theories[i]->decideTheory();
					
				}
			}
			
			if (next == lit_Undef) {
				// New variable decision:
				decisions++;
				next = pickBranchLit();
				// int p = priority[var(next)];
				
				if (next == lit_Undef) {
					
					//solve theories if this solver is completely assigned
					for (int i = 0; i < theories.size(); i++) {
						if (opt_subsearch == 3 && track_min_level < initial_level)
							continue; //Disable attempting to solve sub-solvers if we've backtracked past the super solver's decision level
							
						if (!theories[i]->solveTheory(theory_conflict)) {
							if (!addConflictClause(theory_conflict, confl)) {
								goto conflict;
							} else {
								goto propagate;
							}
						}
						//If propagating one of the sub theories caused this solver to backtrack, then go back to propagation
						if (qhead < trail.size() || nUnassignedVars() > 0)
							goto propagate;
					}
					
					// Model found:
					return l_True;
				}
			}
			//last_dec = var(next);
			// Increase decision level and enqueue 'next'
			newDecisionLevel();
			uncheckedEnqueue(next);
		}
	}
	//Unreachable
	return l_Undef;
}

double Solver::progressEstimate() const {
	double progress = 0;
	double F = 1.0 / nVars();
	
	for (int i = 0; i <= decisionLevel(); i++) {
		int beg = i == 0 ? 0 : trail_lim[i - 1];
		int end = i == decisionLevel() ? trail.size() : trail_lim[i];
		progress += pow(F, i) * (end - beg);
	}
	
	return progress / nVars();
}

/*
 Finite subsequences of the Luby-sequence:

 0: 1
 1: 1 1 2
 2: 1 1 2 1 1 2 4
 3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
 ...


 */

static double luby(double y, int x) {
	
	// Find the finite subsequence that contains index 'x', and the
	// size of that subsequence:
	int size, seq;
	for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1)
		;
	
	while (size - 1 != x) {
		size = (size - 1) >> 1;
		seq--;
		x = x % size;
	}
	
	return pow(y, seq);
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_() {
	cancelUntil(0);
	model.clear();
	conflict.clear();
	if (!ok)
		return l_False;
	
	solves++;
	
	max_learnts = nClauses() * learntsize_factor;
	learntsize_adjust_confl = learntsize_adjust_start_confl;
	learntsize_adjust_cnt = (int) learntsize_adjust_confl;
	lbool status = l_Undef;
	
	if (verbosity >= 1 && !printed_header) {
		//on repeated calls, don't print the header again
		printed_header = true;
		printf("============================[ Search Statistics ]==============================\n");
		printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
		printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
		printf("===============================================================================\n");
	}
	initial_level = 0;
	track_min_level = 0;
	// Search:
	int curr_restarts = 0;
	while (status == l_Undef) {
		double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
		
		if (opt_rnd_phase) {
			for (int i = 0; i < nVars(); i++)
				polarity[i] = irand(random_seed, 1);
		}
		
		status = search(rest_base * restart_first);
		if (!withinBudget())
			break;
		curr_restarts++;
		if (opt_rnd_restart && status == l_Undef) {
			
			for (int i = 0; i < nVars(); i++) {
				activity[i] = drand(random_seed) * 0.00001;
			}
			rebuildOrderHeap();
			
		}
	}
	
	if (status == l_True) {
		// Extend & copy model:
		model.growTo(nVars());
		for (int i = 0; i < nVars(); i++)
			model[i] = value(i);
		
		if (opt_check_solution && theories.size()) {
			Theory * t = theories[0];
			
			if (!t->check_solved()) {
				fprintf(stderr, "Error! Solution doesn't satisfy theory properties!\n");
				fflush(stderr);
				exit(4);
			}
		}
#ifdef DEBUG_SOLVER
		if(dbg_solver)
		assert(dbg_solver->solve(assumptions));
#endif
		
	} else if (status == l_False && conflict.size() == 0) {
		ok = false;
#ifdef DEBUG_SOLVER
		if(dbg_solver)
		assert(!dbg_solver->solve());
#endif
	} else if (status == l_False) {
		assert(ok);
#ifdef DEBUG_SOLVER
		if(dbg_solver)
		assert(!dbg_solver->solve(assumptions));
#endif
	}
	
	assumptions.clear();
	return status;
}

bool Solver::solveTheory(vec<Lit> & conflict_out) {
	initial_level = decisionLevel();
	track_min_level = initial_level;
	lbool status = l_Undef;
	// Search:
	int curr_restarts = 0;
	conflict.clear();
	conflict_out.clear();
	while (status == l_Undef) {
		double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
		status = search(rest_base * restart_first);
		if (!withinBudget() || (opt_subsearch == 0 && track_min_level < initial_level))
			break;
		curr_restarts++;
	}
	cancelUntil(track_min_level);
	if (track_min_level < initial_level) {
		S->cancelUntil(track_min_level);
	}
	if (!ok) {
		return false;
	}
	if (conflict.size()) {
		toSuper(conflict, conflict_out);
		interpolant.push();
		conflict_out.copyTo(interpolant.last());    //Add this clause into the interpolant vector
		return false;
	}
	
	return propagateTheory(conflict_out);
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max) {
	if (map.size() <= x || map[x] == -1) {
		map.growTo(x + 1, -1);
		map[x] = max++;
	}
	return map[x];
}

void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max) {
	if (satisfied(c))
		return;
	
	for (int i = 0; i < c.size(); i++)
		if (value(c[i]) != l_False)
			fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max) + 1);
	fprintf(f, "0\n");
}

void Solver::toDimacs(const char *file, const vec<Lit>& assumps) {
	FILE* f = fopen(file, "wr");
	if (f == NULL)
		fprintf(stderr, "could not open file %s\n", file), exit(1);
	toDimacs(f, assumps);
	fclose(f);
}

void Solver::toDimacs(FILE* f, const vec<Lit>& assumps) {
	// Handle case when solver is in contradictory state:
	if (!ok) {
		fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
		return;
	}
	
	vec<Var> map;
	Var max = 0;
	
	// Cannot use removeClauses here because it is not safe
	// to deallocate them at this point. Could be improved.
	int cnt = 0;
	for (int i = 0; i < clauses.size(); i++)
		if (!satisfied(ca[clauses[i]]))
			cnt++;
	
	for (int i = 0; i < clauses.size(); i++)
		if (!satisfied(ca[clauses[i]])) {
			Clause& c = ca[clauses[i]];
			for (int j = 0; j < c.size(); j++)
				if (value(c[j]) != l_False)
					mapVar(var(c[j]), map, max);
		}
	
	// Assumptions are added as unit clauses:
	cnt += assumptions.size();
	
	fprintf(f, "p cnf %d %d\n", max, cnt);
	
	for (int i = 0; i < assumptions.size(); i++) {
		assert(value(assumptions[i]) != l_False);
		fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max) + 1);
	}
	
	for (int i = 0; i < clauses.size(); i++)
		toDimacs(f, ca[clauses[i]], map, max);
	
	if (verbosity > 0)
		printf("Wrote %d clauses with %d variables.\n", cnt, max);
}

//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to) {
	if (tmp_clause != CRef_Undef) {
		ca[tmp_clause].mark(1);    //is this needed?
		ca.free(tmp_clause);
		tmp_clause = CRef_Undef;
	}
	//Re-allocate the 'theory markers'
	vec<int> marker_theory_tmp;
	
	/*	for(int i = 0;i<markers.size();i++){
	 CRef old_cr = markers[i];
	 assert(old_cr!=CRef_Undef);
	 int old_theory = getTheory(old_cr);
	 assert(old_theory>=0);
	 CRef cr=to.makeMarkerReference();
	 assert(cr==markers[i]);//these should be identical in the current implementation
	 markers[i]=cr;
	 int index = CRef_Undef-cr - 1;
	 marker_theory_tmp.growTo(index+1,-1 );
	 marker_theory_tmp[index] = old_theory;

	 assert( marker_theory_tmp[index] == old_theory);

	 assert(markers[ marker_theory_tmp[index]] == cr);
	 }
	 marker_theory_tmp.copyTo(marker_theory);*/
	// All watchers:
	//
	// for (int i = 0; i < watches.size(); i++)
	watches.cleanAll();
	for (int v = 0; v < nVars(); v++)
		for (int s = 0; s < 2; s++) {
			Lit p = mkLit(v, s);
			// printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
			vec<Watcher>& ws = watches[p];
			for (int j = 0; j < ws.size(); j++)
				ca.reloc(ws[j].cref, to);
		}
	
	// All reasons:
	//
	for (int i = 0; i < trail.size(); i++) {
		Var v = var(trail[i]);
		
		if (ca.isClause(reason(v)) && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
			ca.reloc(vardata[v].reason, to);
	}
	
	// All learnt:
	//
	for (int i = 0; i < learnts.size(); i++)
		ca.reloc(learnts[i], to);
	
	// All original:
	//
	for (int i = 0; i < clauses.size(); i++)
		ca.reloc(clauses[i], to);
}

void Solver::garbageCollect() {
	// Initialize the next region to a size corresponding to the estimated utilization degree. This
	// is not precise but should avoid some unnecessary reallocations for the new region:
	ClauseAllocator to(ca.size() - ca.wasted());
	
	relocAll(to);
	if (verbosity >= 2)
		printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n",
				ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
	to.moveTo(ca);
}

