/*
 * CTLSolver.h
 *
 *  Created on: Mar 26, 2015
 *      Author: tobias k
 *
 *      TODO Implement changes outlined in fbb8e4d. Basically, only have a constant
 *      number of Bitsets, instead of creating and deleting them on the fly.
 */

#ifndef CTL_CTLSOLVER_H_
#define CTL_CTLSOLVER_H_

#include <vector>
#include "mtl/Vec.h"
#include "mtl/Bitset.h"
#include <algorithm>
#include <cassert>
#include "alg/NFATypes.h"
#include "dgl/DynamicGraph.h"
#include "DynamicKripke.h"
#include "CTLFormula.h"
#include "CTLTarjanSCC.h"

using namespace dgl;
namespace Monosat {

class CTLSolver {
public:
	DynamicKripke* k;
	DynamicKripke* otherk;
	DynamicKripke* originalk;
	DynamicKripke* originalotherk;
	DynamicKripke* tmpk;
	vec<Bitset*> bitsets;
	vec<bool> bitsetsAvail;

	bool isover; // indicates that this is the solver for the Overapprox., not the Underapprox.

	int bitsetcounter;
	int bitsetsmax;
	int id;

	// Used for Tarjan's SCC algorithm
	CTLTarjansSCC<int>* tarjan;
	std::vector<int> q;
	std::vector<int> comp;
	std::vector<int> innerFair;
	std::vector<char> seen;

	CTLSolver(int myid, DynamicKripke& myk, DynamicKripke& myotherk, bool myisover){
		id = myid;
		k = &myk;
		otherk = &myotherk;
		originalk = &myk;
		originalotherk = &myotherk;
		tmpk = NULL;
		isover = myisover;
		bitsetcounter = 0;
		bitsetsmax = 0;
		tarjan = new CTLTarjansSCC<int>(myk.g, myk);

	};
	~CTLSolver() {};

	int getFreshBitset() {
		int i = getDirtyBitset();
		bitsets[i]->memset(false);
		return i;
	}

	int getDirtyBitset() {
		bitsetcounter++;
		/*
		printf("%d getBitset", isover);
		for (int j = 0; j < bitsetsAvail.size(); j++) {
			printf("%d: %d | ", j, bitsetsAvail[j]);
		}
		*/
		int i;

		// Look for an available bitset in vector
		for (i = 0; i < bitsetsAvail.size(); i++) {
			if (bitsetsAvail[i]) {
				bitsetsAvail[i] = false; // set to false

				/*
				printf("getFreshBitset() returns existing bitset %d  ", i);
				for (int j = 0; j < bitsetsAvail.size(); j++) {
					printf("%d: %d | ", j, bitsetsAvail[j]);
				}
				printf("\n");
				*/
				return i;
			}
		}
		// Nothing free, create new one
		i = bitsetsAvail.size();
		bitsets.growTo(i + 1);
		bitsetsAvail.growTo(i + 1);
		bitsetsAvail[i] = false; // set to false
		bitsets[i] = new Bitset(k->states());
		bitsetsmax = i;

		/*
		printf("getFreshBitset() creates new bitset %d  ", i);
		for (int j = 0; j < bitsetsAvail.size(); j++) {
			printf("%d: %d | ", j, bitsetsAvail[j]);
		}
		printf("\n");
		*/

		return i;
	}
	void freeBitset(int i) {
		//printf("%d freeBitset %d: ", isover, i);
		//for (int j = 0; j < bitsetsAvail.size(); j++) {
		//	printf("%d: %d | ", j, bitsetsAvail[j]);
		//}

		assert(i < bitsetsAvail.size());
		if (!bitsetsAvail[i])
			bitsetsAvail[i] = true; // set to false
		else
			throw std::runtime_error("Trying to free a bitset that is available!");

		//printf("freeBitset cleaned bitset %d  ", i);
		//for (int j = 0; j < bitsetsAvail.size(); j++) {
		//	printf("%d: %d | ", j, bitsetsAvail[j]);
		//}
		//printf("\n");

	}

	void reset() {
		resetBitsets();
		resetSwap();
	}

	void resetBitsets() {
		//printf("%d getBitset", isover);
		//printf("Resetting Bitsets\n");
		for (int i = 0; i < bitsetsAvail.size(); i++) {
			if (!bitsetsAvail[i]) {
				bitsetsAvail.operator [](i) = true; // set to false
			}
		}
		//printf("freeBitset cleaned all bitsets  ");
		//for (int j = 0; j < bitsetsAvail.size(); j++) {
		//	printf("%d: %d | ", j, bitsetsAvail[j]);
		//}
		//printf("\n");

	}

	void swapKripkes() {
		tmpk = k;
		k = otherk;
		otherk = tmpk;
	}

public:
	void resetSwap() {
		k = originalk;
		otherk = originalotherk;
	}

	// Predecessors with respect to the Kripke structure's transition system.
	// Note that this creates a new Bitset, you might want to reuse an existing bitset and use the next function below

	int pre(int st) {
		int prest = getDirtyBitset(); // we clean later
		pre(st, prest);
		return prest;
	}
	void pre(int st, int prest) {
		bitsets[prest]->memset(false); // clean it, just in case
		for (int i=0; i<k->nEdgeIDs();i++) {
			if (bitsets[st]->operator [](k->getEdge(i).to) && k->transitionEnabled(i)) {
				bitsets[prest]->set(k->getEdge(i).from);
			}
		}
	}

	// Similar to solve(), but returns a bitset pointer// and resets the system beforhand.
	int solveFormula(CTLFormula& f) {
		//printf("%d: Solving ", isover);
		//printFormula(&f);
		//printf("\n");
		//reset();
		int i = solve(f);
		//printf("%d: Done solving\n", isover);
		return i;
	}

	// Main solve function.
	int solve(CTLFormula& f) {
		switch (f.op) {
		case ID : return solveID(f);
		case True : return solveTrue(f);
		case NEG : return solveNEG(f);
		case OR : return solveOR(f);
		case AND : return solveAND(f);
		case EX : return solveEX(f);
		case EG : return solveEG(f);
		case EF : return solveEF(f);
		case EW : return solveEW(f);
		case EU : return solveEU(f);
		case AX : return solveAX(f);
		case AG : return solveAG(f);
		case AF : return solveAF(f);
		case AW : return solveAW(f);
		case AU : return solveAU(f);
		default: throw std::runtime_error("CTLSolver.solve has encountered unsupported CTL operator!");
		}
	}

	int solveID(CTLFormula& f) {
		assert(f.op == ID);
		int st = getFreshBitset();
		for (int i=0; i<k->states();i++) {
			if (k->statelabel[i]->operator [](f.value))
				bitsets[st]->set(i);
		}
		return st;
	}

	int solveTrue(CTLFormula& f) {
		assert(f.op == True);
		int st = getDirtyBitset(); // dirty is OK, we clean it
		bitsets[st]->memset(true);
		return st;
	}

	int solveNEG(CTLFormula& f) {
		assert(f.op == NEG);

		// Be smart about double negations and cancel them out
		if (f.operand1->op == NEG) {
			return solve(*f.operand1->operand1);
		}

		swapKripkes();

		int st = solve(*f.operand1);
		bitsets[st]->NotSelf();

		swapKripkes();

		return st;
	}

	int solveOR(CTLFormula& f) {
		assert(f.op == OR);
		int st1 = solve(*f.operand1);
		int st2 = solve(*f.operand2);
		bitsets[st1]->Or(*bitsets[st2]);
		freeBitset(st2);
		return st1;
	}

	// OK, this is not strictly needed, since we have "negation" and "or", but for performance
	// reasons we implement it directly anyway.
	int solveAND(CTLFormula& f) {
		assert(f.op == AND);
		int st1 = solve(*f.operand1);
		int st2 = solve(*f.operand2);
		bitsets[st1]->And(*bitsets[st2]);
		freeBitset(st2);
		return st1;
	}

	int solveEX(CTLFormula& f) {
		assert(f.op == EX);
		int st = solve(*f.operand1);

		int prest = pre(st);
		freeBitset(st);

		return prest;

		/* Jeff says, ideally something like this
		Bitset prest(k->states());
		prest = std::move(*pre(*st));
		delete st;

		return prest;*/
	}

	/*
	 * This is where it gets interesting. We look for the largest solution of X = μ(p) ∩ pre(X).
	 * Luckily for us, we can simply compute the fixpoint
	 */
	int solveEG(CTLFormula& f) {
		assert(f.op == EG);

		if (f.fairnessConstraints.size() > 0)
			return solveEGwithFairness(f);
		// μ(p)
		int st = solve(*f.operand1);
		// Auxiliary bitsets
		int andst = getFreshBitset();
		int prest = getFreshBitset();

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(st, prest); // prest := pre(st)

			// μ(p) ∩ pre(X).
			bitsets[st]->And(*bitsets[prest], *bitsets[andst]); // andst := st ∩ prest

			if (bitsets[st]->Equiv(*bitsets[andst])) {
				freeBitset(st);
				freeBitset(prest);
				return andst; // fixpoint reached
			}
			bitsets[st]->copyFrom(*bitsets[andst]);
		}
	}

	int solveEGwithFairness(CTLFormula& f) {
		printf("Components: \n");
		DynamicGraph<int> & g = k->g;

		int innerBitset = solve(*f.operand1);
		tarjan->setInnerFormula(*bitsets[innerBitset]);
		printStateSet(*bitsets[innerBitset]);
		printf("\n");

		// This will end up being μ(f.operand1), restricted to states that are on an SCC which satisfies all fairness constraints
		int st = getFreshBitset();

		// A bitset over the set of fairness constraints, representing which fairness constraint is satisfied by a certain SCC.
		Bitset* fairnessConstraintsSat = new Bitset(f.fairnessConstraints.size());

		// Solve formulas of the fairness constraints yielding a set of bitsets, which tell in which states each fairness constraint is satisfied
		innerFair.clear();
		innerFair.resize(f.fairnessConstraints.size());
		for (int i = 0; i < innerFair.size(); i ++) {
			innerFair[i] = solve(*f.fairnessConstraints.operator [](i));
			printf("Fairness constraint and set of states in which it is true: ");
			printFormula(f.fairnessConstraints.operator [](i));
			printStateSet(*bitsets[innerFair[i]]);
		}


		seen.clear();
		seen.resize(g.nodes(),false);
		for (int i = 0; i < tarjan->numComponents();i++ ) {
			// Start out with one arbitrary element of the SCC
			int node = tarjan->getElement(i);
			// Find all other elements of the SCC
			// Makes sure that we skip over lone vertices that have no self loop.
			if (tarjan->getComponentSize(i) > 1 || suitable(k->getEdge(node, node),*bitsets[innerBitset])) {
				printf("Component %d:",i);
				assert(q.size()==0);
				comp.clear();
				q.push_back(node);
				fairnessConstraintsSat->memset(false);
				seen[node]=true;

				//do a DFS to recover the connected nodes
				while(q.size()){
					int u = q.back();
					q.pop_back();
					comp.push_back(u);
					// Go through all fairness constraints. If a fairness constraint is satisfied in this state, then mark it satisfied for the entire SCC in fairnessConstraintsSat
					for (int j = 0; j < innerFair.size(); j ++) {
						if (bitsets[innerFair[j]]->operator [](u)) {
							fairnessConstraintsSat->set(j);
							printf("\nSCC satisfies fairness constraint %d, because state %d satisfies it\n", j, u);
						}
					}
					printf(" %d",u);
					for(int j = 0;j<g.nIncident(u);j++){
						int edgeID = g.incident(u,j).id;
						int to = g.incident(u,j).node;
						if(suitable(edgeID,*bitsets[innerBitset]) && !seen[to]){//only traverse
							seen[to]=true;
							q.push_back(to);
						}
					}
				}
				bool fairSCC = true; // innocent until proven guilty
				for (int j=0; j<fairnessConstraintsSat->size(); j++) {
					if (!fairnessConstraintsSat->operator [](j)) {
						fairSCC = false;
						printf("\nUnfair SCC due to fairness constraint #%d\n", j);
					}
				}
				if (fairSCC) {
					for (int j=0; j<comp.size(); j++) {
						bitsets[st]->set(j);
					}
				}

				printf("\n");
			}
		}

		printf("\n");

		printf("Nodes belonging to an SCCs that we are backtracking from in EGFair: ");
		printStateSet(*bitsets[st]);
		printf("\n\n");

		// Auxiliary bitsets
		int andst = getFreshBitset();
		int prest = getFreshBitset();

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(st, prest); // prest := pre(st)

			// μ(p) ∩ pre(X).
			bitsets[st]->And(*bitsets[prest], *bitsets[andst]); // andst := st ∩ prest

			if (bitsets[st]->Equiv(*bitsets[andst])) {
				freeBitset(st);
				freeBitset(prest);
				for (int i = 0; i < innerFair.size(); i ++)
					freeBitset(innerFair[i]);
				delete fairnessConstraintsSat;

				printf("Solution to EGFair: ");
				printStateSet(*bitsets[andst]);
				printf("\n\n");
				return andst; // fixpoint reached
			}
			bitsets[st]->copyFrom(*bitsets[andst]);
		}


	}
	// An edge is suitable if its origin and destination satisfies the formula and it is enabled
	// This represents a graph restricted to vertices that satisfy the formula
	bool suitable(int edgeID, Bitset& inner) {
		return (k->edgeEnabled(edgeID)) && inner.operator [](k->getEdge(edgeID).from) && inner.operator [](k->getEdge(edgeID).to);
	}

	// X = μ(p) ∪ pre(X)
	int solveEF(CTLFormula& f) {
		assert(f.op == EF);
		// μ(p)
		int st = solve(*f.operand1);
		// Auxiliary bitsets
		int orst = getFreshBitset();
		int prest = getFreshBitset();

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(st, prest); // prest := pre(st)

			// μ(p) ∩ pre(X).
			bitsets[st]->Or(*bitsets[prest], *bitsets[orst]); // andst := st ∩ prest

			if (bitsets[st]->Equiv(*bitsets[orst])) {
				freeBitset(st);
				freeBitset(prest);
				return orst; // fixpoint reached
			}
			bitsets[st]->copyFrom(*bitsets[orst]);
		}
	}

	// φ 1 EW φ 2 ≡ EG φ 1 ∨ (φ 1 EU φ 2 )
	int solveEW(CTLFormula& f) {
		CTLFormula phi1 = {EG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EU, f.operand1, f.operand2, 0};
		CTLFormula phi3 = {OR, &phi1, &phi2, 0};
		return solve(phi3);
	}

	// X = μ(q) ∪ (μ(p) ∩ pre(X)).
	int solveEU(CTLFormula& f) {
		assert(f.op == EU);
		// μ(p)
		int st1 = solve(*f.operand1);
		int st2 = solve(*f.operand2);
		// Auxiliary bitsets
		int andst = getFreshBitset();
		int orst = getFreshBitset();
		int prest = getFreshBitset();
		int x = getFreshBitset();
		bitsets[x]->copyFrom(*bitsets[st2]); // start out with X as μ(q)

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(x, prest);

			// μ(p) ∩ pre(X).
			bitsets[st1]->And(*bitsets[prest], *bitsets[andst]); // andst := st1 ∩ prest

			// μ(q) ∪ (μ(p) ∩ pre(X)).
			bitsets[andst]->Or(*bitsets[st2], *bitsets[orst]); // orst := andset ∪ st2

			if (bitsets[x]->Equiv(*bitsets[orst])) {
				freeBitset(st1);
				freeBitset(st2);
				freeBitset(andst);
				freeBitset(x);
				freeBitset(prest);
				return orst; // fixpoint reached
			}

			bitsets[x]->copyFrom(*bitsets[orst]);
		}
	}

	// AX φ ≡ ¬ EX ¬φ
	int solveAX(CTLFormula& f) {
		CTLFormula phi1 = {NEG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EX, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, &phi2, NULL, 0};
		return solve(phi3);
	}

	// AF φ ≡ ¬ EG ¬φ
	int solveAF(CTLFormula& f) {
		CTLFormula phi1 = {NEG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EG, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, &phi2, NULL, 0};
		return solve(phi3);
	}

	// AG φ ≡ ¬ EF ¬φ
	int solveAG(CTLFormula& f) {
		CTLFormula phi1 = {NEG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EF, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, &phi2, NULL, 0};
		return solve(phi3);
	}

	// φ 1 AW φ 2 ≡ ¬(¬φ 2 EU ¬(φ 1 ∨ φ 2 ))
	int solveAW(CTLFormula& f) {
		CTLFormula phi1 = {OR, f.operand1, f.operand2, 0};
		CTLFormula phi2 = {NEG, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, f.operand2, NULL, 0};
		CTLFormula phi4 = {EU, &phi3, &phi2, 0};
		CTLFormula phi5 = {NEG, &phi4, NULL, 0};
		return solve(phi5);
	}
	// φ 1 AU φ 2 ≡ AF φ 2 ∧ (φ 1 AW φ 2 ))
	int solveAU(CTLFormula& f) {
		CTLFormula phi1 = {AF, f.operand2, NULL, 0};
		CTLFormula phi2 = {AW, f.operand1, f.operand2, 0};
		CTLFormula phi3 = {AND, &phi1, &phi2, 0};
		return solve(phi3);
	}

	void funwithctl() {
		CTLFormula a {ID, NULL, NULL, 0};
		CTLFormula b {ID, NULL, NULL, 1};
		CTLFormula c {ID, NULL, NULL, 2};
		CTLFormula EXa {EX, &a, NULL, 0};
		CTLFormula NEGEXa {NEG, &EXa, NULL, 0};
		CTLFormula EUEXab {EU, &EXa, &b, 0};
		CTLFormula complicated {AF, &EUEXab, NULL, 0};
		CTLFormula NEGEXaORb {OR, &NEGEXa, &b, 0};
		CTLFormula NEGEXaORc {OR, &NEGEXa, &c, 0};
		CTLFormula EGb {EG, &b, NULL, 0};
		CTLFormula EUbc {EU, &b, &c, 0};
		CTLFormula AXb {AX, &b, NULL, 0};
		CTLFormula AFc {AF, &c, NULL, 0};
		CTLFormula AUbc {AU, &b, &c, 0};
		CTLFormula AUac {AU, &a, &c, 0};
		CTLFormula AWbc {AW, &b, &c, 0};
		CTLFormula AUbcOREGb {OR, &AUbc, &EGb, 0};
		CTLFormula AXc {AX, &c, NULL, 0};
		CTLFormula bAUAXc {AU, &b, &AXc, 0};

		//printFormula(complicated); printf("\n");

		Bitset* foo;

		/*
		foo = solve(EXa);
		printStateSet(*foo);

		foo = solve(NEGEXa);
		printFormula(NEGEXa); printStateSet(*foo);

		foo = solve(NEGEXaORb);
		printFormula(NEGEXaORb); printStateSet(*foo);

		foo = solve(NEGEXaORc);
		printFormula(NEGEXaORc); printStateSet(*foo);

		foo = solve(EGb);
		printFormula(EGb); printStateSet(*foo);

		foo = solve(EUbc);
		printFormula(EUbc); printStateSet(*foo);

		printf("All-quantified formulas:\n");

		foo = solve(AXb);
		printFormula(AXb); printStateSet(*foo);
*/

/*
		foo = solve(AUbc);
		printFormula(AUbc); printStateSet(*foo);
		foo = solve(AWbc);
		printFormula(AWbc); printStateSet(*foo);
		foo = solve(AUac);
		printFormula(AUac); printStateSet(*foo);
		foo = solve(AFc);
		printFormula(AFc); printStateSet(*foo);
		*/

		/*
		foo = solve(AUbcOREGb);
		printFormula(AUbcOREGb); printStateSet(*foo);

		foo = solve(AWbc);
		printFormula(AWbc); printStateSet(*foo);

		foo = solve(AXc);
		printFormula(AXc); printStateSet(*foo);

		foo = solve(bAUAXc);
		printFormula(bAUAXc); printStateSet(*foo);
		 */

		delete foo;
	}

public:
	void printStateSet(Bitset& s) {
		if(opt_verb>1){
			printf("{");
			for (int i=0; i<s.size(); i++) {
				if (s[i]) {
					printf("%d, ", i);
				}
			}
			printf("}\n");
		}
	}

};
};

#endif /* CTL_CTLSOLVER_H_ */

