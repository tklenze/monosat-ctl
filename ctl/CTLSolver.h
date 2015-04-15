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
using namespace dgl;
namespace Monosat {

class CTLSolver {
public:
	DynamicKripke k;
	CTLSolver(DynamicKripke myk) {
		k = myk;
	};
	~CTLSolver() {};

	// Predecessors with respect to the Kripke structure's transition system.
	// Note that this creates a new Bitset, you might want to reuse an existing bitset and use the next function below
	Bitset* pre(Bitset& st) {
		Bitset *prest = new Bitset(k.states());
		for (int i=0; i<k.nEdgeIDs();i++) {
			if (st[k.getEdge(i).to]) {
				prest->set(k.getEdge(i).from);
			}
		}
		return prest;
	}

	// Use st as input and store pre(st) in out. Out does not need to be clean
	void pre(Bitset& st, Bitset& out) {
		out.memset(false);
		for (int i=0; i<k.nEdgeIDs();i++) {
			if (st[k.getEdge(i).to]) {
				out.set(k.getEdge(i).from);
			}
		}
	}

	// Main solve function.
	Bitset* solve(CTLFormula& f) {
		switch (f.op) {
		case ID : return solveID(f);
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
		default : return NULL;
		}
	}

	Bitset* solveID(CTLFormula& f) {
		assert(f.op == ID);
		Bitset *st = new Bitset(k.states());
		for (int i=0; i<k.states();i++) {
			if (k.statelabel[i][f.value])
				st->set(i);
		}
		return st;
	}

	// FIXME it's not optimal that we create a new bitset every time we negate
	Bitset* solveNEG(CTLFormula& f) {
		assert(f.op == NEG);

		// Be smart about double negations and cancel them out
		if (f.operand1->op == NEG) {
			return solve(*f.operand1->operand1);
		}

		Bitset *st = solve(*f.operand1);
		Bitset *negst = new Bitset(k.states());
		st->Not(*negst);
		delete st; // TODO does this actually work? Please confirm
		return negst;
	}

	Bitset* solveOR(CTLFormula& f) {
		assert(f.op == OR);
		Bitset *st1 = solve(*f.operand1);
		Bitset *st2 = solve(*f.operand2);
		st1->Or(*st2);
		delete st2;
		return st1;
	}

	// OK, this is not strictly needed, since we have "negation" and "or", but for performance
	// reasons we implement it directly anyway.
	Bitset* solveAND(CTLFormula& f) {
		assert(f.op == AND);
		Bitset *st1 = solve(*f.operand1);
		Bitset *st2 = solve(*f.operand2);
		st1->And(*st2);
		delete st2;
		return st1;
	}

	Bitset* solveEX(CTLFormula& f) {
		assert(f.op == EX);
		Bitset *st = solve(*f.operand1);
		return pre(*st);
	}

	/*
	 * This is where it gets interesting. We look for the largest solution of X = μ(p) ∩ pre(X).
	 * Luckily for us, we can simply compute the fixpoint
	 */
	Bitset* solveEG(CTLFormula& f) {
		assert(f.op == EG);
		// μ(p)
		Bitset *st = solve(*f.operand1);
		// Auxiliary bitsets
		Bitset *andst = new Bitset(k.states());
		Bitset *prest = new Bitset(k.states());

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(*st, *prest); // prest := pre(st)

			// μ(p) ∩ pre(X).
			st->And(*prest, *andst); // andst := st ∩ prest

			if (st->Equiv(*andst)) {
				delete st;
				delete prest;
				return andst; // fixpoint reached
			}
			st->copyFrom(*andst);
		}
	}

	// X = μ(p) ∪ pre(X)
	Bitset* solveEF(CTLFormula& f) {
		assert(f.op == EF);
		// μ(p)
		Bitset *st = solve(*f.operand1);
		// Auxiliary bitsets
		Bitset *orst = new Bitset(k.states());
		Bitset *prest = new Bitset(k.states());

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(*st, *prest); // prest := pre(st)

			// μ(p) ∩ pre(X).
			st->Or(*prest, *orst); // andst := st ∩ prest

			if (st->Equiv(*orst)) {
				delete st;
				delete prest;
				return orst; // fixpoint reached
			}
			st->copyFrom(*orst);
		}
	}

	// φ 1 EW φ 2 ≡ EG φ 1 ∨ (φ 1 EU φ 2 )
	Bitset* solveEW(CTLFormula& f) {
		CTLFormula phi1 = {EG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EU, f.operand1, f.operand2, 0};
		CTLFormula phi3 = {OR, &phi1, &phi2, 0};
		return solve(phi3);
	}

	// X = μ(q) ∪ (μ(p) ∩ pre(X)).
	Bitset* solveEU(CTLFormula& f) {
		assert(f.op == EU);
		// μ(p)
		Bitset *st1 = solve(*f.operand1);
		Bitset *st2 = solve(*f.operand2);
		// Auxiliary bitsets
		Bitset *andst = new Bitset(k.states());
		Bitset *orst = new Bitset(k.states());
		Bitset *prest = new Bitset(k.states());
		Bitset *x = new Bitset(k.states());
		x->copyFrom(*st2); // start out with X as μ(q)

		while (true) { // fixpoint guaranteed to exist, therefore this will terminate
			// pre(X)
			pre(*x, *prest);

			// μ(p) ∩ pre(X).
			st1->And(*prest, *andst); // andst := st1 ∩ prest

			// μ(q) ∪ (μ(p) ∩ pre(X)).
			andst->Or(*st2, *orst); // orst := andset ∪ st2

			if (x->Equiv(*orst)) {
				delete st1;
				delete st2;
				delete andst;
				delete x;
				return orst; // fixpoint reached
			}

			x->copyFrom(*orst);
		}
	}

	// AX φ ≡ ¬ EX ¬φ
	Bitset* solveAX(CTLFormula& f) {
		CTLFormula phi1 = {NEG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EX, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, &phi2, NULL, 0};
		return solve(phi3);
	}

	// AF φ ≡ ¬ EG ¬φ
	Bitset* solveAF(CTLFormula& f) {
		CTLFormula phi1 = {NEG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EG, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, &phi2, NULL, 0};
		return solve(phi3);
	}

	// AG φ ≡ ¬ EF ¬φ
	Bitset* solveAG(CTLFormula& f) {
		CTLFormula phi1 = {NEG, f.operand1, NULL, 0};
		CTLFormula phi2 = {EF, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, &phi2, NULL, 0};
		return solve(phi3);
	}

	// φ 1 AW φ 2 ≡ ¬(¬φ 2 EU ¬(φ 1 ∨ φ 2 ))
	Bitset* solveAW(CTLFormula& f) {
		CTLFormula phi1 = {OR, f.operand1, f.operand2, 0};
		CTLFormula phi2 = {NEG, &phi1, NULL, 0};
		CTLFormula phi3 = {NEG, f.operand2, NULL, 0};
		CTLFormula phi4 = {EU, &phi3, &phi2, 0};
		CTLFormula phi5 = {NEG, &phi4, NULL, 0};
		return solve(phi5);
	}
	// φ 1 AU φ 2 ≡ AF φ 2 ∧ (φ 1 AW φ 2 ))
	Bitset* solveAU(CTLFormula& f) {
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


		foo = solve(AUbc);
		printFormula(AUbc); printStateSet(*foo);
		foo = solve(AWbc);
		printFormula(AWbc); printStateSet(*foo);
		foo = solve(AUac);
		printFormula(AUac); printStateSet(*foo);
		foo = solve(AFc);
		printFormula(AFc); printStateSet(*foo);

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

	void printStateSet(Bitset& s) {
		printf("{");
		for (int i=0; i<s.size(); i++) {
			if (s[i]) {
				printf("%d, ", i);
			}
		}
		printf("}\n");
	}

	void printFormula(CTLFormula &f) {
		switch (f.op) {
		case ID : printf("%d", f.value); break;
		case NEG: printf("not "); printFormula(*f.operand1); break;
		case OR : printf("("); printFormula(*f.operand1); printf(" or "); printFormula(*f.operand2); printf(")"); break;
		case AND : printf("("); printFormula(*f.operand1); printf(" and "); printFormula(*f.operand2); printf(")"); break;

		case EX : printf("EX "); printFormula(*f.operand1); break;
		case EF : printf("EF "); printFormula(*f.operand1); break;
		case EG : printf("EG "); printFormula(*f.operand1); break;
		case EW : printf("("); printFormula(*f.operand1); printf(" EW "); printFormula(*f.operand2); printf(")"); break;
		case EU : printf("("); printFormula(*f.operand1); printf(" EU "); printFormula(*f.operand2); printf(")"); break;

		case AX : printf("AX "); printFormula(*f.operand1); break;
		case AF : printf("AF "); printFormula(*f.operand1); break;
		case AG : printf("AG "); printFormula(*f.operand1); break;
		case AW : printf("("); printFormula(*f.operand1); printf(" AW "); printFormula(*f.operand2); printf(")"); break;
		case AU : printf("("); printFormula(*f.operand1); printf(" AU "); printFormula(*f.operand2); printf(")"); break;

		default : printf("bar");
		}
	}
};
};

#endif /* CTL_CTLSOLVER_H_ */

