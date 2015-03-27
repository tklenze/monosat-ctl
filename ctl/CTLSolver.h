/*
 * CTLSolver.h
 *
 *  Created on: Mar 26, 2015
 *      Author: saibot
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
#include "ctl/DynamicKripke.h"
using namespace dgl;
namespace Monosat {

class CTLSolver {
public:
	DynamicKripke k;
	CTLSolver(DynamicKripke myk) {
		k = myk;
	};
	~CTLSolver() {};


	// This represents the extra labels of states, which are sets of temporarily added atomic propositions
	// They are generated during the evaluation of a CTL formula.
	vector<vector<bool>> extralabel;


	// CTL Formula representation

	// CTL Operators, plus ID (identity operator, for formulas that are atomic) and NEG, OR, AND.
	enum CTLOp { ID, NEG, OR, AND, EX, EF, EG, EU, AX, AF, AG, AU};

	// CTL Formula tree. Everything besides "op" is only considered in some cases.
	// For instance, operand1 and operand2 are needed for AU, but value is only needed for ID
	struct CTLFormula {
		CTLOp op;
		CTLFormula *operand1;
		CTLFormula *operand2;
		int value;
	};

	Bitset* pre(Bitset& st) {
		Bitset *prest = new Bitset(k.states());
		for (int i=0; i<k.nEdgeIDs();i++) {
			if (st[k.getEdge(i).to]) {
				prest->set(k.getEdge(i).from);
			}
		}
		return prest;
	}

	Bitset* solve(CTLFormula& f) {
		switch (f.op) {
		case ID : return solveID(f);
		case NEG : return solveNEG(f);
		case OR : return solveOR(f);
		case AND : return solveAND(f);
		case EX : return solveEX(f);
		case EG : return solveEG(f);
		default : return new Bitset(k.states());
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
			prest->copyFrom(*st);
			prest = pre(*prest);

			// μ(p) ∩ pre(X).
			st->And(*prest, *andst); // andst := st ∩ prest

			if (st->Equiv(*andst))
				return andst; // fixpoint reached

			st->copyFrom(*andst);
		}

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

		//printFormula(complicated); printf("\n");

		Bitset* foo = solve(EXa);
		printStateSet(*foo);

		foo = solve(NEGEXa);
		printFormula(NEGEXa); printStateSet(*foo);

		foo = solve(NEGEXaORb);
		printFormula(NEGEXaORb); printStateSet(*foo);

		foo = solve(NEGEXaORc);
		printFormula(NEGEXaORc); printStateSet(*foo);

		foo = solve(EGb);
		printFormula(EGb); printStateSet(*foo);

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
		case EU : printf("("); printFormula(*f.operand1); printf(" EU "); printFormula(*f.operand2); printf(")"); break;

		case AX : printf("AX "); printFormula(*f.operand1); break;
		case AF : printf("AF "); printFormula(*f.operand1); break;
		case AG : printf("AG "); printFormula(*f.operand1); break;
		case AU : printf("("); printFormula(*f.operand1); printf(" AU "); printFormula(*f.operand2); printf(")"); break;

		default : printf("bar");
		}
	}
};
};

#endif /* CTL_CTLSOLVER_H_ */

