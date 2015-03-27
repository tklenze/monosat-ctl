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
	CTLSolver() {};
	~CTLSolver() {};

	DynamicKripke g;

	// This represents the extra labels of states, which are sets of temporarily added atomic propositions
	// They are generated during the evaluation of a CTL formula.
	vector<vector<bool>> extralabel;


	// CTL Formula representation

	// CTL Operators, plus ID (identity operator, for formulas that are atomic) and NEG (negation).
	enum CTLOp { EX, EF, EG, EU, AX, AF, AG, AU, ID, NEG};

	// CTL Formula tree. Everything besides "op" is only considered in some cases.
	// For instance, operand1 and operand2 are needed for AU, but value is only needed for ID
	struct CTLFormula {
		CTLOp op;
		CTLFormula *operand1;
		CTLFormula *operand2;
		int value;
	};


	void funwithctl() {
		CTLFormula a {ID, NULL, NULL, 3};
		CTLFormula b {ID, NULL, NULL, 7};
		CTLFormula c {ID, NULL, NULL, 8};
		CTLFormula EXa {EX, &a, NULL, 1};
		CTLFormula EUEXab {EU, &EXa, &b, 0};
		CTLFormula complicated {AF, &EUEXab, NULL, 0};

		printFormula(complicated);
	}
	void printFormula(CTLFormula f) {
		switch (f.op) {
		case EX : printf("EX "); printFormula(*f.operand1); break;
		case EF : printf("EF "); printFormula(*f.operand1); break;
		case EG : printf("EG "); printFormula(*f.operand1); break;
		case EU : printf("("); printFormula(*f.operand1); printf(" EU "); printFormula(*f.operand2); printf(")"); break;

		case AX : printf("AX "); printFormula(*f.operand1); break;
		case AF : printf("AF "); printFormula(*f.operand1); break;
		case AG : printf("AG "); printFormula(*f.operand1); break;
		case AU : printf("("); printFormula(*f.operand1); printf(" AU "); printFormula(*f.operand2); printf(")"); break;

		case ID : printf("%d", f.value); break;
		case NEG: printf("NOT "); printFormula(*f.operand1); break;
		default : printf("bar");
		}
	}
};
};

#endif /* CTL_CTLSOLVER_H_ */

