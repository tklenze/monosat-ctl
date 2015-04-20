/*
 * CTLFormula.h
 *
 *  Created on: Apr 13, 2015
 *      Author: tobias k
 *
 */

#ifndef CTL_FORMULA_H_
#define CTL_FORMULA_H_
#include <stdio.h>

	// CTL Formula representation

	// CTL Operators, plus ID (identity operator, for formulas that are atomic) and NEG, OR, AND.
	enum CTLOp { ID, NEG, OR, AND, EX, EF, EG, EW, EU, AX, AF, AG, AW, AU};

	/* CTL Formula tree. Everything besides "op" is only considered in some cases.
	 * For instance, operand1 and operand2 are needed for AU, but value is only meaningful in case op=ID
	 *
	 * Sample formulas:
	 *   {ID, NULL, NULL, 3} -- corresponds to the atomic proposition number 3. Note how operands are ignored.
	 *   {EF, phi, NULL, 0} -- corresponds to "EF phi". Note that "value" and operand2 are ignored.
	 */
	struct CTLFormula {
		CTLOp op;
		CTLFormula *operand1;
		CTLFormula *operand2;
		int value;
	};

	CTLFormula* newCTLFormula() {
		CTLFormula f {ID, NULL, NULL, 0};
		return &f;
	}

	void printFormula(CTLFormula* foo) {
		CTLFormula f = *foo;
		switch (f.op) {
		case ID : printf("%d", f.value); break;
		case NEG: printf("not "); printFormula(f.operand1); break;
		case OR : printf("("); printFormula(f.operand1); printf(" or "); printFormula(f.operand2); printf(")"); break;
		case AND : printf("("); printFormula(f.operand1); printf(" and "); printFormula(f.operand2); printf(")"); break;

		case EX : printf("EX "); printFormula(f.operand1); break;
		case EF : printf("EF "); printFormula(f.operand1); break;
		case EG : printf("EG "); printFormula(f.operand1); break;
		case EW : printf("("); printFormula(f.operand1); printf(" EW "); printFormula(f.operand2); printf(")"); break;
		case EU : printf("("); printFormula(f.operand1); printf(" EU "); printFormula(f.operand2); printf(")"); break;

		case AX : printf("AX "); printFormula(f.operand1); break;
		case AF : printf("AF "); printFormula(f.operand1); break;
		case AG : printf("AG "); printFormula(f.operand1); break;
		case AW : printf("("); printFormula(f.operand1); printf(" AW "); printFormula(f.operand2); printf(")"); break;
		case AU : printf("("); printFormula(f.operand1); printf(" AU "); printFormula(f.operand2); printf(")"); break;

		default : printf("Unknown formula: %d", f.op);
		}
	}


#endif /* CTL_FORMULA_H_ */
