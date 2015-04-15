/*
 * CTLFormula.h
 *
 *  Created on: Apr 13, 2015
 *      Author: tobias k
 *
 */

#ifndef CTL_FORMULA_H_
#define CTL_FORMULA_H_

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



#endif /* CTL_FORMULA_H_ */
