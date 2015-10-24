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

#ifndef CTL_PARSERUTILS_H_
#define CTL_PARSERUTILS_H_

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "ctl/CTLTheory.h"
#include "ctl/CTLFormula.h"
#include "ctl/CTLParser.h"

#include "core/Config.h"
#include "pb/PbTheory.h"
#include "core/Dimacs.h"
#include <gmpxx.h>
#include <set>
#include <string>
#include <sstream>

	static CTLFormula* parseCTL(char* in) {
		skipWhitespace(in);
		CTLFormula* f = newCTLFormula();
		if (match(in, "NEG")) {
			CTLFormula* inside = parseCTL(in);
			f->op = NEG;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "NOT")) { // out of convenience, identical to NEG
			CTLFormula* inside = parseCTL(in);
			f->op = NEG;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "EX")) {
			CTLFormula* inside = parseCTL(in);
			f->op = EX;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "EF")) {
			CTLFormula* inside = parseCTL(in);
			f->op = EF;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "EG")) {
			CTLFormula* inside = parseCTL(in);
			f->op = EG;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "AX")) {
			CTLFormula* inside = parseCTL(in);
			f->op = AX;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "AF")) {
			CTLFormula* inside = parseCTL(in);
			f->op = AF;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "AG")) {
			CTLFormula* inside = parseCTL(in);
			f->op = AG;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "(")) {
			CTLFormula* inside1 = parseCTL(in);
			skipWhitespace(in);
			if (match(in, "EW")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = EW;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "EU")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = EU;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "AW")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AW;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "AU")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AU;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "OR")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = OR;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "AND")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AND;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else {
				printf("Error: Operator missing, remainder is "); std::cout << std::string(in); exit(1);
			}
			return f;
		}
		// No other matches, must be a number.
		int n = parseInt(in);
		f->op = ID;
		f->value = n;

		return f;
	}

#endif /* CTL_PARSERUTILS_H_ */
