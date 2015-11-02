/****************************************************************************************

ATTENTION: This file is not in use, currently. It was planned to put some things from CTLParser into this static
           lump of functions, but it didn't work out. Oct 25, 2015


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