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
namespace Monosat {

template<typename B>

static CTLFormula* parseCTL(B& in);

template<typename B>

static void matchFairness(B& in, CTLFormula* f) {
	skipWhitespace(in);
	if (match(in, "{")) {
		do {
			CTLFormula* c = parseCTL(in);
			f->fairnessConstraints.push_back(c);
			skipWhitespace(in);
		} while (match(in, ","));
		if (match(in, "}")) {
		} else {
			throw std::runtime_error("Syntax Error: Curly bracket not closed!");
		}
	}
}
	template<typename B>

	// Parse the formula

	// We convert formulas with universal quantification to formulas with only existential quantification, this this
	// aids the caching in the solver
	// e.g.: AX φ ≡ ¬ EX ¬φ
	// For now, we only do this for unary predicates, because I can't be bothered and it seems to have little effect anyway
	static CTLFormula* parseCTL(B& in) {
		skipWhitespace(in);

		CTLFormula* f = newCTLFormula();


		// These are hardcoded synonyms for the Clarke/Emerson Mutex example
		if (match(in, "NCS1")) {
			f->op = ID;
			f->value = 0;
			return f;
		}
		if (match(in, "NCS2")) {
			f->op = ID;
			f->value = 3;
			return f;
		}
		if (match(in, "NCS3")) {
			f->op = ID;
			f->value = 6;
			return f;
		}
		if (match(in, "TRY1")) {
			f->op = ID;
			f->value = 1;
			return f;
		}
		if (match(in, "TRY2")) {
			f->op = ID;
			f->value = 4;
			return f;
		}
		if (match(in, "TRY3")) {
			f->op = ID;
			f->value = 7;
			return f;
		}
		if (match(in, "CS1")) {
			f->op = ID;
			f->value = 2;
			return f;
		}
		if (match(in, "CS2")) {
			f->op = ID;
			f->value = 5;
			return f;
		}
		if (match(in, "CS3")) {
			f->op = ID;
			f->value = 8;
			return f;
		}

		// Normal parsing
		if (match(in, "NEG") || match(in, "NOT") || match(in, "~") || match(in, "!")) {
			CTLFormula* inside = parseCTL(in);
			f->op = NEG;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "True") || match(in, "TRUE") || match(in, "true")) {
			f->op = True;
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
			matchFairness(in, f);
			CTLFormula* inside = parseCTL(in);
			f->op = EG;
			f->operand1 = inside;
			return f;
		}
		if (match(in, "AX")) {
			CTLFormula* phi1 = newCTLFormula();
			CTLFormula* phi2 = newCTLFormula();
			phi1->op = NEG;
			phi1->operand1 = parseCTL(in);
			phi2->op = EX;
			phi2->operand1 = phi1;
			f->op = NEG;
			f->operand1 = phi2;
			return f;
		}
		if (match(in, "AF")) {
			CTLFormula* phi1 = newCTLFormula();
			CTLFormula* phi2 = newCTLFormula();
			phi1->op = NEG;
			phi1->operand1 = parseCTL(in);
			phi2->op = EG;
			phi2->operand1 = phi1;
			f->op = NEG;
			f->operand1 = phi2;
			return f;
		}
		if (match(in, "AG")) {
			CTLFormula* phi1 = newCTLFormula();
			CTLFormula* phi2 = newCTLFormula();
			matchFairness(in, phi2);
			phi1->op = NEG;
			phi1->operand1 = parseCTL(in);
			phi2->op = EF;
			phi2->operand1 = phi1;
			f->op = NEG;
			f->operand1 = phi2;
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
			} else if (match(in, "OR") || match(in, "v") || match(in, "|")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = OR;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "AND") || match(in, "^") || match(in, "&")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AND;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket"); exit(1);
				}
			} else if (match(in, "IMP") || match(in, "->")) {
				CTLFormula* inside2 = parseCTL(in);
				skipWhitespace(in);
				if (match(in, ")")) {
					CTLFormula* notinside1 = newCTLFormula();
					notinside1->op = NEG;
					notinside1->operand1 = inside1;
					f->op = OR;
					f->operand1 = notinside1;
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
}
#endif /* CTL_PARSERUTILS_H_ */
