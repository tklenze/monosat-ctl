/*
 * Written by Tobias Klenze and Sam Bayless
 */

#include "CTLFormula.h"
namespace Monosat{
CTLFormula* newCTLFormula() {
		//std::vector<CTLFormula *> fairnessC;
		CTLFormula * f = new CTLFormula(ID, nullptr, nullptr);

		return f;
	}



	void printFairnessConstraints(CTLFormula* f) {
		if (f->fairnessConstraints.size() > 0) {
			printf("{");
			printFormula(f->fairnessConstraints.operator [](0));
			for (int i = 1; i < f->fairnessConstraints.size(); i++) {
				printf(", ");
				printFormula(f->fairnessConstraints.operator [](i));
			}
			printf("} ");
		}
	}

	void printFormula(CTLFormula* foo) {
		CTLFormula f = *foo;
		switch (f.op) {
		case ID : printf("%d", f.value); break;
		case True : printf("True"); break;
		case NEG: printf("not "); printFormula(f.operand1); break;
		case OR : printf("("); printFormula(f.operand1); printf(" or "); printFormula(f.operand2); printf(")"); break;
		case AND : printf("("); printFormula(f.operand1); printf(" and "); printFormula(f.operand2); printf(")"); break;

		case EX : printf("EX "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case EF : printf("EF "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case EG : printf("EG "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case EW : printf("("); printFormula(f.operand1); printf(" EW "); printFairnessConstraints(foo); printFormula(f.operand2); printf(")"); break;
		case EU : printf("("); printFormula(f.operand1); printf(" EU "); printFairnessConstraints(foo); printFormula(f.operand2); printf(")"); break;

		case AX : printf("AX "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case AF : printf("AF "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case AG : printf("AG "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case AW : printf("("); printFormula(f.operand1); printf(" AW "); printFairnessConstraints(foo); printFormula(f.operand2); printf(")"); break;
		case AU : printf("("); printFormula(f.operand1); printf(" AU "); printFairnessConstraints(foo); printFormula(f.operand2); printf(")"); break;

		default : printf("Unknown formula: %d", f.op);
		}
	}

	// FIXME tidy up this mess of printFormula and printFormula2
	void printFormula2(CTLFormula& f) {
		//CTLFormula f = *foo;
		switch (f.op) {
		case ID : printf("%d", f.value); break;
		case True : printf("True"); break;
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

	std::string getFormulaNuSMVFormat(CTLFormula& f) {

		// NuSMV doesn't understand weak until, that's why we convert the formula to use other operators
		CTLFormula * phiEW1 = new CTLFormula(EG, f.operand1, nullptr);
		CTLFormula * phiEW2 = new CTLFormula(EU, f.operand1, f.operand2);
		CTLFormula * phiEW3 = new CTLFormula(OR, phiEW1, phiEW2);

		CTLFormula * phiAW1 = new CTLFormula(OR, f.operand1, f.operand2);
		CTLFormula * phiAW2 = new CTLFormula(NEG, phiAW1, nullptr);
		CTLFormula * phiAW3 = new CTLFormula(NEG, f.operand2, nullptr);
		CTLFormula * phiAW4 = new CTLFormula(EU, phiAW3, phiAW2);
		CTLFormula * phiAW5 = new CTLFormula(NEG, phiAW4, nullptr);

		std::string s;
		switch (f.op) {
		case ID : s = "v" + std::to_string(f.value); break;
		case True : s = "TRUE"; break;
		case NEG: s = "!" + getFormulaNuSMVFormat(*f.operand1); break;
		case OR : s = "(" + getFormulaNuSMVFormat(*f.operand1) + " | " + getFormulaNuSMVFormat(*f.operand2) + ")"; break;
		case AND : s = "(" + getFormulaNuSMVFormat(*f.operand1) + " & " + getFormulaNuSMVFormat(*f.operand2) + ")"; break;

		case EX: s = "EX " + getFormulaNuSMVFormat(*f.operand1); break;
		case EF: s = "EF " + getFormulaNuSMVFormat(*f.operand1); break;
		case EG: s = "EG " + getFormulaNuSMVFormat(*f.operand1); break;
		case EW: s = getFormulaNuSMVFormat(*phiEW3); break;
		case EU: s = "E [" + getFormulaNuSMVFormat(*f.operand1) + " U " + getFormulaNuSMVFormat(*f.operand2) + " ]"; break;

		case AX: s = "AX " + getFormulaNuSMVFormat(*f.operand1); break;
		case AF: s = "AF " + getFormulaNuSMVFormat(*f.operand1); break;
		case AG: s = "AG " + getFormulaNuSMVFormat(*f.operand1); break;
		case AW: s = getFormulaNuSMVFormat(*phiAW5); break;
		case AU : s = "A [" + getFormulaNuSMVFormat(*f.operand1) + " U " + getFormulaNuSMVFormat(*f.operand2) + " ]"; break;

		default : s = "Unknown formula";
		}
		return s;
	}

	std::string getFormulaPctlFormat(CTLFormula& f) {

		// Pctl doesn't understand most of the operators we are using, that's why we convert the formula to use other operators
		CTLFormula * phiEW1 = new CTLFormula(EG, f.operand1, nullptr);
		CTLFormula * phiEW2 = new CTLFormula(EU, f.operand1, f.operand2);
		CTLFormula * phiEW3 = new CTLFormula(OR, phiEW1, phiEW2);

		CTLFormula * phiAW1 = new CTLFormula(OR, f.operand1, f.operand2);
		CTLFormula * phiAW2 = new CTLFormula(NEG, phiAW1, nullptr);
		CTLFormula * phiAW3 = new CTLFormula(NEG, f.operand2, nullptr);
		CTLFormula * phiAW4 = new CTLFormula(EU, phiAW3, phiAW2);
		CTLFormula * phiAW5 = new CTLFormula(NEG, phiAW4, nullptr);

		CTLFormula * phiEG1 = new CTLFormula(True, nullptr, nullptr);
		CTLFormula * phiEG11 = new CTLFormula(NEG, f.operand1, nullptr);
		CTLFormula * phiEG2 = new CTLFormula(AU, phiEG1, phiEG11);
		CTLFormula * phiEG3 = new CTLFormula(NEG, phiEG2, nullptr);

		CTLFormula * phiAF1 = new CTLFormula(NEG, f.operand1, nullptr);
		CTLFormula * phiAF2 = new CTLFormula(EG, phiAF1, nullptr);
		CTLFormula * phiAF3 = new CTLFormula(NEG, phiAF2, nullptr);

		CTLFormula * phiAG1 = new CTLFormula(NEG, f.operand1, nullptr);
		CTLFormula * phiAG2 = new CTLFormula(EF, phiAG1, nullptr);
		CTLFormula * phiAG3 = new CTLFormula(NEG, phiAG2, nullptr);

		// undo double negation
		if (f.op == NEG && f.operand1->op == NEG)
			return getFormulaPctlFormat(*f.operand1->operand1);

		std::string s;
		switch (f.op) {
		case ID : s = "v" + std::to_string(f.value); break;
		case True : s = "true"; break;
		case NEG: s = "not " + getFormulaPctlFormat(*f.operand1); break;
		case OR : s = "(" + getFormulaPctlFormat(*f.operand1) + " or " + getFormulaPctlFormat(*f.operand2) + ")"; break;
		case AND : s = "(" + getFormulaPctlFormat(*f.operand1) + " and " + getFormulaPctlFormat(*f.operand2) + ")"; break;

		case EX: s = "P>0/1 (X " + getFormulaPctlFormat(*f.operand1) + ")"; break;
		case EF: s = "P>0/1 (true U " + getFormulaPctlFormat(*f.operand1) + ")"; break;
		case EG: s = getFormulaPctlFormat(*phiEG3); break;
		case EW: s = getFormulaPctlFormat(*phiEW3); break;
		case EU: s = "P>0/1 (" + getFormulaPctlFormat(*f.operand1) + " U " + getFormulaPctlFormat(*f.operand2) + " )"; break;

		case AX: s = "P>=1/1 (X " + getFormulaPctlFormat(*f.operand1) + ")"; break;
		case AF: s = getFormulaPctlFormat(*phiAF3); break;
		case AG: s = getFormulaPctlFormat(*phiAG3); break;
		case AW: s = getFormulaPctlFormat(*phiAW5); break;
		case AU: s = "P>=1/1 (" + getFormulaPctlFormat(*f.operand1) + " U " + getFormulaPctlFormat(*f.operand2) + " )"; break;

		default : s = "getFormulaPctlFormat: Unknown formula";
		}
		return s;
	}

	std::string getFormulaCTLSATFormat(CTLFormula* foo) {
		CTLFormula f = *foo;

		// Pctl doesn't understand most of the operators we are using, that's why we convert the formula to use other operators
		CTLFormula * phiEW1 = new CTLFormula(EG, f.operand1, nullptr);
		CTLFormula * phiEW2 = new CTLFormula(EU, f.operand1, f.operand2);
		CTLFormula * phiEW3 = new CTLFormula(OR, phiEW1, phiEW2);

		CTLFormula * phiAW1 = new CTLFormula(OR, f.operand1, f.operand2);
		CTLFormula * phiAW2 = new CTLFormula(NEG, phiAW1, nullptr);
		CTLFormula * phiAW3 = new CTLFormula(NEG, f.operand2, nullptr);
		CTLFormula * phiAW4 = new CTLFormula(EU, phiAW3, phiAW2);
		CTLFormula * phiAW5 = new CTLFormula(NEG, phiAW4, nullptr);

		std::string s;
		switch (f.op) {
		case ID : switch (f.value) {
		case  0 : s = "p"; break;
		case  1 : s = "q"; break;
		case  2 : s = "r"; break;
		case  3 : s = "s"; break;
		case  4 : s = "u"; break;
		case  5 : s = "w"; break;
		case  6 : s = "a"; break;
		case  7 : s = "b"; break;
		case  8 : s = "c"; break;
		case  9 : s = "d"; break;
		case 10 : s = "e"; break;
		case 11 : s = "f"; break;
		case 12 : s = "g"; break;
		case 13 : s = "h"; break;
		case 14 : s = "i"; break;
		case 15 : s = "j"; break;
		case 16 : s = "k"; break;
		case 17 : s = "l"; break;
		case 18 : s = "m"; break;
		case 19 : s = "n"; break;
		case 20 : s = "o"; break;
		default : s = "ERROR: too many atomic propositions in formula";
		} break;
		case True :  s+="T"; break;
		case NEG: s = "~ " + getFormulaCTLSATFormat(f.operand1); break;
		case OR : s = "(" + getFormulaCTLSATFormat(f.operand1) + " v " + getFormulaCTLSATFormat(f.operand2) + ")"; break;
		case AND : s = "(" + getFormulaCTLSATFormat(f.operand1) + " ^ " + getFormulaCTLSATFormat(f.operand2) + ")"; break;

		case EX: s = "EX " + getFormulaCTLSATFormat(f.operand1); break;
		case EF: s = "EF " + getFormulaCTLSATFormat(f.operand1); break;
		case EG: s = "EG " + getFormulaCTLSATFormat(f.operand1); break;
		case EW: s = getFormulaCTLSATFormat(phiEW3); break;
		case EU: s = "E(" + getFormulaCTLSATFormat(f.operand1) + " U " + getFormulaCTLSATFormat(f.operand2) + " )"; break;

		case AX: s = "AX " + getFormulaCTLSATFormat(f.operand1); break;
		case AF: s = "AF " + getFormulaCTLSATFormat(f.operand1); break;
		case AG: s = "AG " + getFormulaCTLSATFormat(f.operand1); break;
		case AW: s = getFormulaCTLSATFormat(phiAW5); break;
		case AU: s = "A(" + getFormulaCTLSATFormat(f.operand1) + " U " + getFormulaCTLSATFormat(f.operand2) + " )"; break;

		default : printf("getFormulaCTLSATFormat: Unknown formula: %d", f.op);
		}
		return s;
	}
};
