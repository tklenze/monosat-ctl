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
#include "ctl/CTLParser.h"
#include <initializer_list>

using namespace dgl;
namespace Monosat {

// CTL Formula representation

	// CTL Operators, plus ID (identity operator, for formulas that are atomic) and NEG, OR, AND.
	enum CTLOp { ID, True, NEG, OR, AND, EX, EF, EG, EW, EU, AX, AF, AG, AW, AU};

	/* CTL Formula tree. Everything besides "op" is only considered in some cases.
	 * For instance, operand1 and operand2 are needed for AU, but value is only meaningful in case op=ID
	 *
	 * Sample formulas:
	 *   {ID, nullptr, nullptr, 3} -- corresponds to the atomic proposition number 3. Note how operands are ignored.
	 *   {EF, phi, nullptr, 0} -- corresponds to "EF phi". Note that "value" and operand2 are ignored.
	 */
	struct CTLFormula {
		CTLOp op;
		CTLFormula *operand1;
		CTLFormula *operand2;
		int value=0;

		//use for caching.
		//Optional, must be a unique non-negative integer or -1
		//(but if -1, caching will be disabled for this formula)
		int id = -1;
		std::vector<CTLFormula *> fairnessConstraints;

		CTLFormula(CTLOp op, CTLFormula * op1, CTLFormula * op2=nullptr,int value =0):op(op),operand1(op1),operand2(op2),value(value){

		}
		CTLFormula(CTLOp op, CTLFormula * op1, CTLFormula * op2, int value, std::vector<CTLFormula *> fairnessConstraints):op(op),operand1(op1),operand2(op2),value(value),fairnessConstraints(fairnessConstraints){

		}
		void setValue(int value){
			this->value=value;
		}
		void setID(int id){
			this->id = id;
		}
		int getID(){
			return id;
		}
		void addFairnessConstraint(CTLFormula*f){
			fairnessConstraints.push_back(f);
		}
		std::vector<CTLFormula*> & getFairnessConstraints(){
			return fairnessConstraints;
		}
	};

	CTLFormula* newCTLFormula() {
		//std::vector<CTLFormula *> fairnessC;
		CTLFormula * f = new CTLFormula(ID, nullptr, nullptr);

		return f;
	}

	void printFormula(CTLFormula* foo);

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

		case EX : printf("EX "); printFormula(f.operand1); break;
		case EF : printf("EF "); printFormula(f.operand1); break;
		case EG : printf("EG "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case EW : printf("("); printFormula(f.operand1); printf(" EW "); printFormula(f.operand2); printf(")"); break;
		case EU : printf("("); printFormula(f.operand1); printf(" EU "); printFormula(f.operand2); printf(")"); break;

		case AX : printf("AX "); printFormula(f.operand1); break;
		case AF : printf("AF "); printFormula(f.operand1); break;
		case AG : printf("AG "); printFairnessConstraints(foo); printFormula(f.operand1); break;
		case AW : printf("("); printFormula(f.operand1); printf(" AW "); printFormula(f.operand2); printf(")"); break;
		case AU : printf("("); printFormula(f.operand1); printf(" AU "); printFormula(f.operand2); printf(")"); break;

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

		default : s = "Unknown formula";
		}
		return s;
	}
}
#endif /* CTL_FORMULA_H_ */
