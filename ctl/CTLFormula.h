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
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <vector>

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

	CTLFormula* newCTLFormula() ;

	void printFormula(CTLFormula* foo);

	void printFairnessConstraints(CTLFormula* f) ;



	// FIXME tidy up this mess of printFormula and printFormula2
	void printFormula2(CTLFormula& f) ;

	std::string getFormulaNuSMVFormat(CTLFormula& f) ;

	std::string getFormulaPctlFormat(CTLFormula& f) ;

	std::string getFormulaCTLSATFormat(CTLFormula* foo);


}
#endif /* CTL_FORMULA_H_ */
