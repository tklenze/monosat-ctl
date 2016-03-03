/*
 * CTLFormula.h
 *
 * Written by Tobias Klenze and Sam Bayless
 *
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
		bool isBinaryPredicate() {
			if (this->op == ID || this->op == True || this->op == NEG || this->op == EX || this->op == EF || this->op == EG || this->op == AX || this->op == AF || this->op == AG)
				return false;
			else
				return true;
		}
		bool isPropositional(bool allow_nestedX=false){
			if(allow_nestedX){
				return isPropositionalWithNoNestedX();
			}
			if (this->op == EX || this->op == EF || this->op == EG || this->op == EU || this->op == EW || this->op == AX || this->op == AF || this->op == AG || this->op == AU || this->op == AW)
				return false;
			if (this->op == ID || this->op == True)
				return true;
			if (this->isBinaryPredicate())
				return this->operand1->isPropositional() && this->operand2->isPropositional();
			else
				return this->operand1->isPropositional();
		}
		//True if this formula is propositional, or has no nested EX/AX constraints
		bool isPropositionalWithNoNestedX(){
			return isPropositionalWithNoNestedX_helper(0);
		}
		bool isPropositionalWithNoNestedX_helper(int x_count){
			if (this->op == EF || this->op == EG || this->op == EU || this->op == EW || this->op == AX || this->op == AF || this->op == AG || this->op == AU || this->op == AW)
				return false;
			if (this->op == ID || this->op == True)
				return true;
			else if (this->op==EX){
				if(x_count>0)
					return false;//found nexted ex/ax
				else{
					return this->operand1->isPropositionalWithNoNestedX_helper(x_count+1);
				}
			}else if (this->op==AX && x_count>0){
				if(x_count>0)
					return false;//found nexted ex/ax
				else{
					return this->operand1->isPropositionalWithNoNestedX_helper(x_count+1);
				}
			}else if (this->isBinaryPredicate())
				return this->operand1->isPropositionalWithNoNestedX_helper(x_count) && this->operand2->isPropositionalWithNoNestedX_helper(x_count);
			else
				return this->operand1->isPropositionalWithNoNestedX_helper(x_count);
		}
		void addFairnessConstraint(CTLFormula*f){
			fairnessConstraints.push_back(f);
		}
		std::vector<CTLFormula*> & getFairnessConstraints(){
			return fairnessConstraints;
		}

		//Create a copy of this formula
		CTLFormula * copy(){
			CTLFormula * f = new CTLFormula(this->op,nullptr,nullptr,this->value);

			f->id = this->id;

			if(this->operand1){
				f->operand1=this->operand1->copy();
			}

			if(this->operand2){
				f->operand2=this->operand2->copy();
			}

			for(CTLFormula * fairness:this->fairnessConstraints){
				f->fairnessConstraints.push_back(fairness->copy());
			}
			return f;
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
