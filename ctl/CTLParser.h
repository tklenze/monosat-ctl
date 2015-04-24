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

#ifndef CTL_PARSER_H_
#define CTL_PARSER_H_

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "ctl/CTLTheory.h"
#include "ctl/CTLFormula.h"

#include "core/Config.h"
#include "pb/PbTheory.h"
#include "core/Dimacs.h"
#include <gmpxx.h>
#include <set>
#include <string>
#include <sstream>
namespace Monosat {


//=================================================================================================
// Kripke Parser:
template<class B, class Solver>
class CTLParser: public Parser<B, Solver> {

	vec<CTLTheorySolver*> kripkes;

	vec<vec<SteinerStruct*>> steiners;
	PbTheory * pbtheory = nullptr;

	vec<Lit> lits;
	int edgeCount = 0;
	int nodeCount = 0;
	vec<char> tmp;

	void readKripke(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		int a, g, n, e, ev;

		n = parseInt(in); //num nodes
		e = parseInt(in); //num edges (I'm ignoring this currently)
		//  ev = parseInt(in);//the variable of the first graph edge.
		a = parseInt(in);  //number of APs
		g = parseInt(in);  //id of the graph

		kripkes.growTo(g + 1);

		CTLTheorySolver *kripke = new CTLTheorySolver(&S, g);
		kripke->newNodes(n);
		kripkes[g] = kripke;
		S.addTheory(kripke);

		printf("Kripkes.size: %d\n", kripkes.size());

		//  return ev;
	}

	void readNodeAP(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int kripkeID = parseInt(in);
		int node = parseInt(in);
		int ap = parseInt(in);
		int nodeVar = parseInt(in) - 1;
			// TODO
			// should set it to undecided
		kripkes[kripkeID]->newNodeAP(node, ap, nodeVar);

	}

	
	void readEdge(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		++in;

		int kripkeID = parseInt(in);
		int from = parseInt(in);
		int to = parseInt(in);
		int edgeVar = parseInt(in) - 1;

		if (kripkeID < 0 || kripkeID >= kripkes.size()) {
			printf("DEBUG: kripkeID %d\n", kripkeID);
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", kripkeID, edgeVar), exit(1);
		}
		if (edgeVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar), exit(1);
		}
		while (edgeVar >= S.nVars())
			S.newVar();

		skipWhitespace(in);
		if(*in=='\n' || *in==0){

			if (kripkes[kripkeID]) {
				// TODO not implemented yet
				kripkes[kripkeID]->newTransition(from, to, edgeVar);
			} else {
				printf("PARSE ERROR! Undeclared kripke identifier %d for edge %d\n", kripkeID, edgeVar), exit(1);
				exit(1);
			}

		}

	}

	void readCTL(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		++in;

		int kripkeID = parseInt(in);
		int initialNode = parseInt(in);
		int ctlVar = parseInt(in) - 1;

		CTLFormula* f = parseCTL(in);

		printf("Parsed formula:\n");
		printFormula(f);
		printf("\n");

		kripkes[kripkeID]->setCTL(*f, initialNode);
		kripkes[kripkeID]->newCTLVar(ctlVar);

		return;
	}

	// 	enum CTLOp { ID, NEG, OR, AND, EX, EF, EG, EW, EU, AX, AF, AG, AW, AU};

	CTLFormula* parseCTL(B& in) {
		skipWhitespace(in);
		CTLFormula* f = newCTLFormula();
		if (match(in, "NEG")) {
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
			printf("DEBUG: Parsed first part, "); printFormula(inside1); printf("\n");
			skipWhitespace(in);
			if (match(in, "EW")) {
				CTLFormula* inside2 = parseCTL(in);
				printf("DEBUG: Parsed second part, "); printFormula(inside2); printf("\n");
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = EW;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket");
				}
			} else if (match(in, "EU")) {
				CTLFormula* inside2 = parseCTL(in);
				printf("DEBUG: Parsed second part, "); printFormula(inside2); printf("\n");
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = EU;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket");
				}
			} else if (match(in, "AW")) {
				CTLFormula* inside2 = parseCTL(in);
				printf("DEBUG: Parsed second part, "); printFormula(inside2); printf("\n");
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AW;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket");
				}
			} else if (match(in, "AU")) {
				CTLFormula* inside2 = parseCTL(in);
				printf("DEBUG: Parsed second part, "); printFormula(inside2); printf("\n");
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AU;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket");
				}
			} else if (match(in, "OR")) {
				CTLFormula* inside2 = parseCTL(in);
				printf("DEBUG: Parsed second part, "); printFormula(inside2); printf("\n");
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = OR;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket");
				}
			} else if (match(in, "AND")) {
				CTLFormula* inside2 = parseCTL(in);
				printf("DEBUG: Parsed second part, "); printFormula(inside2); printf("\n");
				skipWhitespace(in);
				if (match(in, ")")) {
					f->op = AND;
					f->operand1 = inside1;
					f->operand2 = inside2;
				} else {
					printf("Error: Expected closing bracket");
				}
			} else {
				printf("Error: Operator missing");
			}
			printf("DEBUG: Parsed bracket formula, "); printFormula(f); printf("\n");
			return f;
		}
		// ... TODO, other matches

		// No other matches, must be a number.
		int n = parseInt(in);
		f->op = ID;
		f->value = n;
		return f;
	}

	// Draws all Kripke structures for debugging purposes
	void readDraw(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		printf("Printing Kripke structures\n");
		for (int i = 0; i<kripkes.size(); i++) {
			printf("Kripke structure, ID: %d, under\n", i);
			kripkes[i]->g_under->draw();
			printf("Kripke structure, ID: %d, over\n", i);
			kripkes[i]->g_over->draw();

			printf("\n\n");
		}

		enum VarType { EDGE, NODEAP, DETECTOR}; // FIXME needs to be synchronized with CTLTheory.VarType
		//enum VarType = kripkes[0].VarType; // does not work

		for (int i = 0; i<kripkes.size(); i++) {
				printf("Kripke structure %d, Vars:\n", i);
			for (int j = 0; j<kripkes[i]->vars.size(); j++) {
				switch (kripkes[i]->vars[j].type) {
				case EDGE :  printf(" Edge, EdgeID: %d, solverVar: %d\n", kripkes[i]->vars[j].detector_node_edge, kripkes[i]->vars[j].solverVar); break;
				case NODEAP :  printf(" NodeAP, NodeID: %d, solverVar: %d\n", kripkes[i]->vars[j].detector_node_edge, kripkes[i]->vars[j].solverVar); break;
				case DETECTOR :  printf(" Detector, DetectorID: %d, solverVar: %d\n", kripkes[i]->vars[j].detector_node_edge, kripkes[i]->vars[j].solverVar); break;
				}
			}
			//kripkes[i]->S->all_theory_vars
			for (int j = 0; j < kripkes[i]->S->all_theory_vars.size(); j++) {
				printf("%d ", kripkes[i]->S->all_theory_vars[j]);
			}
			printf("\n");

		}

		return;
	}

public:
	bool parseLine(B& in, Solver& S) {
		/* DEBUG */
		printf("Parse line: ");

		for (int i = 0; in[i] != '\0'; i++) {
			printf("%c", in[i]);
		}
		printf("\n");
		// */

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (*in == 'c') {
			//just a comment
			return false;
		} else if (match(in, "kripke")) {
			skipWhitespace(in);
			readKripke(in, S);
			skipWhitespace(in);
			return true;
		} else if (match(in, "kedge")) {
			edgeCount++;
			readEdge(in, S);
			return true;
		} else if (match(in, "knodeap")) {
			nodeCount++;
			readNodeAP(in, S);
			return true;
		} else if (match(in, "kctl")) {
			readCTL(in, S);
			return true;
		} else if (match(in, "kdraw")) {
			readDraw(in, S);
			return true;
		}
		return false;
	}
	
	void implementConstraints(Solver & S) {

		//not really implemented, yet!
	/*	for (int gid = 0; gid < steiners.size(); gid++) {
			for (auto & steiner : steiners[gid]) {
				if (steiner) {
					graphs[gid]->steinerTree(steiner->terminals, steiner->id);
					for (auto & weight : steiner->weight_constraints) {
						graphs[gid]->addSteinerWeightConstraint(steiner->id, weight.first, weight.second);
					}
					delete (steiner);
				}
			}
		}*/

		for (int i = 0; i < kripkes.size(); i++) {
			/*
			if (kripkes[i]) {
				// TODO not implemented
				// kripkes[i]->implementConstraints();
			}
			*/
		}
		if (pbtheory)
			pbtheory->implementConstraints();
	}

	
};

//=================================================================================================
}
;

#endif /* CTL_PARSER_H_ */
