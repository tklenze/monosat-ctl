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
#include "core/Dimacs.h"
#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "ctl/CTLTheory.h"
#include "ctl/CTLFormula.h"
#include "ctl/CTLParserUtils.h"

#include "core/Config.h"
#include "pb/PbTheory.h"
#include "core/Dimacs.h"
#include <gmpxx.h>
#include <set>
#include <string>
#include <sstream>
#include <unistd.h>
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
	int currentKripkeID; // when reading multiple lines, we have to keep track of which Kripke structure we are talking about.
	int currentInitialNode;
	CTLFormula* f; // when reading multiple lines, keep track of the already passed subformula


	// kripke <#Nodes> <#Edges> <#APs> <KripkeID>
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
		kripke->g_under->setAPCount(a);
		kripke->g_over->setAPCount(a);
		nodeCount = n;
		kripke->initNodeAPVarLookup(n, a);
		kripkes[g] = kripke;
		S.addTheory(kripke);


		//  return ev;
	}

	// Shorthand form for a simpler input format, where edges and nodeAPs don't have to be explicitly added
	// kctlsimp <KripkeID> <#Nodes> <#APs> <selfloops> <ctlvar> <CTL formula>"
	void readCTLSimp(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}


		int kripkeID = parseInt(in);
		int n = parseInt(in); //num node
		int a = parseInt(in);  //number of APs
		int loops = parseInt(in);  //loops: 0 or 1
		int ctlVar = parseInt(in) - 1;


		kripkes.growTo(kripkeID + 1);

		CTLTheorySolver *kripke = new CTLTheorySolver(&S, kripkeID);
		kripke->newNodes(n);
		kripke->g_under->setAPCount(a);
		kripke->g_over->setAPCount(a);
		nodeCount = n;
		kripke->initNodeAPVarLookup(n, a);
		kripkes[kripkeID] = kripke;
		currentKripkeID = kripkeID;
		currentInitialNode = 0;

		while (n * (n - 1 + loops) + n * a + 2 >= S.nVars())
			S.newVar();

		int var = 0;
		// Add edges
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < n; j ++) {
				if (i!=j || loops == 1) {
					kripkes[kripkeID]->newTransition(i, j, var);
					var++;
				}
			}
		}
		// Add stateAPs
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < a; j ++) {
				kripkes[kripkeID]->newNodeAP(i, j, var);
				var++;
			}
		}
		// if we wanted, we could replace the explicit ctlVar by: (n * (n - 1 + loops) + n * a)
		addCTL(in, S, ctlVar); // set initial node to 0

		S.addTheory(kripke);
	}

	// Even simpler format, which assumes that there are a certain number of processes and each process is in exactly one of a certain number of states
	// The APs will be numbered continuously, e.g. for two processes and three states per process:
	// AP 0: process 0, state 0
	// AP 5: process 1, state 2
	// The advantage of using kctlsinglestate is that certain optimizations can be done outside of the CTL formula
	//c kctlsinglestate <KripkeID> <#Nodes> <#Processes> <States/Process> <selfloops> <ctlvar> <CTL formula>"
	void readCTLSingleState(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}

		int kripkeID = parseInt(in);
		int n = parseInt(in); //num node
		int p = parseInt(in);  //number of Processes
		int s = parseInt(in);  //number of States per Process
		int a = p * s;
		int loops = parseInt(in);  //loops: 0 or 1
		int ctlVar = parseInt(in) - 1;

		kripkes.growTo(kripkeID + 1);

		CTLTheorySolver *kripke = new CTLTheorySolver(&S, kripkeID);
		kripke->newNodes(n);
		kripke->g_under->setAPCount(a);
		kripke->g_over->setAPCount(a);
		kripke->processes = p;
		kripke->statesperprocess = s;
		nodeCount = n;
		kripke->initNodeAPVarLookup(n, a);
		kripkes[kripkeID] = kripke;
		currentKripkeID = kripkeID;
		currentInitialNode = 0;

		while (n * (n - 1 + loops) + n * a + 2 >= S.nVars())
			S.newVar();

		int var = 0;
		// Add edges
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < n; j ++) {
				if (i!=j || loops == 1) {
					kripkes[kripkeID]->newTransition(i, j, var);
					var++;
				}
			}
		}
		// Add stateAPs
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < a; j ++) {
				kripkes[kripkeID]->newNodeAP(i, j, var);
				var++;
			}
		}
		// if we wanted, we could replace the explicit ctlVar by: (n * (n - 1 + loops) + n * a)
		addCTL(in, S, ctlVar); // set initial node to 0

		S.addTheory(kripke);
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
			parse_errorf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", kripkeID, edgeVar);
		}
		if (edgeVar < 0) {
			parse_errorf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar);
		}
		while (edgeVar >= S.nVars())
			S.newVar();

		skipWhitespace(in);
		if(*in=='\n' || *in==0){

			if (kripkes[kripkeID]) {
				// TODO not implemented yet
				kripkes[kripkeID]->newTransition(from, to, edgeVar);
			} else {
				parse_errorf("PARSE ERROR! Undeclared kripke identifier %d for edge %d\n", kripkeID, edgeVar);

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
		currentKripkeID = kripkeID;
		currentInitialNode = initialNode;
		addCTL(in, S, ctlVar);
	}

	// Helper function that assumes everything besides the formula has already been read.
	void addCTL(B& in, Solver& S, int ctlVar) {

		// OK, this is probably a bit controversial, but I will attach AND the formula together with AG (EX True).
		//   properly document this, and make it able to switch it off
		// done in CTLTheory's preprocessing instead now!

		////Sam: Don't do it this way.
		std::string tmp = std::string("(AG EX True AND ")+in+std::string(")");
		char *fSafe = &tmp[0];
		f = parseCTL(fSafe);
		//f = parseCTL(in);

		kripkes[currentKripkeID]->setCTL(*f, currentInitialNode);
		kripkes[currentKripkeID]->newCTLVar(ctlVar);
	}
	void readCTLLine(B& in, Solver& S) {
		++in;
		CTLFormula* fLine = parseCTL(in);
		CTLFormula* fBoth = newCTLFormula();
		fBoth->op = AND;
		fBoth->operand1 = f;
		fBoth->operand2 = fLine;
		f = fBoth;

		kripkes[currentKripkeID]->setCTL(*f, currentInitialNode);
		// printf("\n"); printFormula(f); printf("\n");
	}

	/*CTLFormula* parseCTL(B& in) {
			skipWhitespace(in);
			CTLFormula* f = newCTLFormula();
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
	}*/

	// 	enum CTLOp { ID, NEG, OR, AND, EX, EF, EG, EW, EU, AX, AF, AG, AW, AU};

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
	CTLParser():Parser<B, Solver> ("CTL"){

	}

	bool parseLine(B& in, Solver& S) {

		skipWhitespace(in);
		if (*in == EOF)
			return false;
		else if (*in == 'c') {
			//just a comment
			return false;
		} else if (match(in, "kctlsinglestate")) {
			readCTLSingleState(in, S);
			return true;
		} else if (match(in, "kctlsimp")) {
			skipWhitespace(in);
			readCTLSimp(in, S);
			skipWhitespace(in);
			return true;
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
			readNodeAP(in, S);
			return true;
		} else if (match(in, "kctl")) {
			readCTL(in, S);
			return true;
		} else if (match(in, "and") || match(in, "AND") || match(in, "&") || match(in, "^")) {
			readCTLLine(in, S);
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
