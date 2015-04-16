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
		printf("CTLParser.readKripke: Adding %d nodes...\n", n);  fflush(stdout);
		kripke->newNodes(n);
		printf("CTLParser.readKripke: Success\n");  fflush(stdout);
		kripkes[g] = kripke;
		S.addTheory(kripke);

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
		/*
		if (isNumber(in)) {
			// TODO
			// Do nothing for now, should set Node-AP to true or false, depending on input
		} else if (in.match('v')) {
			int nodeVar = parseInt(in);
			// TODO
			// should set it to undecided
		}
		*/
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
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", kripkeID, edgeVar), exit(1);
		}
		if (edgeVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", edgeVar), exit(1);
		}
		while (edgeVar >= S.nVars())
			S.newVar();

		skipWhitespace(in);
		if(*in=='\n' || *in==0){
			/*
			if (kripkes[kripkeID]) {
				// TODO not implemented yet
				//kripkes[kripkeID]->newEdge(from, to, edgeVar);
			} else {
				printf("PARSE ERROR! Undeclared kripke identifier %d for edge %d\n", kripkeID, edgeVar), exit(1);
				exit(1);
			}
			*/
		}
	}

	void readReach(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		//reach grachID u w var is a reach query: var is true if can u reach w in graph g, false otherwise
		
		++in;
		
		int kripkeID = parseInt(in);
		int from = parseInt(in);
		// int steps = parseInt(in);
		int to = parseInt(in);
		int reachVar = parseInt(in) - 1;
		if (kripkeID < 0 || kripkeID >= kripkes.size()) {
			printf("PARSE ERROR! Undeclared graph identifier %d for edge %d\n", kripkeID, reachVar), exit(1);
		}
		if (reachVar < 0) {
			printf("PARSE ERROR! Edge variables must be >=0, was %d\n", reachVar), exit(1);
		}

		while (reachVar >= S.nVars())
			S.newVar();
		
		if (kripkes[kripkeID]) {
			// FIXME this, of course, does not work right now, since we have not implemented reaches
			//kripkes[kripkeID]->reaches(from, to, reachVar);
		} else {
			printf("PARSE ERROR! Undeclared graph identifier %d\n", kripkeID), exit(1);
			exit(1);
		}
	}

	void readCTL(B& in, Solver& S) {
		if (opt_ignore_theories) {
			skipLine(in);
			return;
		}
		// TODO

		// this is just test code
		printf("Test code! Printing Kripke structures\n");
		for (int i = 0; i<kripkes.size(); i++) {
			printf("\n\nKripke structure:\n");
			kripkes[i]->g_over->draw();
		}

		return;
	}

public:
	bool parseLine(B& in, Solver& S) {
		printf("Parse line: ");

		for (int i = 0; in[i] != '\0'; i++) {
			printf("%c", in[i]);
		}
		printf("\n");


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
		} else if (match(in, "edge")) {
			edgeCount++;
			readEdge(in, S);
			return true;
		} else if (match(in, "nodeap")) {
			nodeCount++;
			readNodeAP(in, S);
			return true;
		}else if (match(in, "ctl")) {
			readCTL(in, S);
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
