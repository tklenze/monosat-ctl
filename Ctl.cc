/*****************************************************************************************[Main.cc]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless
 Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 Copyright (c) 2007-2010, Niklas Sorensson

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
#include <cstddef>
#include <gmpxx.h>
#include <fstream>
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <signal.h>
#include <zlib.h>
#include <sstream>
#include <string>
#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "graph/GraphParser.h"
#include "fsm/FSMParser.h"
#include "core/Dimacs.h"
#include "core/AssumptionParser.h"
#include "fsm/LSystemParser.h"
#include "core/Solver.h"
#include "core/Config.h"
#include <unistd.h>
#include <sys/time.h>
#include <algorithm>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include "simp/SimpSolver.h"
#include "pb/PbTheory.h"
#include "pb/PbParser.h"
#include "mtl/Map.h"
#include "graph/GraphTheory.h"
#include "geometry/GeometryTheory.h"
#include "geometry/GeometryParser.h"
#include "ctl/DynamicKripke.h"



using namespace Monosat;
using namespace std;

int main(int argc, char** argv) {


	// git test change


	DynamicKripke myKripke;
	std::vector<int> states;
	states.push_back(true);
	states.push_back(false);

	//statelabel[stateid][labelid] = false;

	int node0 = myKripke.addState(states);
	states.clear();
	states.push_back(false);
	states.push_back(true);
	int node1 = myKripke.addState(states);
	states.clear();
	states.push_back(true);
	states.push_back(true);
	int node2 = myKripke.addState(states);
	states.clear();
	states.push_back(false);
	states.push_back(false);
	int node3 = myKripke.addState(states);
	myKripke.addTransition(node1, node2, -1, true);
	myKripke.addTransition(node2, node3, -1, true);
	myKripke.addTransition(node1, node3, -1, true);
	myKripke.disableTransition(2); // node1 -> node3
	myKripke.draw(node1, node3);


	int nodecount = myKripke.nodes();


}
