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

#ifndef CTLTarjansSCC_H_
#define CTLTarjansSCC_H_

#include <vector>
#include "dgl/alg/Heap.h"
#include "dgl/DynamicGraph.h"
#include "core/Config.h"
#include "dgl/ConnectedComponents.h"
#include "dgl/alg/DisjointSets.h"
#include <limits>
#include <algorithm>

namespace dgl {
template<typename Weight, class Status = ConnectedComponents::NullConnectedComponentsStatus>
class CTLTarjansSCC: public ConnectedComponents {
public:

	DynamicGraph<Weight> & g;
	DynamicKripke & k;
	Status & status;
	int last_modification;
	int last_k_modification;

	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;
	int INF;

	Bitset* inner;

	struct Component {
		int id = -1; //the unique id of the component
		int next = -1; //pointer to the next node in the component
		Component() {
		}
		Component(int id, int next) :
				id(id), next(next) {
		}
	};

	std::vector<Component> scc; //the components of each node
	struct SCC {
		int sz; //size of the component
		int element; //pointer to an arbitrary node in the scc.
	};
	std::vector<SCC> scc_set;//this scc_set treats individual nodes as sccs
	std::vector<int> strict_scc_set;//this lists only the 'real' strongly connected components, excluding lone vertices or vertices in acyclic components
	struct QNode{
		int node;
		int edgeID;
	};
	std::vector<QNode> q;

	std::vector<int> indices;
	std::vector<int> lowlink;
	std::vector<bool> in_q;
	int last_edgeID = -1;
	int stats_full_updates = 0;
	int stats_fast_updates = 0;
	int stats_fast_failed_updates = 0;
	int stats_skip_deletes = 0;
	int stats_skipped_updates = 0;
	int stats_num_skipable_deletions = 0;
	double mod_percentage = 0.2;

	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;

public:
	CTLTarjansSCC(DynamicGraph<Weight> & graph, DynamicKripke & kripke) :
		g(graph), k(kripke), status(nullConnectedComponentsStatus), last_k_modification(-1), last_modification(-1), last_addition(-1), last_deletion(
					-1), history_qhead(0), last_history_clear(0), INF(0) {

	}

	CTLTarjansSCC(DynamicGraph<Weight> & graph, Status & _status, DynamicKripke & kripke) :
		g(graph), k(kripke), status(_status), last_modification(-1), last_k_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
					0), INF(0) {

	}

	void setNodes(int n) {
		q.reserve(n);
		in_q.resize(n);
		indices.resize(n, -1);
		lowlink.resize(n, -1);
		scc.resize(n);
		INF = std::numeric_limits<int>::max();
	}

	void strongConnect(int node, int fromEdge, int & index, std::vector<int> * scc_out = nullptr) {
		//Following wikipedia's pseudocode, from http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
		assert(indices[node] < 0);
		assert(lowlink[node] < 0);
		indices[node] = index;
		lowlink[node] = index;
		index++;
		assert(!in_q[node]);
		q.push_back({node,fromEdge});
		in_q[node] = true;

		//for(auto & edge:g.adjacency[node]){
		for (int i = 0; i < g.nIncident(node); i++) {
			auto & edge = g.incident(node, i);
			int edgeID = edge.id;
			if (suitable(edgeID)) {
				int to = edge.node;

				if (indices[to] < 0) {
					strongConnect(to,edgeID, index, scc_out);
					lowlink[node] = std::min(lowlink[node], lowlink[to]);
				} else if (in_q[to]) {
					last_edgeID=edgeID;
					lowlink[node] = std::min(lowlink[node], indices[to]);
				}
			}
		}

		// If v is a root node, pop the stack and generate an SCC
		// Do not so if fromEdge = -1, i.e. if this is the highest call in the hierarchy. This would mean that no edge leads here, i.e. a SCC with just the initial node and no incoming transition from itself.
		if (opt_verb > 2) {
		//	if (fromEdge < 0)
		//		printf("Tarjan: Ignoring this node %d, since fromEdge is negative: %d\n", node, fromEdge);
		//	else
				printf("Tarjan: node %d, fromEdge %d (%d -> %d). Working on overapprox? %d\n", node, fromEdge, k.getEdge(fromEdge).from, k.getEdge(fromEdge).to, k.isover);
		}
		if (lowlink[node] == indices[node]) {// && fromEdge >= 0
			int sccID = scc_set.size();
			int sz = 0;
			assert(q.size());
			int first = q.back().node;
			int curEdge = q.back().edgeID;

			if (scc_out) {
				scc_out->clear();
			}
			do {
				int n = q.back().node;
				int curEdge = q.back().edgeID;
				sz++;
				q.pop_back();
				assert(in_q[n]);
				in_q[n] = false;
				int next;
				if (q.size() && n != node) {
					next = q.back().node;
				} else {
					next = first;
				}
				scc[n]= {sccID,next};
				if (scc_out && (n != node)) {
					if(curEdge!=-1 && suitable(curEdge))
						scc_out->push_back(curEdge);
					else{
						int a=1;
					}
				}
				if(n == node){
					if (scc_out&& sz >1){
						//this is a bit of a hack, but because I haven't figured out how to properly keep track of the last connected edgeID,
						//I need to do a quick search for it here (fix this!).

						for(int i = 0;i<g.nIncoming(n);i++){
							int edgeID = g.incoming(n,i).id;
							if(suitable(edgeID) && g.incoming(n,i).node ==first){
								scc_out->push_back(edgeID);
								break;
							}
						}
					}
					break;
				}

			} while (q.size());
			scc_set.push_back( { sz, first });

			if(sz>1){
				strict_scc_set.push_back(sccID);
			}
		}

	}

	// An edge is suitable if its origin and destination satisfies the formula and it is enabled
	bool suitable(int edgeID) {
		if ((k.edgeEnabled(edgeID)) && inner->operator [](k.getEdge(edgeID).from) && inner->operator [](k.getEdge(edgeID).to)) {
			if (opt_verb > 2) {
				printf("  Tarjan: %d -> %d is suitable (isover? %d) (enabled? %d) (edgeID: %d)\n", k.getEdge(edgeID).from, k.getEdge(edgeID).to, k.isover ,k.edgeEnabled(edgeID), edgeID);
				printf("Enabled transitions:");
				for (int i=0; i<k.transitions.size(); i++)
					if (k.transitions[i])
						printf(" %d (%d -> %d),", i, k.getEdge(i).from, k.getEdge(i).to);
			}
			return true;
		} else {
			if (opt_verb > 2)
				printf("  Tarjan: %d -> %d is NOT suitable (isover? %d) (enabled? %d) (edgeID: %d)\n", k.getEdge(edgeID).from, k.getEdge(edgeID).to, k.isover , k.edgeEnabled(edgeID), edgeID);
			return false;
		}
	}

	void setInnerFormula (Bitset& in) {
		inner = &in;
	}

	// Force update, no caching
	void updateForced() {

		static int iteration = 0;
		int local_it = ++iteration;

		stats_full_updates++;

		if (last_deletion == g.deletions) {
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes());
		scc_set.clear();
		strict_scc_set.clear();
#ifndef NDEBUG
		for(int i = 0;i<in_q.size();i++)
			assert(!in_q[i]);
#endif
		q.clear();
		for (int i = 0; i < g.nodes(); i++) {
			indices[i] = -1;
			lowlink[i]=-1;
		}
		int index = 0;
		for (int i = 0; i < g.nodes(); i++) {
			if (indices[i] < 0 && inner->operator [](i)) { // restrict graph to states satisfying the inner formula
				strongConnect(i,-1, index);
			}
		}


#ifndef NDEBUG
		for(int i = 0;i<in_q.size();i++)
			assert(!in_q[i]);
#endif
		status.setComponents(scc_set.size());
		last_modification = g.modifications;
		last_k_modification = k.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;

		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		;
	}

	// makes use of caching
	void update() {
		if (last_modification > 0 && k.modifications > 0 && g.modifications == last_modification && k.modifications == last_k_modification) {
			stats_skipped_updates++;
			return;
		}
		updateForced();
	}

	bool connected(int from, int to) {
		update();
		return scc[from].id == scc[to].id;
	}

	int numStrictSCCs(){
		update();
		return strict_scc_set.size();
	}

	std::vector<int> & getStrictSCCs(){
		return strict_scc_set;
	}

	int numComponents() {
		update();
		return scc_set.size();
	}
	int getComponentUnsafe(int node) {
			return scc[node].id;
		}
	int getComponent(int node) {
		update();
		return scc[node].id;
	}
	int getElement(int sccID) {
		update();
		assert(sccID < scc_set.size());
		return scc_set[sccID].element;
	}
	//Returns the next node in the cyclic linked list representing this scc (in arbitrary order)
	int nextNode(int node){
		update();
		return scc[node].next;
	}
	//Returns the next node in the cyclic linked list representing this scc (in arbitrary order). Update() must have been called before this method.
	int nextNodeUnsafe(int node){
		return scc[node].next;
	}
	//get all the nodes strongly connected to this node (including this one)
	void getConnectedComponent(int forNode, std::vector<int> & store){
		update();
		int sccID = getComponentUnsafe(forNode);
		store.clear();
		store.push_back(forNode);
		int n = nextNodeUnsafe(forNode);
		while(n!=forNode){
			assert(sccID == getComponentUnsafe(n));
			store.push_back(n);
			n = nextNodeUnsafe(n);
		}
	}
	int getComponentSize(int sccID) {
		update();
		assert(sccID < scc_set.size());
		return scc_set[sccID].sz;
	}

	bool dbg_uptodate() {
		return true;
	}
	;

	//Compute the SCC for a single node.
	//This is faster for one-shot scc computations from a single source than doing a full update().
	static void getSCC(int node, DynamicGraph<Weight> & graph, std::vector<int>&scc) {
		CTLTarjansSCC<Weight> s(graph);
		s.setNodes(graph.nodes());
		int index = 0;
		s.strongConnect(node,-1, index, &scc);
		assert(scc.size());

	}

};

}
;

#endif /* CTLTarjansSCC_H_ */
