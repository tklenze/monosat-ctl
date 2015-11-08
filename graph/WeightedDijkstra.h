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

#ifndef WEIGHTED_DIJKSTRA_H_
#define WEIGHTED_DIJKSTRA_H_

#include <vector>
#include "dgl/alg/Heap.h"
#include "dgl/DynamicGraph.h"
#include "dgl/Reach.h"
#include <limits>

//Removed this from the dgl, because it is obsoleted by the (now weighted) Dijkstra implementation.

using namespace dgl;

template<class GraphWeight, class Weight = double, bool undirected = false>
class WeightedDijkstra: public Distance<Weight> {
public:
	DynamicGraph<GraphWeight> & g;
	std::vector<Weight> & weights;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	Weight INF;

	std::vector<Weight> old_dist;
	std::vector<int> changed;

	std::vector<Weight> dist;
	std::vector<int> prev;
	struct DistCmp {
		std::vector<Weight> & _dist;
		bool operator()(int a, int b) const {
			return _dist[a] < _dist[b];
		}
		DistCmp(std::vector<Weight> & d) :
				_dist(d) {
		}
		;
	};
	dgl::Heap<DistCmp> q;

public:
	//stats
	
	long stats_full_updates = 0;
	long stats_fast_updates = 0;
	long stats_fast_failed_updates = 0;
	long stats_skip_deletes = 0;
	long stats_skipped_updates = 0;
	long stats_num_skipable_deletions = 0;
	double mod_percentage = 0;

	double stats_full_update_time = 0;
	double stats_fast_update_time = 0;
	WeightedDijkstra(int s, DynamicGraph<GraphWeight> & graph, std::vector<Weight> & weights) :
			g(graph), weights(weights), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
					0), source(s), INF(0), q(DistCmp(dist)) {
		
		mod_percentage = 0.2;
		
		INF = std::numeric_limits<Weight>::infinity();
	}
	
	void setSource(int s) {
		source = s;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	int getSource() {
		return source;
	}
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	void updateFast() {
		stats_fast_updates++;
		assert(last_deletion == g.deletions);
		num_updates++;
		last_modification = g.modifications;
		last_addition = g.additions;
		
		dist.resize(g.nodes());
		prev.resize(g.nodes());
		while (weights.size() <= g.nodes()) {
			weights.push_back(1);
		}
		
		q.clear();
		if (last_history_clear != g.historyclears) {
			history_qhead = 0;
			last_history_clear = g.historyclears;
		}
		//ok, now check if any of the added edges allow for a decrease in distance.
		for (int i = history_qhead; i < g.historySize(); i++) {
			assert(g.getChange(i).addition); //NOTE: Currently, this is glitchy in some circumstances - specifically, ./modsat -rinc=1.05 -rnd-restart  -conflict-shortest-path  -no-conflict-min-cut   -rnd-init -rnd-seed=01231 -rnd-freq=0.01 /home/sam/data/gnf/unit_tests/unit_test_17_reduced.gnf can trigger this assertion!
			int edgeID = g.getChange(i).id;
			int u = g.all_edges[edgeID].from;
			int v = g.all_edges[edgeID].to;
			Weight alt = dist[u] + weights[u];
			if (alt < dist[v]) {
				
				if (dist[v] >= INF) {
					//this was changed
					changed.push_back(v);
				}
				
				dist[v] = alt;
				prev[v] = edgeID;
				
				if (!q.inHeap(v))
					q.insert(v);
				else
					q.decrease(v);
			} else if (undirected) {
				int v = g.all_edges[edgeID].from;
				int u = g.all_edges[edgeID].to;
				Weight alt = dist[u] + weights[u];
				if (alt < dist[v]) {
					
					if (dist[v] >= INF) {
						//this was changed
						changed.push_back(v);
					}
					
					dist[v] = alt;
					prev[v] = edgeID;
					
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
			}
		}
		history_qhead = g.historySize();
		
		while (q.size()) {
			int u = q.removeMin();
			if (dist[u] == INF)
				break;
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				if (!g.edgeEnabled(g.incident(u, i, undirected).id))
					continue;
				
				int edgeID = g.incident(u, i, undirected).id;
				int v = g.incident(u, i, undirected).node;
				Weight alt = dist[u] + weights[u];
				if (alt < dist[v]) {
					if (dist[v] >= INF) {
						//this was changed
						changed.push_back(v);
					}
					
					dist[v] = alt;
					prev[v] = edgeID;
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
				
			}
		}
		
	}
	std::vector<int> & getChanged() {
		return changed;
	}
	void clearChanged() {
		changed.clear();
	}
	
	void drawFull() {
		
	}
	
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;
		if (local_it == 17513) {
			int a = 1;
		}
		if (last_modification > 0 && g.modifications == last_modification)
			return;
		assert(weights.size() >= g.edges());
		/*		while(weight.size()<=g.nodes()){
		 weight.push_back(1);
		 }*/
		/*if (last_addition==g.additions && last_modification>0){
		 //if none of the deletions were to edges that were the previous edge of some shortest path, then we don't need to do anything
		 if(last_history_clear!=g.historyclears){
		 history_qhead=0;
		 last_history_clear=g.historyclears;
		 }
		 bool need_recompute = false;
		 //ok, now check if any of the added edges allow for a decrease in distance.
		 for (int i = history_qhead;i<g.history.size();i++){
		 assert(!g.getChange(i).addition);
		 int u=g.getChange(i).u;
		 int v=g.getChange(i).v;
		 if(prev[v]==u){
		 history_qhead = i-1;
		 need_recompute=true;
		 //this deletion matters, so we need to recompute.
		 break;
		 }
		 }
		 if(!need_recompute){
		 //none of these deletions touched any shortest paths, so we can ignore them.

		 last_modification=g.modifications;
		 last_deletion = g.deletions;
		 last_addition=g.additions;

		 history_qhead=g.history.size();
		 last_history_clear=g.historyclears;

		 assert(dbg_uptodate());
		 stats_skip_deletes++;
		 return;
		 }
		 }*/

		/*if(last_deletion==g.deletions && last_modification>0  ){
		 //Don't need to do anything at all.
		 if(last_addition==g.additions){
		 last_modification = g.modifications;
		 stats_skipped_updates++;
		 assert(dbg_uptodate());
		 return;
		 }
		 //we can use a faster, simple dijkstra update method
		 updateFast();
		 assert(dbg_uptodate());
		 return;
		 }*/

		stats_full_updates++;
		
		dist.resize(g.nodes(), INF);
		prev.resize(g.nodes());
		old_dist.resize(g.nodes(), INF);
		q.clear();
		for (int i = 0; i < g.nodes(); i++) {
			old_dist[i] = last_modification > 0 ? dist[i] : INF; //this won't work properly if we added nodes...
			dist[i] = INF;
			prev[i] = -1;
		}
		dist[source] = 0;
		q.insert(source);
		while (q.size()) {
			int u = q.peakMin();
			if (dist[u] == INF)
				break;
			if (old_dist[u] >= INF) {
				changed.push_back(u);
			}
			q.removeMin();
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				int edgeID = g.incident(u, i, undirected).id;
				if (!g.edgeEnabled(edgeID))
					continue;
				
				int v = g.incident(u, i, undirected).node;
				
				Weight alt = dist[u] + weights[edgeID];
				if (alt < dist[v]) {
					dist[v] = alt;
					prev[v] = edgeID;
					if (!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
			}
		}
		
		for (int u = 0; u < g.nodes(); u++) {
			//while(q.size()){
			//iterate through the unreached nodes and check which ones were previously reached
			
			if (last_modification <= 0 || (old_dist[u] < INF && dist[u] >= INF)) {
				changed.push_back(u);
			}
		}
		//}
		assert(dbg_uptodate());
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
		;
	}
	bool dbg_path(int to) {
#ifdef DEBUG_DIJKSTRA
		assert(connected(to));
		if(to == source) {
			return true;
		}
		int p = prev[to];

		if(p<0) {
			return false;
		}
		if(p==to) {
			return false;
		}

		return dbg_path(p);

#endif
		return true;
	}
	bool dbg_uptodate() {
		
		return true;
	}
	
	bool connected_unsafe(int t) {
		return t < dist.size() && dist[t] < INF;
	}
	bool connected_unchecked(int t) {
		assert(last_modification == g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t) {
		if (last_modification != g.modifications)
			update();
		
		assert(dbg_uptodate());
		
		return dist[t] < INF;
	}
	Weight & distance(int t) {
		if (last_modification != g.modifications)
			update();
		return dist[t];
	}
	Weight & distance_unsafe(int t) {
		if (connected_unsafe(t))
			return dist[t];
		else
			return INF;
	}
	int incomingEdge(int t) {
		assert(t >= 0 && t < prev.size());
		assert(prev[t] >= -1);
		return prev[t];
	}
	int previous(int t) {
		if (prev[t] < 0)
			return -1;
		if (undirected && g.getEdge(incomingEdge(t)).from == t) {
			return g.getEdge(incomingEdge(t)).to;
		}
		assert(g.getEdge(incomingEdge(t)).to == t);
		return g.getEdge(incomingEdge(t)).from;
	}
	
};

#endif /* DIJKSTRA_H_ */
