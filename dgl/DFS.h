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

#ifndef DFS_REACHABILITY_H_
#define DFS_REACHABILITY_H_

#include <vector>
#include "alg/Heap.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "Reach.h"

namespace dgl {

template<typename Weight,class Status = Reach::NullStatus, bool undirected = false>
class DFSReachability: public Reach {
public:
	
	DynamicGraph<Weight> & g;
	Status & status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	int INF;

	bool opt_skip_deletions = false;
	bool opt_skip_additions = false;
	bool opt_inc_graph = false;
	int opt_dec_graph = 0;

	std::vector<int> q;
	std::vector<int> check;
	const int reportPolarity;

	//std::vector<char> old_seen;
	std::vector<char> seen;
//	std::vector<int> changed;
	
	std::vector<int> prev;
	
	long stats_full_updates=0;
	long stats_fast_updates=0;
	long stats_fast_failed_updates=0;
	long stats_skip_deletes=0;
	long stats_skipped_updates=0;
	long stats_num_skipable_deletions=0;


	double stats_full_update_time=0;
	double stats_fast_update_time=0;
public:


	DFSReachability(int s, DynamicGraph<Weight> & graph, Status & _status= Reach::nullStatus, int _reportPolarity = 0) :
			g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
					0), source(s), INF(0), reportPolarity(_reportPolarity) {
	}
	DFSReachability(DynamicGraph<Weight> & graph,Status & _status= Reach::nullStatus, int _reportPolarity = 0):			g(graph), status(_status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
			0), source(0), INF(0), reportPolarity(_reportPolarity) {
	}
	//Connectivity(const Connectivity& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),mod_percentage(0.2),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};
	
	void setSource(int s) {
		source = s;
		last_modification = -1;
		last_addition = -1;
		last_deletion = -1;
	}
	int getSource() {
		return source;
	}
	
	/*void updateFast(){
	 stats_fast_updates++;
	 

	 assert(last_deletion==g.deletions);
	 last_modification=g.modifications;
	 last_addition=g.additions;
	 INF=g.nodes()+1;
	 seen.resize(g.nodes());
	 prev.resize(g.nodes());

	 if(lastaddlist!=g.addlistclears){
	 addition_qhead=0;
	 lastaddlist=g.addlistclears;
	 }
	 int start = q.size();
	 //ok, now check if any of the added edges allow for new connectivity
	 for (int i = addition_qhead;i<g.addition_list.size();i++){
	 int u=g.addition_list[i].u;
	 int v=g.addition_list[i].v;

	 if(!seen[v]){
	 q.push_back(v);
	 seen[v]=1;
	 prev[v]=u;
	 }
	 }
	 addition_qhead=g.addition_list.size();

	 for(int i = start;i<q.size();i++){
	 int u = q[i];
	 assert(seen[u]);
	 for(int i = 0;i<g.nIncident(u);i++){
	 int v = g.incident(u,i);

	 if(!seen[v]){
	 //this was changed
	 changed.push_back(v);
	 seen[v]=1;
	 prev[v]=u;
	 q.push_back(v);
	 }
	 }
	 }
	 
	 }*/
	/*	std::vector<int> & getChanged(){
	 return changed;
	 }
	 void clearChanged(){
	 changed.clear();
	 }*/

	/*
	 * WARNING: THIS FUNDAMENTALLY WONT WORK if there are any cycles in the graph!
	 * inline void delete_update(int to){
	 q.clear();
	 q.push_back(to);
	 seen[to]=0;
	 //Is this really safe? check it very carefully, it could easily be wrong
	 while(q.size()){
	 int u = q.back();
	 q.pop_back();
	 assert(!seen[u]);
	 for(int i = 0;i<g.nIncoming(u);i++){
	 int v = g.incoming(u,i);
	 if(seen[v]){
	 seen[v]=1;
	 //Then since to is still seen, we are up to date
	 break;
	 }
	 }
	 if(!seen[u]){
	 for(int i = 0;i<g.nIncident(u);i++){
	 int v = g.incident(u,i);
	 if(seen[v] && prev[v]==to){
	 seen[v]=0;
	 }
	 }
	 }else{
	 #ifdef GRAPH_DEBUG
	 for(int i = 0;i<g.nIncident(u);i++){
	 int v = g.incident(u,i);
	 assert(seen[v]);
	 }
	 #endif
	 }
	 }
	 }*/

	void setNodes(int n) {
		q.reserve(n);
		check.reserve(n);
		seen.resize(n);
		prev.resize(n);
		INF = g.nodes() + 1;
	}
	
	inline void add_update(int to, bool update) {
		q.clear();
		q.push_back(to);
		//while(q.size()){
		//auto & adjacency = undirected? g.adjacency_undirected:g.adjacency;
		while (q.size()) {
			int u = q.back();
			q.pop_back();
			assert(seen[u]);
			if (update)
				status.setReachable(u, seen[u]);
			//status.setReachable(u,true);
			//if(!old_seen[u]){
			//	changed.push_back(u);
			//}
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				if (!g.edgeEnabled(g.incident(u, i, undirected).id))
					continue;
				int v = g.incident(u, i, undirected).node;
				int edgeID = g.incident(u, i, undirected).id;
				if (!seen[v]) {
					seen[v] = 1;
					prev[v] = edgeID;
					q.push_back(v);
				}
			}
		}
	}
	
	inline void delete_update(int to) {
		q.clear();
		q.push_back(to);
		seen[to] = false;
		prev[to] = -1;
		check.clear();
		check.push_back(to);
		//while(q.size()){
		//auto & adjacency = undirected? g.adjacency_undirected:g.adjacency;
		while (q.size()) {
			int u = q.back();
			q.pop_back();
			assert(!seen[u]);
			
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				/*				if(!g.edgeEnabled( g.incident(u,i).id))
				 continue;*/
				int v = g.incident(u, i, undirected).node;
				int edgeID = g.incident(u, i, undirected).id;
				if (seen[v] && previous(v) == u) {
					seen[v] = 0;
					prev[v] = -1;
					check.push_back(v);
					q.push_back(v);
					
				}
			}
		}
		
		q.clear();
		for (int i = 0; i < check.size(); i++) {
			int u = check[i];
			if (!seen[u]) {
				for (int i = 0; i < g.nIncoming(u, undirected); i++) {
					if (g.edgeEnabled(g.incoming(u, i, undirected).id)) {
						int from = g.incoming(u, i, undirected).node;
						int to = u;
						int edgeID = g.incoming(u, i, undirected).id;
						if (seen[from]) {
							seen[to] = 1;
							prev[to] = edgeID;
							add_update(to, false);
							break;
						}
					}
				}
			}
		}
		
	}
	bool update_additions() {
		
		if (g.historyclears != last_history_clear) {
			last_history_clear = g.historyclears;
			history_qhead = 0;
		}
		
		assert(INF > g.nodes());
		assert(seen.size() >= g.nodes());
		//old_seen.resize(g.nodes());
		q.clear();
		
		for (int i = history_qhead; i < g.historySize(); i++) {
			int edgeid = g.getChange(i).id;
			int from = g.getEdge(edgeid).from;
			int to = g.getEdge(edgeid).to;
			if (g.getChange(i).addition) {
				//incrementally add edge
				if (seen[from] && !seen[to]) {
					seen[to] = 1;
					prev[to] = edgeid;
					add_update(to, true);
				} else if (undirected) {
					if (seen[to] && !seen[from]) {
						seen[from] = 1;
						prev[from] = edgeid;
						add_update(from, true);
					}
				}
				
			} else if (!undirected
					&& (to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))) {
				//then deleting this edge has no impact on connectivity, so don't need to do anything
				
			} else if (undirected
					&& (to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))
					&& (from == source || !seen[to] || (seen[from] && seen[previous(to)] && previous(to) != to))) {
				
			} else {
				stats_fast_failed_updates++;
				
				return false;
			}
			
		}
		
		stats_fast_updates++;
		

		assert(dbg_uptodate());
		
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
		return true;
	}
	
	bool incrementalUpdate() {
		
		if (g.historyclears != last_history_clear) {
			last_history_clear = g.historyclears;
			history_qhead = 0;
		}
		
		assert(INF > g.nodes());
		assert(seen.size() >= g.nodes());
		//old_seen.resize(g.nodes());
		q.clear();
		
		for (int i = history_qhead; i < g.historySize(); i++) {
			int edgeid = g.getChange(i).id;
			int from = g.getEdge(edgeid).from;
			int to = g.getEdge(edgeid).to;
			
			if (g.getChange(i).addition && g.edgeEnabled(g.getChange(i).id)) {
				//incrementally add edge
				
				if (seen[from] && !seen[to]) {
					seen[to] = 1;
					prev[to] = edgeid;
					add_update(to, false);
				} else if (undirected) {
					if (seen[to] && !seen[from]) {
						seen[from] = 1;
						prev[from] = edgeid;
						add_update(from, false);
					}
				}
				
			} else if (!g.getChange(i).addition && !g.edgeEnabled(g.getChange(i).id)) {
				
				if ((to == source || !seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from))
						&& (!undirected
								|| (from == source || !seen[to]
										|| (seen[from] && seen[previous(to)] && previous(to) != to)))) {
					//then deleting this edge has no impact on connectivity, so don't need to do anything
				} else {
					delete_update(to);
				}
			}
		}
		
		for (int u = 0; u < g.nodes(); u++) {
			status.setReachable(u, seen[u]);
		}
		
		stats_fast_updates++;

		assert(dbg_uptodate());
		
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		last_history_clear = g.historyclears;
		
		return true;
	}
	
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;
		
		if (last_modification > 0 && g.modifications == last_modification) {
			stats_skipped_updates++;
			return;
		}
		
		if (last_modification > 0 && last_deletion == g.deletions) {
			stats_num_skipable_deletions++;
			if (opt_skip_deletions && reportPolarity < 1) {
				return;		//I don't trust the correctness of these shortcuts
			}
		}
		if (last_modification > 0 && last_addition == g.additions) {
			if (opt_skip_additions && reportPolarity > -1) {
				return;		//I don't trust the correctness of these shortcuts
			}
		}
		
		setNodes(g.nodes());
		
		if (g.historyclears != last_history_clear) {
			last_history_clear = g.historyclears;
			history_qhead = 0;
		} else if (opt_inc_graph && last_modification > 0 && (g.historyclears <= (last_history_clear + 1))) {// && (g.history.size()-history_qhead < g.edges()*mod_percentage)){
			if (opt_dec_graph == 2) {
				if (incrementalUpdate())
					return;
			} else {
				if (opt_dec_graph == 1 && last_deletion < g.deletions) {
					
					//scan through the deletions and check if any of them matter..
					bool safe = true;
					for (int i = history_qhead; i < g.historySize(); i++) {
						int edgeid = g.getChange(i).id;
						int from = g.getEdge(edgeid).from;
						int to = g.getEdge(edgeid).to;
						if (g.getChange(i).addition) {
							//safe
						} else if (!seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from)) {
							//then deleting this edge has no impact on connectivity, so don't need to do anything
						} else {
							safe = false;
							break;
						}
						if (undirected) {
							int from = g.getEdge(edgeid).to;
							int to = g.getEdge(edgeid).from;
							if (g.getChange(i).addition) {
								//safe
							} else if (!seen[from] || (seen[to] && seen[previous(to)] && previous(to) != from)) {
								//then deleting this edge has no impact on connectivity, so don't need to do anything
							} else {
								safe = false;
								break;
							}
						}
						
					}
					if (safe) {
						last_deletion = g.deletions;
					}
					
				}
				
				if (last_deletion == g.deletions) {
					if (update_additions())
						return;
				}
			}
			/**/
		}
		
		stats_full_updates++;
		
		q.clear();
		for (int i = 0; i < g.nodes(); i++) {
			seen[i] = 0;
			prev[i] = -1;
		}
		seen[source] = 1;
		//auto & adjacency = undirected? g.adjacency_undirected:g.adjacency;
		q.push_back(source);
		while (q.size()) {
			int u = q.back();
			q.pop_back();
			assert(seen[u]);
			if (reportPolarity == 1)
				status.setReachable(u, true);
			
			for (int i = 0; i < g.nIncident(u, undirected); i++) {
				if (!g.edgeEnabled(g.incident(u, i, undirected).id))
					continue;
				int v = g.incident(u, i, undirected).node;
				int edgeid = g.incident(u, i, undirected).id;
				if (!seen[v]) {
					seen[v] = 1;
					prev[v] = edgeid;
					q.push_back(v);
				}
			}
		}
		
		if (reportPolarity < 1) {
			for (int u = 0; u < g.nodes(); u++) {
				if (!seen[u]) {
					status.setReachable(u, false);
				} else if (reportPolarity == 0) {
					status.setReachable(u, true);
				}
			}
		}
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
		int p = previous(to);

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
	void drawFull() {
		printf("digraph{\n");
		for (int i = 0; i < g.nodes(); i++) {
			
			if (seen[i]) {
				printf("n%d [fillcolor=blue style=filled]\n", i);
			} else {
				printf("n%d \n", i);
			}
			
		}
		
		for (int i = 0; i < g.nodes(); i++) {
			for (int j = 0; j < g.nIncident(i); j++) {
				int id = g.incident(i, j).id;
				int u = g.incident(i, j).node;
				const char * s = "black";
				if (g.edgeEnabled(id))
					s = "blue";
				else
					s = "red";
				
				printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
			}
		}
		
		printf("}\n");
	}
	bool dbg_uptodate() {
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
		return true;
		UnweightedDijkstra<> d(source,g);
		d.update();
		//drawFull();
		for(int i = 0;i<g.nodes();i++) {

			int dbgdist = d.dist[i];
			if(!seen[i]) {
				if(reportPolarity<1)
				assert(dbgdist==d.INF );
			} else {
				if(reportPolarity>-1) {
					if(!(dbgdist<d.INF)) {
						drawFull();
					}

					assert(dbgdist<d.INF);
				}
			}
		}
#endif
		return true;
	}
	
	bool connected_unsafe(int t) {
		return t < seen.size() && seen[t];
	}
	bool connected_unchecked(int t) {
		assert(last_modification == g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t) {
		if (last_modification != g.modifications)
			update();
		
		assert(dbg_uptodate());
		
		return seen[t];
	}
	int distance(int t) {
		if (connected(t))
			return 1;
		else
			return -1;
	}
	int distance_unsafe(int t) {
		if (connected_unsafe(t))
			return 1;
		else
			return -1;
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
}
;
#endif
