/*
 * DynamicKripke.h
 *
 * Created on: March, 2015
 * 	   Author: Tobias
 *
 * Based on DynamicFSM.h
 *  Created on: Dec 15, 2014
 *      Author: sam
 *
 */

#ifndef DynamicKripke_H_
#define DynamicKripke_H_

#include <vector>
#include "mtl/Vec.h"
#include "mtl/Bitset.h"
#include <algorithm>
#include <cassert>
#include "alg/NFATypes.h"
#include "dgl/DynamicGraph.h"
#include <string>
using namespace dgl;
namespace Monosat {

/*
 * setAPCount(int) has to be called BEFORE adding any states during initialization.
 *
 */

class DynamicKripke{

public:
	DynamicGraph<int> g;
	//std::vector<Bitset> edge_status;
	int id;
	bool is_changed = true;
	bool is_generator=true;
	bool is_acceptor=true;
	bool is_linear=true;
public:
	std::vector<bool> transitions;

	// This represents the labels of states, which are sets of atomic propositions
	// <s>Access via statelabel[stateID][AP]</s>
	std::vector<Bitset*> statelabel;
	int apcount = 0;

	bool adaptive_history_clear = false;
	long historyClearInterval = 1000;
	int modifications=0;
	int additions=0;
	int deletions=0;
	long historyclears=0;
	struct EdgeChange {
		bool addition;
		bool APchange;

		int id; // Either EdgeID or StateID, depending on whether it is an
				// Edge change (APchange==false) or a change in the AP assignment.
		int ap; // Ignored if it is not an APchange, if it is an APchange, this will be the AP id.
		int mod;
		// int prev_mod; // not needed, apparently

	};
	std::vector<EdgeChange> history;

public:
	DynamicKripke(int id = -1):id(id) {

	}

	~DynamicKripke() {

	}
	bool isGenerator()const{
		return is_generator;
	}
	bool isAcceptor()const{
		return is_acceptor;
	}
	bool isLinear()const{
		return is_linear;
	}

	int getID(){
		return id;
	}

	bool transitionEnabled(int edgeID)const{
			return transitions[edgeID];
	}


	void setStateLabel(int state, Bitset& label) {
		statelabel[state] = &label;
	}

	Bitset* getStateLabel(int state) {
		return statelabel[state];
	}

	int isAPinStateLabel(int state, int ap) {
		assert(state < statelabel.size());
		assert(ap < statelabel[state]->size());
		return statelabel[state]->operator [](ap);
	}

	void enableAPinStateLabel(int state, int ap) {
		//printf("Called enableAP on state %d, ap %d\n", state, ap);

		if (state >= statelabel.size()) {
			statelabel.resize(state+1);
		}
		if (ap >= statelabel[state]->size()) {
			statelabel[state]->growTo(ap+1);
		}
		assert(state < statelabel.size());
		assert(ap < statelabel[state]->size());
		if (!statelabel[state]->operator [](ap)) {
			statelabel[state]->set(ap);
			modifications++;
			additions = modifications;
			history.push_back( { true, true, state, ap, modifications });
		}
	}

	void disableAPinStateLabel(int state, int ap) {
		//printf("Called disableAP on state %d, ap %d\n", state, ap);

		if (state >= statelabel.size()) {
			statelabel.resize(state+1);
		}
		if (ap >= statelabel[state]->size()) {
			statelabel[state]->growTo(ap+1);
		}
		assert(state < statelabel.size());
		assert(ap < statelabel[state]->size());
		if (statelabel[state]->operator [](ap)) {
			statelabel[state]->clear(ap);
			modifications++;
			deletions = modifications;
			history.push_back( { true, true, state, ap, modifications });
		}
	}

	// Note that this is not reliable information, vectors might have different
	// lengths if initialization messes up
	void setAPCount(int apc) {
		apcount = apc;
	}
	int nProperties(){
		return this->apcount;
	}

	// Is ap in the label of state?
	// 1: Yes, 0: No, -1: AP not present
	//// Commented out, we use booleans in statelabel
	/*
	int stateAPAssignment(int state, int ap) {
		if (statelabel[state].size() > ap) {
			if (statelabel[state][ap])
				return 1;
			else
				return 0;
		}
		return -1;
	}*/


	int addTransition(int from, int to,int edgeID, bool defaultEnabled=true){
		while(from>=g.nodes() || to>=g.nodes()) {
			assert(false); // for now, we just don't want this to happen
			addEmptyState();
		}
		if(edgeID==-1){
			edgeID = g.addEdge(from, to, edgeID);
		}
		transitions.resize(edgeID+1);
		transitions[edgeID] = defaultEnabled;

		return edgeID;
	}

	void enableTransition(int edgeID) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		if (!transitions[edgeID]) {
			transitions[edgeID] = true;
			modifications++;
			additions = modifications;
			history.push_back( { true, false, edgeID, 0, modifications });
		}
	}
	void disableTransition(int edgeID) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		if (transitions[edgeID]) {
			transitions[edgeID] = false;
			modifications++;
			history.push_back( { false, false, edgeID, 0, modifications });
			deletions = modifications;
		}
	}
	int addEmptyState(){
		Bitset *empty = new Bitset(apcount);
		empty->zero();
		return addNode(*empty);
	}
	int addState(Bitset& v){
		return addNode(v);
	}
	int addNode(Bitset& v) {
		g.addNode();
		modifications++;
		additions = modifications;
		deletions = modifications;
		statelabel.push_back(&v);
		markChanged();
		clearHistory(true);

		return g.nodes()-1;
	}

	bool edgeEnabled(int edgeID) const {
		return transitions[edgeID];
	}
	bool isEdge(int edgeID) const {
		return g.isEdge(edgeID);
	}
	bool hasEdge(int edgeID) const {
		return isEdge(edgeID);
	}
	//Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
	int addEdge(int from, int to, int nid = -1) { //, int weight=1
		int id = g.addEdge(from,to,nid);

		modifications++;
		additions = modifications;
		markChanged();

		if(transitions.size()<=id){
			transitions.resize(id+1);
		}

		return id;
	}

	int nEdgeIDs() {
		return g.nEdgeIDs();
	}
	inline int states() const {
		return g.nodes();
	}
	inline int nodes() const {
		return g.nodes();
	}
	inline int edges() const {
		return g.edges();
	}

	bool hasEdge(int from, int to)const{

		return g.hasEdge(from,to);
	}
	int getEdge(int from,int to)const{
		return g.getEdge(from,to);
	}
	inline int nIncident(int node, bool undirected = false) {
		return g.nIncident(node,undirected);
	}

	inline int nIncidentEnabled(int node) {
		assert(node >= 0);
		assert(node < nodes());
		int count = 0;
		Edge e;
		for (auto & e : g.getIncidentDirectedEdges(node)) {
			if (transitions[e.id]) {
				count++;
			}
		}
		return count;
	}

	inline int nDirectedEdges(int node, bool incoming) {
		return g.nDirectedEdges(node,incoming);
	}
	inline DynamicGraph<int>::Edge & directedEdges(int node, int i, bool is_incoming) {
		return g.directedEdges(node,i,is_incoming);
	}

	inline int nIncoming(int node, bool undirected = false) {
		return g.nIncoming(node,undirected);
	}

	inline DynamicGraph<int>::Edge & incident(int node, int i, bool undirected = false) {
		return g.incident(node,i,undirected);
	}
	inline DynamicGraph<int>::Edge & incoming(int node, int i, bool undirected = false) {
		return g.incoming(node,i,undirected);
	}

	DynamicGraph<int>::FullEdge & getEdge(int id)  {
		return g.getEdge(id);
	}

	int getCurrentHistory() {
		return modifications;
	}

	void clearHistory(bool forceClear = false) {
		//long expect=std::max(1000,historyClearInterval*edges());
		if (history.size()
				&& (forceClear
						|| (history.size()
								> (adaptive_history_clear ?
										std::max(1000L, historyClearInterval * edges()) : historyClearInterval)))) {//){
			history.clear();
			historyclears++;

		}
		g.clearHistory();
	}
	//force a new modification
	void invalidate() {
		modifications++;
		additions = modifications;
		modifications++;
		deletions = modifications;
		is_changed = true;

	}

	void markChanged() {
		is_changed = true;

	}
	bool changed() {
		return is_changed;
	}

	void clearChanged() {
		is_changed = false;
		g.clearChanged();
	}

	// Turns a state label into a human readable string
	std::string statelabelToString(int state) {
		std::string s = "";
		for (int i = 0; i < statelabel[state]->size(); i++) {
			if (statelabel[state]->operator [](i) == true) {
				s += "1";
			} else if (statelabel[state]->operator [](i) == false) {
				s += "0";
			} else {
				s += "E";
			}
		}
		return s;
	}

public:
	void draw(int source=-1, int dest=-1, bool forceprint=false){ // forceprint will make it print even if opt_verb=0
		if(opt_verb>1 || forceprint){
			printf("digraph{\n");

			for (int i = 0; i < g.nodes(); i++) {
				std::string s = "node [label=\"" + std::to_string(i) + std::string(": {") + statelabelToString(i) + "}\"] ";
				if(i == dest){
					std::cout << s << "[shape=doublecircle] " << i << ";\n";
				} else {
					std::cout << s << "[shape=circle] " << i << ";\n";
				}
			}


			if(source>=0){
				printf("node [label=\"start\",shape=plaintext] start;\n");
				printf("start->%d\n",source);
			}
			for(int i = 0;i<transitions.size();i++){
				if (transitions[i]){
					printf("%d->%d\n", g.getEdge(i).from,g.getEdge(i).to);
				}
			}

			printf("}\n");
		}
	}


	void clear(){
		g.clear();
		is_changed=true;
		transitions.clear();


		invalidate();
		clearHistory(true);
	}
	void copyTo(DynamicKripke & to){
		to.clear();
		g.copyTo(to.g);
		// TODO Does this actually make a copy? This suggests so: http://www.cplusplus.com/reference/vector/vector/operator=/
		std::vector<bool> transcopy = transitions;
		to.transitions = transcopy;
		std::vector<Bitset*> statelabelcopy = statelabel;
		to.statelabel = statelabelcopy;
	}

};

}
;



#endif /* DynamicKripke_H_ */
