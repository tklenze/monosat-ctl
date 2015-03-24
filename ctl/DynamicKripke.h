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
using namespace dgl;
namespace Monosat {



class DynamicKripke{
	DynamicGraph g;
	//std::vector<Bitset> edge_status;
	int id;
	bool has_epsilon=true;
	bool is_changed = true;
	bool is_generator=true;
	bool is_acceptor=true;
	bool is_linear=true;
public:
	vec<Bitset> transitions;

	// This represents the labels of states, which are sets of atomic propositions
	// Access via statelabel[stateID][AP]
	vector<vector<int>> statelabel;

	bool adaptive_history_clear = false;
	long historyClearInterval = 1000;
	int modifications=0;
	int additions=0;
	int deletions=0;
	int in_alphabet =1;
	int out_alphabet=1;
	long historyclears=0;
	struct EdgeChange {
		bool addition;

		int id;
		int input;
		int output;
		int mod;
		int prev_mod;

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
	void setEmovesEnabled(bool enabled){
		has_epsilon=enabled;
	}

	bool emovesEnabled()const{
		return has_epsilon;
	}

/*
	bool emove(int edgeID)const{
		return emovesEnabled() && transitions[edgeID][0];
	}
*/

	int inAlphabet()const{
		return in_alphabet;
	}
	int outAlphabet()const{
		return out_alphabet;
	}
	void addInCharacter(){
		in_alphabet++;
	}
	void addOutCharacter(){
		out_alphabet++;
	}
	bool transitionEnabled(int edgeID, int input, int output)const{
		assert(input<inAlphabet());
		assert(output<outAlphabet());
		if(output<0){
			for(int i = 0;i<out_alphabet;i++){
				if(transitionEnabled(edgeID,input,i))
					return true;
			}
			return false;
		}else if (input<0){
			for(int i = 0;i<in_alphabet;i++){
				if(transitionEnabled(edgeID,i,output))
					return true;
			}
			return false;
		}else{
			int pos = input +output*inAlphabet();
			return transitions[edgeID][pos];
		}
	}



/*
	void setStateLabel(int state, vec<bool> label) {
		statelabel[state]. label;
	}
	*/

	/*
	void setAP(int apcount) {
		statelabel[state] = label;
	}*/

	// Is ap in the label of state?
	// 1: Yes, 0: No, -1: AP not present
	int stateAPAssignment(int state, int ap) {
		if (statelabel[state].size() > ap) {
			if (statelabel[state][ap])
				return 1;
			else
				return 0;
		}
		return -1;
	}












	int addTransition(int from, int to,int edgeID, int input,int output, bool defaultEnabled=true){
		assert(input<inAlphabet());
		assert(output<outAlphabet());
		while(from>=g.nodes() || to>=g.nodes())
			g.addNode();
		if(edgeID==-1){
			edgeID = g.addEdge(from, to, edgeID);
		}
		transitions.growTo(edgeID+1);
		transitions[edgeID].growTo(inAlphabet()*outAlphabet());
		int pos = input +output*inAlphabet();
		if(defaultEnabled)
			transitions[edgeID].set(pos);
		return edgeID;
	}

	void enableTransition(int edgeID, int input,int output) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		int pos = input +output*inAlphabet();
		if (!transitions[edgeID][pos]) {
			transitions[edgeID].set(pos);
			//edge_status.setStatus(id,true);
			modifications++;
			additions = modifications;
			history.push_back( { true, edgeID,input,output, modifications, additions });
		}
	}
	void disableTransition(int edgeID, int input,int output) {
		assert(edgeID >= 0);
		assert(edgeID < g.edges());
		assert(isEdge(edgeID));
		int pos = input +output*inAlphabet();
		if (transitions[edgeID][pos]) {
			transitions[edgeID].clear(pos);
			modifications++;
			history.push_back( { false, edgeID,input,output, modifications, deletions });
			deletions = modifications;
		}
	}
	int addEmptyState(){
		vector<int> empty;
		return addNode(empty);
	}
	int addState(vector<int>& v){
		return addNode(v);
	}
	int addNode(vector<int> v) {

		g.addNode();
		modifications++;
		additions = modifications;
		deletions = modifications;
		statelabel.push_back(v);
		markChanged();
		clearHistory(true);

		return g.nodes()-1;
	}

	bool edgeEnabled(int edgeID) const {
		return g.edgeEnabled(edgeID);
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
			transitions.growTo(id+1);
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

	inline int nDirectedEdges(int node, bool incoming) {
		return g.nDirectedEdges(node,incoming);
	}
	inline DynamicGraph::Edge & directedEdges(int node, int i, bool is_incoming) {
		return g.directedEdges(node,i,is_incoming);
	}

	inline int nIncoming(int node, bool undirected = false) {
		return g.nIncoming(node,undirected);
	}

	inline DynamicGraph::Edge & incident(int node, int i, bool undirected = false) {
		return g.incident(node,i,undirected);
	}
	inline DynamicGraph::Edge & incoming(int node, int i, bool undirected = false) {
		return g.incoming(node,i,undirected);
	}

	DynamicGraph::FullEdge getEdge(int id) const {
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
	string statelabelToString(int state) {
		string s = "";
		for (int i = 0; i < statelabel[state].size(); i++) {
			if (statelabel[state][i] == true) {
				s += "1";
			} else if (statelabel[state][i] == false) {
				s += "0";
			} else {
				s += "E";
			}
		}
		return s;
	}

	void draw(int source=-1, int dest=-1){
		printf("digraph{\n");

		for (int i = 0; i < g.nodes(); i++) {
			string s = "node [label=\"" + std::to_string(i) + std::string(": {") + statelabelToString(i) + "}\"] ";
			if(i == dest){
				cout << s << "[shape=doublecircle] " << i << ";\n";
			} else {
				cout << s << i << ";\n";
			}
		}


		if(source>=0){
			printf("node [label=\"start\",shape=plaintext] start;\n");
			printf("start->%d\n",source);
		}
		for(int i = 0;i<transitions.size();i++){
			bool any_enabled=false;
			for(int l= 0;l<transitions[i].size();l++){
				if(transitions[i][l]){
					any_enabled=true;
					break;
				}
			}
			if (any_enabled){
				// Ignore transition labels, we don't need them for kripke structures.
				printf("%d->%d\n", g.getEdge(i).from,g.getEdge(i).to);

				/*
				printf("%d->%d [label=\"", g.getEdge(i).from,g.getEdge(i).to);

				for(int in = 0;in<inAlphabet();in++){
					for(int out = 0;out<outAlphabet();out++){
						int pos = in + inAlphabet()*out;
						if(transitions[i][pos]){
							if(out==0){
								if(in==0){
									printf("{},");
								}else{
									printf("%c:,",'a'+in-1);
								}
							}else{
								if(in==0){
									printf(":%c,",'a'+out-1);
								}else{
									printf("%c:%c,",'a'+in-1,'a'+out-1);
								}
							}

						}
					}
				}

				printf("\"]\n");
				*/

			}
		}


		printf("}\n");

	}


	bool buildPrefixTable(int startState, int finalState, vec<int> & string, vec<Bitset> & table){
		table.growTo(string.size());
		for(int i = 0;i<table.size();i++){
			table[i].clear();
			table[i].growTo(this->states());
		}
		if(string.size()==0)
			return startState==finalState;//this isn't quite correct, because there may be emoves connecting start to final state...
		static vec<int> curStates;
		static vec<int> nextStates;
		nextStates.clear();
		curStates.clear();
		curStates.push(startState);
		int pos = 0;
		table[pos].set(startState);

		//initial emove pass:
		if(emovesEnabled()){
			for(int i = 0;i<curStates.size();i++){
				int s = curStates[i];
				for(int j = 0;j<nIncident(s);j++){
					//now check if the label is active
					int edgeID= incident(s,j).id;
					int to = incident(s,j).node;
					if(!table[pos][to] && transitionEnabled(edgeID,0,-1)){
						table[pos].set(to);
						curStates.push(to);
					}

				}
			}
		}

		for(;pos<string.size();pos++)
		{
			int l = string[pos];
			assert(l>0);
			for(int i = 0;i<curStates.size();i++){
				int s = curStates[i];
				for(int j = 0;j<nIncident(s);j++){
					//now check if the label is active
					int edgeID= incident(s,j).id;
					int to = incident(s,j).node;
					if(!table[pos][to] && transitionEnabled(edgeID,0,-1)){
						table[pos].set(to);
						curStates.push(to);
					}

					if (pos+1<string.size() && !table[pos+1][to] && transitionEnabled(edgeID,l,-1)){
						table[pos+1].set(to);
						nextStates.push(to);
					}
				}
			}

			nextStates.swap(curStates);
			nextStates.clear();
		}

		pos = string.size()-1;
		//final emove pass:
		if(emovesEnabled()){
			for(int i = 0;i<curStates.size();i++){
				int s = curStates[i];
				for(int j = 0;j<nIncident(s);j++){
					//now check if the label is active
					int edgeID= incident(s,j).id;
					int to = incident(s,j).node;
					if(!table[pos][to] && transitionEnabled(edgeID,0,-1)){
						table[pos].set(to);
						curStates.push(to);
					}
				}
			}
		}
		return table[pos][finalState];
	}

	bool buildSuffixTable(int startState, int finalState, vec<int> & string, vec<Bitset> & table){

		table.growTo(string.size()+1);
		for(int i = 0;i<table.size();i++){
			table[i].clear();
			table[i].growTo(this->states()+1);
		}
		if(string.size()==0)
			return startState==finalState;//this isn't quite correct, because there may be emoves connecting start to final state...
		static vec<int> curStates;
		static vec<int> nextStates;
		nextStates.clear();
		curStates.clear();
		curStates.push(finalState);
		int pos = string.size();
		table[pos].set(finalState);

		//initial emove pass:
		if(emovesEnabled()){
			for(int i = 0;i<curStates.size();i++){
				int s = curStates[i];
				for(int j = 0;j<nIncoming(s);j++){
					//now check if the label is active
					int edgeID= incoming(s,j).id;
					int to = incoming(s,j).node;
					if(!table[pos][to] && transitionEnabled(edgeID,0,-1)){
						table[pos].set(to);
						curStates.push(to);
					}

				}
			}
		}

		for(;pos>0;pos--)
		{
			int l = string[pos-1];
			assert(l>0);
			for(int i = 0;i<curStates.size();i++){
				int s = curStates[i];
				for(int j = 0;j<nIncoming(s);j++){
					//now check if the label is active
					int edgeID= incoming(s,j).id;
					int to = incoming(s,j).node;
					if(!table[pos][to] && transitionEnabled(edgeID,0,-1)){
						table[pos].set(to);
						curStates.push(to);
						//status.reaches(str,to,edgeID,0);
					}

					if (pos>0 && !table[pos-1][to] && transitionEnabled(edgeID,l,-1)){
						//status.reaches(str,to,edgeID,l);
						table[pos-1].set(to);
						nextStates.push(to);
					}
				}
			}

			nextStates.swap(curStates);
			nextStates.clear();

		}
		pos = 0;
		//final emove pass:
		if(emovesEnabled()){
			for(int i = 0;i<curStates.size();i++){
				int s = curStates[i];
				for(int j = 0;j<nIncoming(s);j++){
					//now check if the label is active
					int edgeID= incoming(s,j).id;
					int to = incoming(s,j).node;
					if(!table[pos][to] && transitionEnabled(edgeID,0,-1)){
						table[pos].set(to);
						curStates.push(to);
					}
				}
			}
		}

		return table[0][startState];
	}
private:
	vec<int> next;
	vec<int> cur;

	vec<bool> next_seen;
	vec<bool> cur_seen;

	bool generates_path_rec(int s,int final,int emove_count, vec<NFATransition> & path){
			if(s==final){
				return true;
			}
			if (emove_count>=states()){
				return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
			}


			for(int j = 0;j<nIncident(s);j++){
				//now check if the label is active
				int edgeID= incident(s,j).id;
				int to = incident(s,j).node;
				if( transitionEnabled(edgeID,0,0)){


						path.push({edgeID,0,0});
						if(generates_path_rec(to,final,emove_count+1,path)){//str_pos is NOT incremented!

							return true;
						}else{

							path.pop();
						}

				}
				for(int l = 0;l<outAlphabet();l++){
					if (transitionEnabled(edgeID,0,l)){
						bool set_transition=false;


						path.push({edgeID,0,l});
						if(generates_path_rec(to,final,0,path)){//str_pos is incremented

							return true;
						}else{

							path.pop();
						}

					}
				}

			}
			return false;
		}

public:
	bool generates(int source, int final, vec<NFATransition> & path){
		return generates_path_rec(source,final,0,path);
	}


	bool accepts(int source, int final, vec<int> & string){
		return accepts_prefix(source,final,string)==string.size();
	}
	int accepts_prefix(int source, int final, vec<int> & string){
		cur_seen.growTo(states());
		next_seen.growTo(states());
		for(int s:cur){
			assert(cur_seen);
			cur_seen[s]=false;
		}
		cur.clear();
		assert(next.size()==0);
		cur_seen[source]=true;
		cur.push(source);

		int largest_prefix=0;

		//initial emove pass:
		if(emovesEnabled()){
			for(int i = 0;i<cur.size();i++){
				int s = cur[i];
				for(int j = 0;j<nIncident(s);j++){
					//now check if the label is active
					int edgeID= incident(s,j).id;
					int to = incident(s,j).node;
					if(!cur_seen[to] && transitionEnabled(edgeID,0,-1)){
						cur_seen[to]=true;
						cur.push(to);

					}

				}
			}
		}
		if(string.size()){
			for(int k = 0;k<string.size();k++)
			{
				int l = string[k];
				assert(l>0);
				for(int i = 0;i<cur.size();i++){
					int s = cur[i];
					for(int j = 0;j<nIncident(s);j++){
						//now check if the label is active
						int edgeID= incident(s,j).id;
						int to = incident(s,j).node;
						if(!cur_seen[to] && transitionEnabled(edgeID,0,-1)){
							cur_seen[to]=true;
							cur.push(to);
							//status.reaches(str,to,edgeID,0);
						}

						if (!next_seen[to] && transitionEnabled(edgeID,l,-1)){
							//status.reaches(str,to,edgeID,l);
							next_seen[to]=true;
							next.push(to);
						}
					}
				}

				next.swap(cur);
				next_seen.swap(cur_seen);

				for(int s:next){
					assert(next_seen[s]);
					next_seen[s]=false;
				}
				next.clear();

				if(cur.size()){
					largest_prefix=k;
				}
			}

			//final emove pass:
			if(emovesEnabled()){
				for(int i = 0;i<cur.size();i++){
					int s = cur[i];
					for(int j = 0;j<nIncident(s);j++){
						//now check if the label is active
						int edgeID= incident(s,j).id;
						int to = incident(s,j).node;
						if(!cur_seen[to] && transitionEnabled(edgeID,0,-1)){
							cur_seen[to]=true;
							cur.push(to);
						}

					}
				}
			}
		}
		if( cur_seen[final]){
			return string.size();
		}else{
			return largest_prefix;
		}

	}
	void clear(){
		g.clear();
		has_epsilon=true;
		is_changed=true;
		transitions.clear();

		in_alphabet =1;
		out_alphabet=1;


		next.clear();
		cur.clear();
		next_seen.clear();
		cur_seen.clear();

		invalidate();
		clearHistory(true);
	}
	void copyTo(DynamicKripke & to){
		to.clear();
		g.copyTo(to.g);
		to.has_epsilon=has_epsilon;
		to.in_alphabet=in_alphabet;
		to.out_alphabet=out_alphabet;
		for(Bitset & b:transitions){
			to.transitions.push();
			b.copyTo(to.transitions.last());
		}
	}

};

}
;



#endif /* DynamicKripke_H_ */
