#The MIT License (MIT)
#
#Copyright (c) 2014, Sam Bayless
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
#associated documentation files (the "Software"), to deal in the Software without restriction,
#including without limitation the rights to use, copy, modify, merge, publish, distribute,
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or
#substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
#NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
#OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import monosat.monosat_c
from monosat.logic import *

from monosat.manager import Manager
import sys
debug=False   

#Collects a set of kripkes to encode together into a formula
class CTLManager(metaclass=Manager):
    
 
    
    def  __init__(self):
        self.kripkes=[]        
    
    def clear(self):
        self.kripkes=[]
    
    def addKripke(self, g):
        self.kripkes.append(g)
    
    def getKripke(self,gid):
        return self.kripkes[gid]
    


class KripkeStructure():
        
    def __init__(self,nStates, nProperties):
        self._monosat = monosat.monosat_c.Monosat()
        manager = CTLManager()
        manager.addKripke(self)
        self.kripke = self._monosat.newKripkeStructure(nStates, nProperties)
        
        self.id = len(manager.kripkes)

        self.nodes=0
        self.numedges=0
        self.names=dict()
        self.out_edges=[]   
        self.in_edges=[]     
        self.alledges=[] 
  
        self.all_undirectededges = []
        
        self.edgemap=dict()
        for n in range(nStates):
            self._addState()
    
    def assertCTL(self,formula, start_state=None):
        if start_state is None:
            start_state = 0            
        self._monosat.assertKripkeFormula(self.kripke,start_state,formula)
        
    
    def _addState(self, name=None):  
        n = self.nodes
        #n= self._monosat.newKripke_State(self.kripke)        
        self.nodes=n+1
        self.out_edges.append([]);
        self.in_edges.append([]);        
        
        if name is None:
            name = str(n)
        self.names[n] = str(name)

        return n
    
    def getSymbol(self,node):
        return self.names[node]
 
    def getProperty(self, state, property):
        p= self._monosat.getKripkePropertyLit(self.kripke,state,property)
        return Var(p)
    
    
    def getTransition(self,f,t):
        for (v,w,var,weight) in self.out_edges[f]:
            if(w==t):
                return var;
           
        return None; 
    
    def getAllTransitions(self, undirected=False):
        if undirected:
            return self.all_undirectededges;
        else:
            return self.alledges;

    
    def hasTransition(self,f,t):
        for (v,w,var,w) in self.out_edges[f]:
           if(w==t):
               return True;
           
        return False; 

    def addTransition(self,v,w):
        while(v>=self.numStates() or w>=self.numStates()):
            self.addState()

        var = Var(self._monosat.newKripke_Transition(self.kripke,v,w))
        e=(v,w,var)
        self.alledges.append(e)
        self.numedges=self.numedges+1
        self.out_edges[v].append(e)
        self.in_edges[w].append(e)
        self.edgemap[e[2].getLit()] =e
        return e[2]
    

    
    def numStates(self):
        return self.nodes
    
    def getStates(self):
        return range(self.nodes)
    
    def getTransitionFromVar(self,var):
        return self.edgemap[var.getLit()]
    
    def getTransitions(self,node=-1, undirected=False):
        if(node>=0):
            for edge in  self.out_edges[node]:
                yield edge
            if undirected:
                for edge in  self.in_edges[node]:
                    yield edge                
        else:
            for node in self.out_edges:
                for edge in node:
                    yield edge
                        
    def getIncomingTransitions(self,node=-1):
        if(node>=0):
            for edge in  self.in_edges[node]:
                yield edge
        else:
            for node in self.in_edges:
                for edge in node:
                    yield edge

    def getTransitionVars(self,node=-1):
        if(node>=0):
            for edge in  self.out_edges[node]:
                yield edge[2]
        else:
            for node in self.out_edges:
                for edge in node:
                    yield edge[2]    

 
        
    def draw(self):
        
        print("digraph{")

        for n in range(self.nodes):
            print("n%d"%( n))
        
        for (v,w,var) in self.getAllTransitions():

            print("""n%d->n%d [label="%d"]"""%(v, w, var.getVar()))
            

        print("}") 
        