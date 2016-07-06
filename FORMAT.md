# MonoSAT File Format

(see below for CTL syntax)

MonoSAT is a SAT Modulo Theory solver for Boolean Monotonic Theories, featuring support for a wide set of graph properties, as well as a number of other theories, including some limited geometric properties (currently, for convex hulls), and finite state machine synthesis (still experimental).

This file documents MonoSAT's __*.GNF__ file format, which is a superset of the common [DIMACS CNF format](http://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps) (specifically, a superset of DIMACS as parsed by Minisat).

Just like a DIMACS CNF, a GNF instance begins with a `p cnf <#Vars> <#Clauses>` header line, followed by a number of clauses in 0-terminated format, one per line:

```
p cnf 4 3
1 3 -4 0
4 0
2 -3 0
```

In addition to clauses, a GNF instance may declare any number of graphs. A graph is declared with a line `digraph <weight type> <# nodes> <# edges> <GraphID>`, followed by at most "# edges" edge declarations. GraphID is a non-negative integer that must be unique for each graph; <weight type> may be one of 'int', 'float', 'rational', or may be ommited, in which case 'int' is assumed. After the graph is declared, edges may be declared for that graph on subsequent lines, as `edge <GraphID> <from> <to> <CNF Variable> [weight]`, where 'weight' is optional. If weight is specified, its type must match the graph's declared weight type. Rational values are read using the GMP [mpq_set_str] function. If weight is ommited, the edge has unit weight.

```
p cnf 4 3
1 3 -4 0
4 0
2 -3 0
digraph int 3 4 0
edge 0 0 1 1
edge 0 1 0 2
edge 0 1 2 3 
edge 0 0 2 4 4
```

This creates an integer-weighted graph (with ID '0'), with 3 nodes (numbered 0..2) and 4 directed edges, three of them unit weight, and the last with edge weight 4. Each edge is paired with a unique varaible from the CNF; any variables in the CNF may be used (that is, any variable <= #Vars in the CNF header). The first edge connects node 0 to node 1, IFF variable 1 is assigned true. The second edge connects node 1 to node 0, IFF varible 2 is assigned true. The third edge connects node 1 to node 2, IFF variable 3 is assigned true. The fourth edge connects node 0 to node 2, IFF variable 4 is assigned true. 

Two important restrictions to note are that 
  - edge variables are variables, *not* literals (that is, they must be positive integers)
  - every edge must be assigned a variable 
  - no two edges (or any other theory element from *any* graph in the instance) may share the same variable.

These restrictions reflect implementation choices in MonoSAT, and may be lifted in the future. 

You can specify constraints on the graph using any of the graph properties that MonoSAT supports.  
* Reachability Properties (*e.g.,* node `a` must reach node `b`, IFF `variable` is true): 
    *  `reach <GraphID> <a> <b> <variable>`
* Shortest Path Properties (*e.g.,* shortest path from node `a` to node `b` must have length <= `distance`, IFF `variable` is true):
    * Unweighted:  `distance_leq <GraphID> <a> <b> <variable> <distance>`, with `distance` a non-negative integer. All edges are treated as having unit weight, regardless of their specified weight.
    * Weighted:  `weighted_distance_leq <GraphID> <a> <b> <variable> <distance>`, where `distance `is parsed according to the graph's weight type.
* Maximum Flow Properties (*e.g.,* the maximum flow from node `s` to node `t` must be >= `flow`, IFF `variable` is true):
    *  `maximum_flow_geq <GraphID> <s> <t> <variable> <flow>`
* Minimum Weight Spanning Tree Properties (*e.g.,* the minimum weight spanning tree of graph `GraphID` must be <= `mstweight`, IFF `variable` is true)
    *  `mst_weight_leq <GraphID> <variable> <mstweight>`
* Acyclicity Properties (*e.g.,* `variable` is true IFF graph `GraphID` is a directed acyclic graph (or, for undirected acyclicity, is a forest)
    *  `acyclic <GraphID> <variable>`
    *  `forest <GraphID> <variable>`
    
Variants of these properties  with strict comparison operators are also supported, *e.g.,* `distance_lt`,  `maximum_flow_gt`, `mst_weight_lt`, .

Each of these graph properties is associated with a (unique) Boolean variable in the CNF, which must be true if and only if that graph property holds. This means that by asserting those variables to be true or false in the CNF, one can force any of these graph properties to be either true or false; but one can also trigger complex conditions in the CNF based on the truth-value of the graph property.

Graph reachability properties are specified as ```reach <GraphID> <from> <to> <variable>```. The variable (as with edge variables, a unique positive integer) must be assigned true if from reaches to, and must be false otherwise. For example, to specify that the there must exist a path from node 0 to node 2, we can create a reach property attached to variable 5, and assert that variable 5 must be true in the CNF:

```
p cnf 5 3
1 3 -4 0
4 0
2 -3 0
5 0
digraph int 3 4 0
edge 0 0 1 1
edge 0 1 0 2
edge 0 1 2 3 
edge 0 0 2 4 4
reach 0 0 2 5
```

By instead adding the unit clause '-5 0', you could have forced the graph to *not* have a path from node 0 to node 2. (This would be UNSAT in this example).

Finally, you can combine multiple graph properties, and combine them togehter with arbitrary Boolean logic. For example, you could assert that assert that *either* node 2 must not be reachable from node 0, *or* the weighted shortest path from node 0 to node 2 must be <= 3:

```
p cnf 5 3
5 6 0
-5 -6 0
digraph int 3 4 0
edge 0 0 1 1
edge 0 1 0 2
edge 0 1 2 3 
edge 0 0 2 4 4
reach 0 0 2 5
weighted_distance_leq 0 0 2 6 3 
```

These are the graph properties that are currently well-supported by MonoSAT; many other useful graph properties are Boolean monotonic with respect to the edges in a graph, and could be supported in the future. Interesting possibilities include planarity detection, connected components, global minimum cuts, and many variatons of network flow properties. 



*CTL*


The format for synthesis of Kripke structures by CTL formulas is described in the following. We use the terms "node" and "state" interchangeably, equally "edge" and "transition". In order to synthesize a CTL formula, a "dynamic" Kripke structure has to be created, which serves as a blank piece for the resulting structure. The dynamic Kripke structure can turn off/on transitions and remove/add atomic propositions in states. However, since this is a bounded approach, the number of states, and the number of atomic propositions of the formula has to be fixed in advance.


There are a number of ways to define Kripke structures and CTL formulas, ranging from the most powerful with explicit construction of variables for each atomic proposition and transition, to simplified versions, where a Kripke structure's variables are created automatically given a certain number of states.

Warning: the input is generally NOT checked for correctness/completeness.


== INPUT METHOD 0: The Performant ==
     NOTE: This is the method referred to in the paper as "MonoSAT-structural"

The most performant option for CTL synthesis with an application to synchronization. You specify a number of processes, a number of "wildcard" states (effectively the upper bound of the problem) and a CTL formula which should hold. The input method enforces that each global state specifies the local states for each process and that only one process' local state changes on each local transition.


This method uses MonoSAT's python API, though more out of convenience and not due to performance benefits (MonoSAT is written in C++).
To use it, you first have to compile and install MonoSAT's python API (see the README file for instructions).

You then create a formula file, with the CTL formula as specified below. You then create a script file with the following contents:

```
#       script                         formula file      #processes, #wildcards, enforceLocalStructure, usePBConstraints
python3 path/to/perprocessAutomated.py path/to/CTL-input 2           1           True                   True
```

#processes specifies (unsurprisingly) the number of processes, and #wildcards the number of "wildcard" states, whose local states and transitions are free to assign by the solver. Leave enforceLocalStructure, and usePBConstraints set to "True".

For the formula file, simple use the same format as described in INPUT METHOD 1 but be advised that everything except for the formula will be ignored.


== INPUT METHOD 1: The Non-Invasive Process Encoding ==
     NOTE: This approach offers little benefits towards the one described above ("The Performant"), while having similar restrictions. It is generally slower, so there is no good reason to use it, unless your instances are small enough and you prefer the touch of a pure text input over installing pythons APIs.

For the synthesis of Kripke structures in parallel processes, the following provides some simple format. It creates a certain number of processes, each with its own local transition system of fixed (or rather: bound) size (<States/Process>).
You might want to enforce the following, either via CNF or CTL: In the global transition system that's being synthesized, each transition can only change one process' local state (called process-state from here on). Local states are identified by a <States/Process> atomic propositions per process, which act as indicator variables. Each global state makes for each process exactly one of the atomic proposition belonging to this process true, i.e. its atomic propositions determine the local states of all processes.

```
c     <#Vars> <#Clauses>
p cnf 1000     1
1000 0
c kctlsinglestate <KripkeID> <#Nodes> <#Processes> <States/Process> <selfloops> <ctlvar> <CTL formula>"
kctlsinglestate   0          9        2            3                0           1000     (0 AND EF 4)
AND (AG 0 OR EF 1)
AND EG EX 3
```

- <KripkeID> is always set to 0
- <#Nodes> is the number of states in the resulting structure (can be thought of as an upper bound, since superflous states can be made unreachable, and unreachable states will be automatically removed from the result)
- <#Processes> is the number of parallel processes. They DO NOT have to be identical processes.
- <States/Process> is the maximum number of internal states per process. Note that in total <States/Process>*<#Processes> states will be created, so setting this too high will make the problem much harder to solve.
- <selfloops> set to 1 if you want loops from states to themselves
- <ctlvar> set to a high enough number so that it does not interfere with the variable naming scheme (1000, 10000, ...). This variable is tied to "CTL formula is being satisfied by the Kripke structure".
- <CTL formula> the CTL formula in infix notation. A number in the formula refers to an atomic proposition, which start at 0 and go up to (<#Processes>*<States/Process> - 1). An atomic proposition in a global state corresponds to a certain process being in a certain process-state. The atomic propositions are numbered in the following way: "process x (range 0 to <#Processes>-1) is in process-state y (range 0 to <States/Process> - 1)" is the atomic proposition x * <States/Process> + y. For instance, the second process (process 1) is in the third state (state 2) in a global state iff the atomic proposition 5 is true in this global state.
    Concerning the syntax of the CTL formula: All binary operators must be bracketed in the form (phi1 op phi2). To make things more readable, you can use a newline and start it with "AND", omitting the brakets for the "AND". Note that to enforce the CTL formula on the Kripke structure, you have to add its <ctlvar> as a unit clause to the set of CNF clauses.


For the special case that you have three states per process, you might want to use the following shortcuts, which were created due to the Mutex problem (in which you have a Non-critical section, a Try section and a Critical section): NCS1, TRY1, CS1, NCS2, TRY2, CS2, ... are shortcuts for the APs 0, 1, 2, 3, 4, 5, ... (only to a certain constant). So "(NCS1 AND NCS2) -> EX TRY1" means that if both processes are in the first initial process-state, the first process can transition to the second process-state.

For a realisitc scenario of how this is being used, see the specification of a 2-process Mutex (benchmarks/generic-ctl/mutex/2-mutex-Original/monosat/input).



== INPUT METHOD 2: The Pure ==
     NOTE: This is the method referred to in the paper as "MonoSAT-generic"
kctlsimp is a more powerful input format than kctlsinglestate, in which the instance is oblivious to processes and process-states. You can thus synthesize based on any CTL formula, not knowing about processes, and simply require an upper bound of states.

In contrast to the "kripke" input format specified below, state properties (nodeap) and transitions (edges) are created automatically. It thus really is quite easy to understand and the closest to pure CTL synthesis without looking at any application-specific optimizations and simplifications in the input format.

```
c     <#Vars> <#Clauses>
p cnf 1000     1
1000 0
c kctlsimp <KripkeID> <#Nodes> <#APs> <selfloops> <ctlvar> <CTL formula>"
kctlsimp   0          13       9      0           1000     (0 AND EF 4)
AND (AG 0 OR EF 1)
AND EG EX 3
```


== INPUT METHOD 3: The Powerful ==
Create state properties and transitions manually. If you have some kind of insight on how the resulting structure should look like, you might be able to use this approach more efficiently than a naive creating of N^2 transitions for N states. (The more simplified input formats simply do: "create N states, have M atomic propositions, thus create N^2 transitions and N*M state properties".)

c kripke <#Nodes> <#Edges> <#APs> <KripkeID>
kripke   3        4        2      0            -- means three states, there are at maximum four transitions in the resulting structure (which the solver can turn off/on), there are at maximum two atomic propositions in the formula, and this structure is assigned the ID 0 (note that currently, only synthesis of a single structure within one SMT instance is supported, so leave this at 0).

c kedge <KripkeID> <from> <to> <edgevar>
kedge   0          1      2    7               -- kripke structure ID, from, to, transition SAT solver variable. The transition SAT solver variable must begin with 1 and be increased with every following edge. Consider renaming other variables if they conflict. You can reference the transition SAT solver in CNF clauses: a variable is true iff its transition is enabled (i.e. present) in the structure.

c knodeap <KripkeID> <node> <ap> <nodeapvar>
knodeap   0          2      1    5             -- means in kripke structure 0, node 2, the second AP is tied to the CNF (SAT solver) variable 5. You can reference the <nodeapvar> SAT solver variable in CNF clauses: a variable is true iff for the corresponding (state,ap) pair, the AP is part of the state.

c kctl <KripkeID> <initialnode> <ctlvar> <CTL formula>
kctl   0          0             9        EX (4 EW EG 2) -- the first int after "kctl" is the kripke structure ID, then comes the initial state and then the variable bound to the CTL formula (just pick any variable high enough so that it isn't used elsewhere). After those three integers comes the CTL formula in infix notation. A number in the formula refers to an atomic proposition, which start at 0 and go up to (<#APs> - 1). All binary operators must be bracketed in the form (phi1 op phi2). Note that to enforce the CTL formula on the Kripke structure, you have to add its <ctlvar> as a unit clause to the set of CNF clauses.

```
c     <#Vars> <#Clauses>
p cnf 10      6
-1 -2 0
-4 -3 0
-6 -5 0
-9 -8 0
-6 0
10 0
c kripke <#Nodes> <#Edges> <#APs> <KripkeID>
kripke   3        6        2      0
c kedge <KripkeID> <from> <to> <edgevar>
kedge   0          0      1    1
kedge   0          1      0    2
kedge   0          0      2    3
kedge   0          1      2    4
kedge   0          2      2    5
kedge   0          2      0    6
c knodeap <KripkeID> <node> <ap> <nodeapvar>
knodeap   0          0      0    7
knodeap   0          1      0    8
knodeap   0          2      0    9
c kctl <KripkeID> <initialnode> <ctlvar> <CTL formula>
kctl   0          0             10     EG 0
```





[mpq_set_str]:https://gmplib.org/manual/Initializing-Rationals.html#Initializing-Rationals
