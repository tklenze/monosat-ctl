import sys
from monosat import *
################## RUNNING THIS SCRIPT
# 
# Just configure it as below and let it run it without arguments.
# Append the CTL formula you wish to use as specification.
# You may run the scrip with the argument "--stats", which will output useful information about the generated instance
# 
################## CONFIGURATION
input_file = "../sample-inputs/mutex/original-clarke--kctlsimp--7-states-sat"
processes = 3
localstatesPerProcess = 3 # We assume that each process has the same number of local states (though they might behave differently!). If this does not suit you, re-write the script :)
wildcardStates = 7

# We assume that processes loops through their local states, i.e. 0->1, 1->2, ... n-1->0
enforceLocalStructure = False

Monosat().init("-use-symmetry-reduction=0 -only-one-process-moves=0")
################## SCRIPT
displayStatistics = False
if len(sys.argv)>1:
    if sys.argv[1] == "--stats":
        displayStatistics = True

APs = processes * localstatesPerProcess
states = localstatesPerProcess ** processes + wildcardStates + 1 # The +1 is the initial state with the state ID 0
kripkeID = 0 # Kripke ID

#varnum = states * (states - 1) + states * APs + 1

k= KripkeStructure(states,APs) #Kripke structure with 3 states, 2 properties
clausesCount = dict()
clausesCount['Processes local state cannot jump'] = 0
clausesCount['Each process is at least in one local state'] = 0
clausesCount['Each process is at most in one local state'] = 0
clausesCount['Only one process moves at a time, for undetermined states'] = 0
clausesCount['Determined state APs (Unit clauses)'] = 0

# Determine State APs for fully determined states (i.e. not the initial state or the wildcards). Saved in the variable localState[state][processID].
localStateInit = dict()
# Reset counters
for pID in range(0, processes):
    localStateInit[pID] = 0
localState = dict()
for s in range(1, states - wildcardStates): # We skip the initial state (state 0), and the wildcard states
    localState[s] = dict()
    #print "state " + str(s) + " "
    for pID in range(0, processes):
        localState[s][pID] = localStateInit[pID]
        #print str(localState[s][pID]) + " "
    #print "\n"
    # Increase counter
    for pID in range(processes-1, -1, -1):
        if (localStateInit[pID] + 1) != localstatesPerProcess:
            localStateInit[pID] = (localStateInit[pID] + 1) % localstatesPerProcess
            break
        localStateInit[pID] = (localStateInit[pID] + 1) % localstatesPerProcess

edgeVar = dict()
localStateVar = dict()

def initEdgeCounter():
    global edgeCount
    for x in range(0, states):
       edgeVar[x] = dict()
       for y in range(1, states): # skip state 0, the initial state
           if x!= y and oneProcessMoves(x,y):
               #print "Adding " + str(x) + " -> " + str(y) + " as var " + str(z)
               edgeVar[x][y] = k.addTransition(x, y)


def initStateCounter():
    for s in range(0, states):
        localStateVar[s] = dict()
        property = 0
        for pID in range(0, processes):
            localStateVar[s][pID] = dict()
            for localS in range(0, localstatesPerProcess):
                localStateVar[s][pID][localS] = k.getProperty(s, property)
                property+=1


def determined(s):
    if not s in localState:
        return False
    for pID in range(0, processes):
        if not pID in localState[s]:
            return False
    return True

def addEachProcessInExactlyOneLocalStateConstraint():
    for s in range(0, states):
        if not determined(s):
            for pID in range(0, processes):
                oneOfEvery = list()
                for localS in range(0, localstatesPerProcess):
                    oneOfEvery.append(localStateVar[s][pID][localS])
                clausesCount['Each process is at least in one local state'] += 1
                AssertClause(oneOfEvery)
                #clauses.append(oneOfEvery)
                for localS1 in range(0, localstatesPerProcess):
                    for localS2 in range(localS1, localstatesPerProcess):
                        if localS1 != localS2:
                            notBoth = list()
                            notBoth.append(~localStateVar[s][pID][localS1])
                            notBoth.append(~localStateVar[s][pID][localS2])
                            AssertClause(notBoth)
                            clausesCount['Each process is at most in one local state'] += 1

# Add a constraint that the edge between s1 and s2 should only exist if exactly one process changes
# We assume that at least one of the states is not fully determined
def addOneProcessMovesConstraint(s1, s2):
    for p1 in range(0, processes):
        for p2 in range(p1, processes):
            if p1 != p2:
                for localS1 in range(0, localstatesPerProcess):
                    for localS2 in range(0, localstatesPerProcess):
                        #print "addOneProcessMovesConstraint: " + str(s1) + " -> " + str(s2) + " p1: " + str(p1) + "_" + str(localS1) + "  " + str(p2) + "_" + str(localS2)
                        onlyOneMoves = list()
                        onlyOneMoves.append(~localStateVar[s1][p1][localS1])
                        onlyOneMoves.append(~localStateVar[s1][p2][localS2])
                        onlyOneMoves.append(localStateVar[s2][p1][localS1])
                        onlyOneMoves.append(localStateVar[s2][p2][localS2])
                        onlyOneMoves.append(~edgeVar[s1][s2])
                        AssertClause(onlyOneMoves)
                        #clauses.append(onlyOneMoves)
                        clausesCount['Only one process moves at a time, for undetermined states'] += 1
                        
    return True

# Processes' local state can't jump. E.g. you can't go from local state 0 to local state 2, you have to go via local state 1. (for instance, in the mutex, don't go from NCS directly into CS).
def addLocalStateMovesByOneConstraint(s1, s2):
    for p in range(0, processes):
        for l in range(0, localstatesPerProcess):
            for offset in range(2, localstatesPerProcess): # "jump distance"
                # If one of the two states has been determined, and their local state is different from l (and l+offset mod localstatesPerProcess respectively), then we can safely ignore this clause and not add it, since it will be trivially satisfied
                if determined(s1) and localState[s1][p]!= l:
                    break
                if determined(s2) and localState[s2][p]!=(l+offset) % localstatesPerProcess :
                    break
                #print "addLocalStateMovesByOneConstraint: " + str(s1) + " -> " + str(s2) + " p: " + str(p) + " can't move from local state " + str(l) + " to local state " + str((l+offset) % localstatesPerProcess)
                movesByOne = list()
                movesByOne.append(~localStateVar[s1][p][l])
                movesByOne.append(~localStateVar[s2][p][(l+offset) % localstatesPerProcess])
                movesByOne.append(~edgeVar[s1][s2])
                #clauses.append(movesByOne)
                AssertClause(movesByOne)
                clausesCount['Processes local state cannot jump'] += 1
                
        
def oneProcessMovesUndetermined():
    for s1 in range(0, states):
        for s2 in range(1, states): # skip state 0, the initial state
            if s1 != s2 and (not determined(s1) or not determined(s2)):
                addOneProcessMovesConstraint(s1, s2)
                if enforceLocalStructure:
                    addLocalStateMovesByOneConstraint(s1, s2)

# Decides if an edge between s1 and s2 is possible, or if it is forbidden, since not exactly one process moves between the two states. If one of the states contains undetermined APs, we will return True and later add a clause via addOneProcessMovesConstraint
def oneProcessMoves(s1,s2):
    if not determined(s1) or not determined(s2): # If one of the state's APs are undetermined, then always return true, since we can't check
        return True
    numProcessMove = 0
    for pID in range(0, processes):
        if localState[s1][pID] != localState[s2][pID]:
            if enforceLocalStructure and (localState[s1][pID]+1) % localstatesPerProcess != localState[s2][pID]:
                #print "Not adding edge " + str(s1) + " -> " + str(s2) + " p: " + str(pID) + " local state " + str(localState[s1][pID]) + " moves to " + str(localState[s2][pID]) + " which is different from " + str(((localState[s1][pID]+1) % localstatesPerProcess))
                return False
            numProcessMove += 1
    return numProcessMove == 1
    
def addClausesDeterminedState(): # this encodes as unit clauses, the determined states' APs
    for s in range(0, states):
        if determined(s):
            for pID in range(0, processes):
                for localS in range(0, localstatesPerProcess):
                    determinedAP = list()
                    if localS == localState[s][pID]:
                        determinedAP.append(localStateVar[s][pID][localS])
                    else:
                        determinedAP.append(~localStateVar[s][pID][localS])
                    #clauses.append(determinedAP)
                    AssertClause(determinedAP)
                    clausesCount['Determined state APs (Unit clauses)'] += 1

initEdgeCounter()
initStateCounter()

addClausesDeterminedState()
addEachProcessInExactlyOneLocalStateConstraint()
oneProcessMovesUndetermined()

#set CTL formula here
for line in open(input_file,"r"):
    full_line=line
    line = line.strip()
    if len(line) == 0:
        continue
    if line.startswith("c"):
        continue
    
    if line.startswith("p cnf"):
        continue
    if line.startswith("p cnf"):
        continue
    tokens = line.split()
    all_ints=True
    for f in tokens:
        if not f.isdigit():
            all_ints = False
            break
    if all_ints:
        continue  #this is a dimacs clause, skip it

    if line.startswith("kctlsinglestate"):
        words = line.split()
        line = " ".join(words[7:])            
    elif line.startswith("kctl"):
        words = line.split()
        line = " ".join(words[6:])

    if line.startswith("AND"):
        line = line[4:]
        line = line.strip()
    print("asserting CTL formula: " + line)
    k.assertCTL(line)
    

print("Solving...")
result = Solve()
print(result)
print("Done")
