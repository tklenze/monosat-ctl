import sys
from monosat import *
################## RUNNING THIS SCRIPT
# 
# Just configure it as below and let it run it without arguments.
# Append the CTL formula you wish to use as specification.
# You may run the scrip with the argument "--stats", which will output useful information about the generated instance
# 
################## CONFIGURATION
if(len(sys.argv)<1):
    print("Please specify ctl formula to encode\n")
    sys.exit(1)
input_file = sys.argv[1]


#Monosat().init("-use-symmetry-reduction=3  -no-ctl-single-state-per-process ")
################## SCRIPT
displayStatistics = False
if len(sys.argv)>2:
    if sys.argv[2] == "--stats":
        displayStatistics = True

APs = 6
states = 7#
kripkeID = 0 # Kripke ID

#varnum = states * (states - 1) + states * APs + 1

k= KripkeStructure(states,APs) #Kripke structure with 3 states, 2 properties


for x in range(0, states):
   edgeVar[x] = dict()
   edges =[]
   for y in range(0, states): 
       if x!= y:
           #print "Adding " + str(x) + " -> " + str(y) + " as var " + str(z)
           e = k.addTransition(x, y)
           edges.append(e)

   AssertClause(edges) #every state must have a next state. 




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
        continue  #this is a dimacs clause, skip it for now (or we could read it in as a clause)

    if line.startswith("kctlsinglestate"):
        words = line.split()
        line = " ".join(words[7:])     
    elif line.startswith("kctlsimp"):
        words = line.split()
        line = " ".join(words[6:])                    
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
