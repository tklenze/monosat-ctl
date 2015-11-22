states = 6
APs = 6
selfloops = 0 # set to 0 or 1
z = 0 # Kripke ID

varnum = states * (states - 1 + selfloops) + states * APs + 1

print "c     <#Vars> <#Clauses>"
print "p cnf %d      1" % (varnum)
print "%d 0" % (varnum)
print "c kripke <#Nodes> <#Edges> <#APs> <KripkeID>"
print "kripke   %d        %d       %d      0" % (states, (states*(states-1)), APs)
print "c kedge <KripkeID> <from> <to> <edgevar>"
for x in range(0, states):
    for y in range(0, states):
        if x!= y or selfloops == 1:
            z += 1
            print "kedge   0          %d      %d    %d" % (x,y,z)

print "c knodeap <KripkeID> <node> <ap> <nodeapvar>"
            
for x in range(0, states):
    for y in range(0,APs):
        z += 1
        print "knodeap   0          %d      %d    %d" % (x,y,z)

print "c kctl <KripkeID> <initialnode> <ctlvar> <CTL formula>"
print("kctl   0          0             %d       " % (varnum)),
#with open("formula", 'r') as fin:
#    print fin.read()
