from monosat import *
k= KripkeStructure(2,2) #Kripke structure with 3 states, 2 properties

e1 = k.addTransition(0,1) 
e2 = k.addTransition(1,0)
#e3 = k.addTransition(0,0)
Assert(Not(And(e1,e2)))

k.assertCTL("((0 AU 1) AND NOT 1)")

result = Solve()
print(result)