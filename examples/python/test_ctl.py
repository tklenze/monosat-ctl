from monosat import *
k= KripkeStructure(3,2) #Kripke structure with 3 states, 2 properties

e1 = k.addTransition(0,1) 
e2 = k.addTransition(1,2)
e3 = k.addTransition(0,2)
Assert(Not(And(e1,e3)))

k.assertCTL("((0 AU 1) AND NOT 1)")

result = Solve()
print(result)