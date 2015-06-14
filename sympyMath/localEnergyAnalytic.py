import sympy as sp

sp.init_printing()


#Defining general stuff needed to solve for the local energy

x1, y1, z1= sp.symbols('x1 y1 z1', positive = True)
x2, y2, z2= sp.symbols('x2 y2 z2', positive = True)

alpha, beta, Z = sp.symbols('alpha beta Z', positive = True)

r1 = sp.sqrt(x1*x1 + y1*y1 + z1*z1)
r2 = sp.sqrt(x2*x2 + y2*y2 + z2*z2)

r12 = sp.sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )


phi1 = sp.exp(-alpha*(r1 + r2))
phi2 = sp.exp(-alpha*(r1 + r2))*sp.exp(r12/(2*(1+beta*r12)))

R1 = sp.Symbol('r1') #Creates a symbolic equivalent of r1
R2 = sp.Symbol('r2')
R12 = sp.Symbol('r12')


def calculateLocalEnergy(phi):
	return 	(-sp.diff(phi, x1,2) - sp.diff(phi, y1, 2) - sp.diff(phi, z1,2) +
			( - sp.diff(phi, x2,2) - sp.diff(phi, y2, 2) - sp.diff(phi, z2, 2 )))/(2*phi)
			# - Z/r1 - Z/r2 + 1 /r12




#Calculating the localenergy1
localEnergy1 = calculateLocalEnergy(phi1)

# #This is the part where the expression for the local energy is simplified to be readable for humans
localEnergy1 = localEnergy1.subs(r12,R12).factor().subs(r1, R1).subs(r1**2,R1**2).subs(r2,R2).subs(r2**2,R2**2).simplify().collect(Z).collect(alpha)

print localEnergy1 
print 
print sp.printing.latex(localEnergy1)
# print sp.printing.ccode(localEnergy1)

#I was not able to simplify the expression for the local energy with the more complicated trialfunction

# localEnergy2 = calculateLocalEnergy(phi2)
# localEnergy2 = localEnergy2#Add lots of simplifying here

# print sp.printing.latex(localEnergy2)