import sympy as sp

sp.init_printing()

alpha, beta, Z = sp.symbols('alpha beta Z', positive = True)
x_i, y_i, z_i= sp.symbols('x_i y_i z_i', positive = True)
x_ij, y_ij, z_ij= sp.symbols('x_ij y_ij z_ij', positive = True)

r_i = sp.sqrt(x_i*x_i + y_i*y_i + z_i*z_i)


psi1S = sp.exp(-alpha *r_i)

psi2S = (1-alpha*r_i/2)*sp.exp(-alpha *r_i /2 )

psi2P = alpha*x_i*sp.exp(-alpha*r_i/2)

R_i = sp.symbols('r_i')


def printDerivatives(psi):
	print "Calculating the derivatives for the wavefunction "
	print "psi : ", psi.subs(r_i,R_i)
	print "d/dx : ", (sp.diff(psi,x_i) ).subs(r_i,R_i).factor()#+ sp.diff(psi,y_i) + sp.diff(psi, z_i))

						
	print "d2/dx2 : " , ((sp.diff(psi,x_i,2) + sp.diff(psi,y_i,2) + sp.diff(psi,z_i,2))
		.factor().subs(r_i,R_i).collect(alpha**2).subs(r_i**2,R_i**2))
	print 
	print "Now in latex format "
	print "psi : ", sp.printing.latex(psi.subs(r_i,R_i))
	print "d/dx : ", sp.printing.latex((sp.diff(psi,x_i)  ).subs(r_i,R_i).factor()  + (sp.diff(psi,y_i)  ).subs(r_i,R_i).factor()  + (sp.diff(psi,z_i)  ).subs(r_i,R_i).factor())
					
	print "d2/dx2 : " , sp.printing.latex((sp.diff(psi,x_i,2) + sp.diff(psi,y_i,2) + sp.diff(psi,z_i,2))
		.factor().subs(r_i,R_i).collect(alpha**2).subs(r_i**2,R_i**2))
	print 
	print "Now in c format "
	print "psi : ", sp.printing.ccode(psi.subs(r_i,R_i))
	print "d/dx : ", sp.printing.ccode((sp.diff(psi,x_i)).subs(r_i,R_i).factor()) # + sp.diff(psi,y_i) + sp.diff(psi, z_i)).subs(r_i,R_i).factor())
	print "d2/dx2 : " , sp.printing.ccode((sp.diff(psi,x_i,2) + sp.diff(psi,y_i,2) + sp.diff(psi,z_i,2))
		.factor().subs(r_i,R_i).collect(alpha**2).subs(r_i**2,R_i**2))
	print 
	print "#############################################################################3"
	print

# printDerivatives(psi1S)
# printDerivatives(psi2S)
printDerivatives(psi2P)

#Calculating some derivatives for the correlated part of the wavefunction

def printCorrelationDerivatives(psi):
	print "Calculating the correlation derivatives for the wavefunction "
	print "psi : ", psi.subs(r_ij,R_ij)
	print "d/drij : ", ( sp.diff(psi,x_ij) + sp.diff(psi, y_ij) + sp.diff(psi, z_ij) ).subs(r_ij,R_ij).factor()
	print "d/drij : ", (( sp.diff(psi,x_ij, 2) + sp.diff(psi, y_ij, 2) + sp.diff(psi, z_ij, 2) ).subs(r_ij,R_ij).factor().subs(r_ij**2,R_ij**2))


r_ij = sp.sqrt( (x_ij)**2 +  (y_ij ) **2 + (z_ij )**2 )
psi_C = sp.exp(r_ij / (2*(1+ r_ij*beta)))

R_ij = sp.symbols('r_ij')

# printCorrelationDerivatives(psi_C)

test = 6*beta**2*r_ij**4 - 6*beta**2*r_ij**2*x_ij**2 - 6*beta**2*r_ij**2*y_ij**2 - 6*beta**2*r_ij**2*z_ij**2 + 12*beta*r_ij**3

print "Testing"
print sp.printing.latex(test.subs(r_ij,R_ij).subs(r_ij**2,R_ij**2).subs(r_ij,R_ij).factor(r_ij**2))



