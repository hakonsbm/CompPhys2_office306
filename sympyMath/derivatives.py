import sympy as sp

sp.init_printing()

alpha, beta, Z = sp.symbols('alpha beta Z', positive = True)
xi, yi, zi= sp.symbols('xi yi zi', positive = True)

ri = sp.sqrt(xi*xi + yi*yi + zi*zi)

psi1S = sp.exp(-alpha *ri)

psi2S = (1-alpha*ri/2)*sp.exp(-alpha *ri)

psi2P = alpha*ri*sp.exp(-alpha*ri/2)

Ri = sp.symbols('ri')

def printDerivatives(psi):
	print "Calculating the derivatives for the wavefunction "
	print "psi : ", psi
	print "d/dx : ", (sp.diff(psi,xi) + sp.diff(psi,yi) + sp.diff(psi, zi)).subs(ri,Ri).factor()
	print "d2/dx2 : " , (sp.diff(psi,xi,2) + sp.diff(psi,yi,2) + sp.diff(psi,zi,2)).factor().subs(ri,Ri)
	print 

printDerivatives(psi1S)
printDerivatives(psi2S)
printDerivatives(psi2P)


# print 
# print "Psi2S: ", psi2S
# print "d/dx : ", sp.diff(psi2S,ri)
# print "d2/dx2 : " , sp.diff(psi2S,ri,2)

# print 
# print "Psi2P: ", psi2P
# print "d/dx : ", sp.diff(psi2P,ri)
# print "d2/dx2 : " , sp.diff(psi2P,ri,2)
