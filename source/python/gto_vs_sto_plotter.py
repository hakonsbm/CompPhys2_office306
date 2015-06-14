import numpy as np
import matplotlib.pyplot as plt

def GTO(x, alpha=1., c=0.3):
	return (2.*alpha/np.pi)**(3./4.)*c*np.exp(-alpha*x*x)

def STO(x, alpha=1.):
	ri = np.sqrt(x*x)
	return np.exp(-alpha*ri)
	


x = np.linspace(-10, 10,501)


plt.figure(figsize=(15,6))
plt.subplot(1, 2, 1)
plt.plot(x, STO(x))
plt.title("STO")
plt.axis([-5,5,0,1])

plt.subplot(1, 2, 2)
plt.plot(x, GTO(x))
plt.title("GTO")
plt.axis([-5,5,0,max(GTO(x))])

plt.savefig("../../Report/figures/GTO_vs_STO_plot.png")
#plt.show()


H_gto_1 = 0.4579*GTO(x, 13.62670, 0.175230)
H_gto_2 = 0.4579*GTO(x, 1.999350, 0.893483)
H_gto_3 = 0.6573*GTO(x, 0.382993, 1.000000)
H_sto = STO(x, 1.843)

plt.figure(figsize=(15,6))
plt.subplot(1,2,1)
plt.plot(x, H_gto_1, "g", label="primitive")
#plt.hold("on")
plt.plot(x, H_gto_2, "g", label="primitive")
plt.plot(x, H_gto_3, "g", label="primitive")
plt.plot(x, H_gto_1 + H_gto_2 + H_gto_3, "r", label="contracted")
#plt.plot(x, H_sto, "r", label="STO")
plt.title("Contracted primitives")
plt.legend()
plt.axis([-5,5,0,max(H_gto_1 + H_gto_2 + H_gto_3)])

plt.subplot(1,2,2)
plt.plot(x, H_gto_1 + H_gto_2 + H_gto_3, "r", label="contracted")
plt.plot(x, H_sto, "b", label="STO")

plt.axis([-5,5,0,max([max(H_sto), max(H_gto_1 + H_gto_2 + H_gto_3)])])
plt.title("STO and contracted GTOs")
plt.legend()


plt.savefig("../../Report/figures/Primitives_vs_STO_plot.png")
plt.show()