#
# testing the lswt algorithm (https://arxiv.org/abs/1402.6069) for a 1d ferromagnetic chain
# @author Tobias Weber <tweber@ill.fr>
# @date 24-oct-2024
# @license GPLv3, see 'LICENSE' file
#

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


# magnetic site
site_S = 1.
site_u = np.array([ 1, 1j, 0 ])
site_v = np.array([ 0,  0, 1 ])


# coupling
term_J    = -1.
term_DMIx = 0.
term_DMIy = 0.
term_DMIz = 0.
term_dist = np.array([ 1, 0, 0 ])


# real interaction matrix
J_real = np.array([
	[     term_J,   term_DMIz,  -term_DMIy ],
	[ -term_DMIz,      term_J,   term_DMIx ],
	[  term_DMIy,  -term_DMIx,      term_J ] ])

#print("\nJ_real =")
#print(J_real)


def get_energies(Qvec):
	# fourier transformed interaction matrix
	phase = -1j * 2.*np.pi * np.dot(term_dist, Qvec)

	J_fourier = J_real*np.exp(phase) + J_real.transpose()*np.exp(-phase)
	J0_fourier = J_real + J_real.transpose()

	#print("\nJ_fourier =")
	#print(J_fourier)
	#print("\nJ0_fourier =")
	#print(J0_fourier)


	# hamilton
	H00 = 0.5 * site_S * np.dot(site_u, np.dot(J_fourier, site_u.conj())) - \
		site_S * np.dot(site_v, np.dot(J0_fourier, site_v))
	HNN = 0.5 * site_S * np.dot(site_u.conj(), np.dot(J_fourier, site_u)) - \
		site_S * np.dot(site_v, np.dot(J0_fourier, site_v))
	H0N = 0.5 * site_S * np.dot(site_u, np.dot(J_fourier, site_u))

	H = np.array([
		[ H00,                    H0N ],
		[ H0N.transpose().conj(), HNN ] ])
	#print("\nH =")
	#print(H)


	C = la.cholesky(H)
	H_trafo = np.dot(C.transpose().conj(), np.dot(np.diag([1, -1]), C))

	#print("\nC =")
	#print(C)
	#print("\nH_trafo =")
	#print(H_trafo)


	Es = np.real(la.eigvals(H_trafo))
	return Es


hs = []
Es = []
for h in np.linspace(-1, 1, 128):
	try:
		Qvec = np.array([ h, 0, 0 ])
		Es.append(get_energies(Qvec))
		hs.append(h)
	except la.LinAlgError:
		pass

plt.plot()
plt.plot(hs, Es)
plt.show()
