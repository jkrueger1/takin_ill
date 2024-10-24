#
# testing the lswt algorithm (https://arxiv.org/abs/1402.6069) for a 1d ferromagnetic chain
# @author Tobias Weber <tweber@ill.fr>
# @date 24-oct-2024
# @license GPLv3, see 'LICENSE' file
#

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


# magnetic sites
sites = [
	{ "S" : 1., "u" : np.array([ 1, 1j, 0 ]), "v" : np.array([ 0, 0, 1 ]) },
]

# magnetic couplings
couplings = [
	{ "sites" : [ 0, 0 ], "J" : -1., "DMI" : [ 0, 0, 0 ], "dist" : [ 1, 0, 0 ] },
]


# calculate real interaction matrices
for coupling in couplings:
	coupling["J_real"] = np.array([
		[       coupling["J"],   coupling["DMI"][2],  -coupling["DMI"][1] ],
		[ -coupling["DMI"][2],        coupling["J"],   coupling["DMI"][0] ],
		[  coupling["DMI"][1],  -coupling["DMI"][0],      coupling["J"] ] ])

	#print("\nJ_real =\n%s" % coupling["J_real"])


def get_energies(Qvec):
	#print("\n\nQ = %s" % Qvec)

	# fourier transform interaction matrices
	num_sites = len(sites)
	J_fourier = np.zeros((num_sites, num_sites, 3, 3), dtype = complex)
	J0_fourier = np.zeros((num_sites, num_sites, 3, 3), dtype = complex)

	for coupling in couplings:
		dist = np.array(coupling["dist"])
		J_real = coupling["J_real"]
		site1 = coupling["sites"][0]
		site2 = coupling["sites"][1]

		phase = -1j * 2.*np.pi * np.dot(dist, Qvec)

		J_fourier[site1, site2] += J_real * np.exp(phase)
		J_fourier[site2, site1] += J_real.transpose() * np.exp(-phase)
		J0_fourier[site1, site2] += J_real
		J0_fourier[site2, site1] += J_real.transpose()

	#print("\nJ_fourier =\n%s\n\nJ0_fourier =\n%s" % (J_fourier, J0_fourier))


	# hamilton
	H = np.zeros((2*num_sites, 2*num_sites), dtype = complex)

	for i in range(num_sites):
		S_i = sites[i]["S"]
		u_i = sites[i]["u"]
		v_i = sites[i]["v"]

		for j in range(num_sites):
			S_j = sites[j]["S"]
			u_j = sites[j]["u"]
			v_j = sites[j]["v"]
			S = 0.5 * np.sqrt(S_i * S_j)

			H[i, j] += \
				S * np.dot(u_i, np.dot(J_fourier[i][j], u_j.conj())) - \
				S_j * np.dot(v_i, np.dot(J0_fourier[i][j], v_j))
			H[num_sites + i, num_sites + j] += \
				S * np.dot(u_i.conj(), np.dot(J_fourier[i][j], u_j)) - \
				S_j * np.dot(v_i, np.dot(J0_fourier[i][j], v_j))
			H[i, num_sites + j] += \
				S * np.dot(u_i, np.dot(J_fourier[i][j], u_j))
			H[num_sites + i, j] += \
				(S * np.dot(u_j, np.dot(J_fourier[j][i], u_i))).conj()

	#print("\nH =\n%s" % H)


	# trafo
	C = la.cholesky(H)
	signs = np.diag(np.concatenate((np.repeat(1, num_sites), np.repeat(-1, num_sites))))
	H_trafo = np.dot(C.transpose().conj(), np.dot(signs, C))

	#print("\nC =\n%s\n\nH_trafo =\n%s" % (C, H_trafo))


	# eigenvalues
	Es = np.real(la.eigvals(H_trafo))
	return Es


# plot
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
plt.xlabel("h (rlu)")
plt.ylabel("E (meV)")
plt.plot(hs, Es)
plt.show()
