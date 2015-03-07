# -- 1. Import libraries
import matplotlib.pyplot as plt
from scipy import *
import usadel1 as u

# -- 2. Specify the geometry

geometry = u.Geometry(nwire=1, nnode=2)

geometry.t_type = [u.NODE_CLEAN_S_TERMINAL, u.NODE_CLEAN_S_TERMINAL]
geometry.w_type = [u.WIRE_TYPE_N]

geometry.w_ends[0,:] = [0, 1]

geometry.t_delta = [100, 100]
geometry.t_phase = [-.25*pi, .25*pi]

geometry.w_length = 1
geometry.w_conductance = 1

# -- 3. Solve the DOS

solver = u.CurrentSolver(geometry)

solver.solve_spectral()
solver.save('sns-spectral.h5')

a, b = solver.spectral.a, solver.spectral.b
E, x = solver.spectral.E, solver.spectral.x
dos = real((1 + a*b)/(1 - a*b))

# Plot it

j = x[::2] * 101
plt.plot(E[:,None] - 0.5*j[None,:], dos[:,0,::2] + 0.08*j[None,:], 'k-')
plt.xlabel('$E/E_T$'); plt.ylabel('$n/n_N$')
plt.ylim(0, 15); plt.xlim(-50, 300)
plt.savefig('sns-dos.eps')

# -- 4. Compute T=0 supercurrent vs. Delta

Deltas = logspace(0,2,15)
I_S = empty_like(Deltas)

for j, Delta in enumerate(Deltas):
    geometry.t_delta = Delta
    geometry.t_t = 1e-5
    geometry.t_phase = array([-.5, .5]) * pi/2 * 1.26
    solver.solve()
    I_S[j] = solver.get_currents(ix=0)[0]

savetxt('sns-I_S-vs-Delta.dat', c_[Deltas, I_S])

geometry.t_delta = 100

# Plot it

plt.clf()
plt.plot(Deltas, I_S, 'k-')
plt.xlabel('$\Delta/E_T$'); plt.ylabel('$e R_N I_S / E_T$')
plt.savefig('sns-I_S-vs-Delta.eps')

# -- 5. Solve the current-phase relation

phi = linspace(0, pi, 50)
I_S = zeros([50])

geometry.t_t = 1e-6 # Zero temperature

solver.set_solvers(sp_solver=u.SP_SOLVER_TWPBVP)

for j, p in enumerate(phi):
    geometry.t_phase = array([-.5, .5]) * p
    solver.solve()
    Ic, Ie = solver.get_currents(ix=0)
    I_S[j] = Ic[0]

savetxt('sns-I_S-vs-phi.dat', c_[phi, I_S])

# Plot it

plt.clf()
plt.plot(phi/pi, I_S, 'k-')
plt.xlabel('$\phi/\pi$'); plt.ylabel('$e R_N I_S / E_T$')
plt.savefig('sns-I_S-vs-phi.eps')
