"""
MOST UP TO DATE CODE
Made as the Final for Statistical Mechanics II, coding a MD simulation of Argon from scratch

"""

import numpy as np
import matplotlib.pyplot as plt
import random
import time

###########################################################################################
# set initial conditions#
###########################################################################################

size = 4  # must be integer number >= 2

nmol = size ** 3  # number of atoms

print("\n\nnumbe of atoms = {}".format(nmol))

rmass = 39.9  # mass of Argon

# conversion of mass from amu to KPsA
rmass = rmass * 10. / 6.022169 / 1.380662

eps = 119.8  # LJ epsilon for Ar
eps6 = eps ** 6
sigma = 3.405  # LJ sigma for Ar
sig6 = sigma ** 6

rhos = 0.9  # rho star
rho = rhos / (sigma) ** 3  # density

T_star = 1.64  # Tstar
Temp_K = T_star * eps  # Temp

onethird = 1 / 3

cell = (nmol / rho) ** (onethird)  # cell length in angstroms
rcut = cell / 2  # cutoff distance
rcut6 = rcut ** 6

sigorcut = sigma / rcut

# value of the potential energy at r = rcut
ecut = 4 * eps * (sig6 / rcut6) * ((sig6 / rcut6) - 1)

dx = cell / size  

print("cell length = ", cell, "(A)")
print("cutoff = ", rcut, "(A)")
print("atoms are __{}(A)__ away\n".format(dx))
print("\nT* = ", T_star)
print("T(K) = ", Temp_K)

print("\np* = ", rhos)
print("p = ", rho)

# long range correction to the energy #
enerlrc = (8 / 3) * eps * np.pi * nmol
enerlrc = enerlrc * rhos * (onethird * sigorcut ** 9 - sigorcut ** 3)

print("Long Range Correction (LRC) to the energy = {}\n".format(enerlrc))
print("Reduced LRC to Energy = {}\n".format(enerlrc / eps))
print("Reduced LRC to Energy per particle = {}\n".format(enerlrc / eps / nmol))

###########################################################################################
# set up cubic lattice #
###########################################################################################
atoms = []
for i in range(size):
    for j in range(size):
        for k in range(size):
            pos = [i * dx, j * dx, k * dx]
            atoms.append(pos)

atoms = np.array(atoms)

# for center of mass scaling, is 3N bc 3DOF per atom
Nf = 3 * len(atoms)

# print("initial position array is \n", atoms, "\n")

###########################################################################################
# create gaussian random initial velocities #
###########################################################################################
atoms_vel = np.zeros_like(atoms)
np.random.seed(12343)
for i in range(len(atoms)):
    for j in range(3):
        atoms_vel[i][j] = np.random.normal(0, 0.5)


# print("initial velocities before scaling are:\n{} \n".format(atoms_vel))

def scale_temp(atoms_vel, target_temp):
    vsum = 0.
    for i in range(len(atoms_vel)):
        vsum += np.dot(atoms_vel[i], atoms_vel[i])
    instant_temp = rmass * vsum / Nf
    scale = np.sqrt(target_temp / instant_temp)
    atoms_vel = atoms_vel * scale

    return atoms_vel


atoms_vel = scale_temp(atoms_vel, Temp_K)
# print("initial scaled velocities are:\n {}\n".format(atoms_vel))

# need to subtract center of mass velocity
sumvx = 0.
sumvy = 0.
sumvz = 0.
for i in range(len(atoms_vel)):
    sumvx += atoms_vel[i][0]
    sumvy += atoms_vel[i][1]
    sumvz += atoms_vel[i][2]
comvx = sumvx / nmol
comvy = sumvy / nmol
comvz = sumvz / nmol
com_vel = np.array([comvx, comvy, comvz])
# print("center of mass velocity is \n",com_vel, "\n")

for i in range(len(atoms_vel)):
    atoms_vel[i] = atoms_vel[i] - com_vel
# print("velocities after subtracting COM_vel:\n{}\n".format(atoms_vel))

# now that we have lost 3 DOF we need to scale the new temperature with Nf = 3N -3
Nf = Nf - 3

atoms_vel = scale_temp(atoms_vel, Temp_K)


# print("final scaled velocities are:\n {}\n".format(atoms_vel))

###########################################################################################
# calculate forces/get energy/ get energy#
###########################################################################################

def calculate_forces(atoms, atoms_vel):
    # start force at 0
    forces = np.zeros_like(atoms)

    # Forces calculation on pair potential
    for i in range(0, len(atoms) - 1):
        for j in range(i + 1, len(atoms)):

            dx = np.array(atoms[i] - atoms[j])

            # PBC
            dx = dx - cell * np.round(dx / cell)

            rdiss = np.sqrt(np.dot(dx, dx))

            # only count forces w/in the cutoff
            if rdiss <= rcut:
                # LJ variables
                sor = sigma / rdiss
                sor2 = sor * sor
                sor6 = sor2 * sor2 * sor2
                sor12 = sor6 * sor6
                # force calculation
                ff = -24 * eps / rdiss * sor6 * (1 - 2 * sor6)
                forces[i] += ff * (dx / rdiss)
                forces[j] -= ff * (dx / rdiss)
    return forces


atom_forces = calculate_forces(atoms, atoms_vel)

# print("initial atom forces are:\n{}\n".format(atom_forces))

def get_temp(atoms_vel):
    sumv = 0.
    # kinetic energy
    for i in range(len(atoms_vel)):
        sumv += np.dot(atoms_vel[i], atoms_vel[i])

    instatemp = rmass * sumv / Nf
    return (instatemp)

print("initial starting temperature: {}".format(get_temp(atoms_vel)))


def get_ener(atoms, atoms_vel):
    KE = 0.
    PE = 0.
    vsum = 0.
    # for kinetic energy
    for i in range(len(atoms_vel)):
        vsum += np.dot(atoms_vel[i], atoms_vel[i])
        KE = 0.5 * rmass * vsum

    # for Potential Energy
    for i in range(0, len(atoms) - 1):
        for j in range(i + 1, len(atoms)):

            dx = np.array(atoms[i] - atoms[j])

            # PBC
            dx = dx - cell * np.round(dx / cell)

            rdiss = np.sqrt(np.dot(dx, dx))

            # only count forces w/in the cutoff
            if rdiss <= rcut:
                # LJ variables
                sor = sigma / rdiss
                sor2 = sor * sor
                sor6 = sor2 * sor2 * sor2
                sor12 = sor6 * sor6
                # force calculation
                PE += 4 * eps * (sor12 - sor6)
    enertot = KE + PE
    return (KE, PE, enertot)

ener_info = np.array(get_ener(atoms, atoms_vel))

print("initial kinetic energy: {}".format(ener_info[0]))
print("initial potential energy: {}".format(ener_info[1]))
print("initial TOTAL energy: {}".format(ener_info[2]))


###############################################################################################
# integrate position over time and write to file#
###############################################################################################
def integrate(atoms, atoms_vel, atom_forces, dt):
    # set old forces to zero
    old_forces = np.zeros_like(atoms)

    # set forces now as old forces
    old_forces = calculate_forces(atoms, atoms_vel)

    # 1. update all positions
    for i in range(len(atoms)):
        atoms[i] = atoms[i] + atoms_vel[i] * dt + 0.5 * (dt ** 2) * old_forces[i] / rmass

    # turn on/off Periodic Boundary conditions
    for k in range(len(atoms)):
        atoms[k] = atoms[k] - cell * np.round(atoms[k] / cell)

    # 2. calculate new forces at new positions
    new_forces = calculate_forces(atoms, atoms_vel)

    # 3. calculate new velocities based on new forces
    for i in range(len(atoms)):
        atoms_vel[i] = atoms_vel[i] + 0.5 * dt * (old_forces[i] + new_forces[i]) / rmass


def write_xyz(filename, atoms):
    out = open(filename, 'a')
    out.write("{}\n\n".format(len(atoms)))
    for i in range(len(atoms)):
        out.write("Ar {} {} {}\n".format(atoms[i][0], atoms[i][1], atoms[i][2]))
    out.close()


###############################################################################################
# MD LOOP#
###############################################################################################

print("###############################################################################################")
print("					START OF MD")
print("###############################################################################################\n")

nsteps = 5000 
vel_dist = []
ener_dist = []
pe_dist = []
ke_dist = []
temp_dist = []

start_time = time.time()

for step in range(nsteps):
    # get info every so often
    if step != 0 and step % 1000 == 0:
        # print("######################################")
        # print("step {} of {}".format(step, nsteps))
        print("{}% complete".format(round(step / (nsteps - 1) * 100)))

    # make sure to scale velocities every 10 steps for (fake) NVT
    # other wise, NVE
    if step != 0 and step % 5 == 0 and step < nsteps / 2:
        atoms_vel = scale_temp(atoms_vel, Temp_K)

    # integrate atoms over timestep and write to file
    integrate(atoms, atoms_vel, atom_forces, 0.01)
    # write_xyz("size3_60k.xyz", atoms)

    # collect energy info in last half
    if step >= (nsteps / 2):
        ener_info = np.array(get_ener(atoms, atoms_vel))
        instant_ke = ener_info[0]
        instant_pe = ener_info[1]
        instant_energy = ener_info[2]
        instant_temp = get_temp(atoms_vel)

        ke_dist.append(instant_ke)
        pe_dist.append(instant_pe)
        ener_dist.append(instant_energy)
        temp_dist.append(instant_temp)

###############################################################################################
# END OF SIMULATION INFO#
###############################################################################################
# plotting energy and temperature over time
plt.plot(range(round(nsteps / 2)), ener_dist, alpha=0.5, label="total energy")
plt.plot(range(round(nsteps / 2)), ke_dist, alpha=0.5, label="kinetic energy")
plt.plot(range(round(nsteps / 2)), pe_dist, alpha=0.5, label="potential energy")
plt.plot(range(round(nsteps / 2)), temp_dist, alpha=0.5, label="temperature")
plt.legend(bbox_to_anchor=(1, 1), loc='upper right')
plt.title("{} Ar atoms\nT* = {}\np* = {}".format(nmol, T_star, rhos))

# plt.hist(vel_dist, bins = 40)

print("\n\n######################################\nEND OF SIMULATION INFO\n######################################")
ave_pe = (sum(pe_dist) / len(pe_dist))

tim = (time.time() - start_time) / 60
tim  = "{:.2f}".format(tim)

udivn = ((ave_pe + enerlrc) / eps / nmol)
udivn = "{:.2f}".format(udivn)

lrcdivn = (enerlrc / eps / nmol)
lrcdivn = "{:.2f}".format(lrcdivn)


print("\nrhos\tT*\tN\tU*/N\tLRC/N\ttime (minutes)")

print(
    "{}\t{}\t{}\t{}\t{}\t{}".format(rhos, T_star, nmol, udivn, lrcdivn, tim))
print("######################################\nEND OF SIMULATION INFO\n######################################")

# plt.show()    
