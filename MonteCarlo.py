import numpy as np
import matplotlib.pyplot as plt

# System parameters
N = 1000              # Number of particles. Larger N improves statistical stability; smaller N causes more energy fluctuations.
m = 1.67e-27          # Mass of a hydrogen atom (in SI units)
k_B = 1.38e-23        # Boltzmann constant (SI)
v_limit = 10          # Initial velocity range. Affects initial distribution, but not long-term behavior under Metropolis.
T = 300               # Temperature in Kelvin. Controls thermal energy.
beta = 1 / (k_B * T)  # Inverse thermal energy
dimensions = 2        # Dimensionality of the system (e.g. 1D, 2D, 3D)

# Function to compute total kinetic energy of a given microstate
def total_energy(momentum, mass):
    return np.sum(momentum**2) / (2 * mass)

# Initial random velocities for each particle and direction
velocities = np.random.uniform(low=-v_limit, high=v_limit, size=(N, dimensions))

# Initialize list to store the total energy of each microstate
momentum = m * velocities
energy_list = [total_energy(momentum, m)]

# Track iteration count for plotting purposes
iteration_list = [0]

# Metropolis algorithm: generate new microstates via small random perturbations
for step in range(5_000_000):
    idx = np.random.randint(0, N)  # Choose a random particle

    # Propose a small velocity change with Gaussian distribution
    dv = np.random.normal(loc=0, scale=100, size=dimensions)

    velocities[idx, :] += dv
    momentum = m * velocities
    new_E = total_energy(momentum, m)
    old_E = energy_list[-1]

    # Metropolis acceptance criterion
    if new_E > old_E:
        rand = np.random.uniform()
        if rand <= np.exp(-beta * (new_E - old_E)):
            energy_list.append(new_E)
        else:
            velocities[idx, :] -= dv
            energy_list.append(old_E)
    else:
        energy_list.append(new_E)

    iteration_list.append(step)

# Theoretical internal energy for an ideal monoatomic gas
U_theoretical = N * k_B * T * dimensions / 2

# Accept only microstates close enough to equilibrium energy
valid_energies = []
record = False
error_margin = 0.05
first_valid_index = 0

for E in energy_list:
    rel_error = abs(E - U_theoretical) / abs(U_theoretical)
    if not record and rel_error <= error_margin:
        record = True
        valid_energies.append(E)
        first_valid_index = energy_list.index(E)
    elif record:
        valid_energies.append(E)

# Compute average internal energy and its relative error
U_exp = np.mean(valid_energies)
U_rel_error = abs(U_exp - U_theoretical) / abs(U_theoretical) * 100

# Compute experimental and theoretical heat capacity and error
C_v_exp = (np.mean(np.array(valid_energies)**2) - np.mean(valid_energies)**2) / (k_B * T**2)
C_v_theoretical = U_theoretical / T
C_v_rel_error = abs(C_v_exp - C_v_theoretical) / C_v_theoretical * 100

# Print results
print(f"Experimental internal energy: {U_exp:.3e} J")
print(f"Theoretical internal energy: {U_theoretical:.3e} J")
print(f"Relative error in internal energy: {U_rel_error:.4f}%")
print(f"Experimental heat capacity: {C_v_exp:.3e} J/K")
print(f"Theoretical heat capacity: {C_v_theoretical:.3e} J/K")
print(f"Relative error in heat capacity: {C_v_rel_error:.4f}%")

# Plot energy evolution and mark valid states
plt.plot(iteration_list, energy_list, label="Total energy of each microstate")
valid_iterations = iteration_list[first_valid_index:]
plt.plot(valid_iterations, valid_energies, color='green', label="Microstates used for C$_v$")
plt.axhline(U_theoretical, linestyle='--', color='red', label="Theoretical internal energy")
plt.xlabel("Number of iterations")
plt.ylabel("Energy (J)")
plt.grid()
plt.legend(loc='lower right')
plt.show()
