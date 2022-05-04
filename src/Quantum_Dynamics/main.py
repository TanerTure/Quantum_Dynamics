import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('./')
import potentials
import matrix_numerov

x_vals = np.linspace(-10, 10, 1000)
y_vals = potentials.abs(x_vals)

fig = (plt.figure(figsize=(6, 3)))
plt.plot(x_vals, y_vals)
plt.show()
fig.savefig('output/abs_potential.png')
plt.clf()

evals, evecs = matrix_numerov.matrix_numerov(x_vals, y_vals)
for i in range(5):
    plt.plot(evecs[:,i])
fig.savefig('output/matrix_numerov_1.png')
with open('output/numerov_energies_1.txt', 'w') as output_file:
    output_file.write(evals.tostring())
plt.clf()

evals, evecs = matrix_numerov.matrix_numerov(x_vals, y_vals, order=2)
for i in range(5):
    plt.plot(evecs[:,i])
fig.savefig('output/matrix_numerov_2.png')
with open('output/numerov_energies_2.txt', 'w') as output_file:
    output_file.write(evals.tostring())
plt.clf()

evals, evecs = matrix_numerov.matrix_numerov(x_vals, y_vals, order=3)
for i in range(5):
    plt.plot(evecs[:,i])
fig.savefig('output/matrix_numerov_3.png')
with open('output/numerov_energies_3.txt', 'w') as output_file:
    output_file.write(evals.to_string())
plt.clf()


y_vals_2 = potentials.step(x_vals)

plt.plot(x_vals, y_vals_2)
plt.show()
fig.savefig('output/step_potential.png')
plt.clf()

y_vals_3 = potentials.bathtub(x_vals) 
plt.plot(x_vals,y_vals_3)
fig.savefig('output/bathtub_potential.png')
plt.clf()


y_vals_4=potentials.harmonic(x_vals)
plt.plot(x_vals,y_vals_4)
fig.savefig('output/harmonic_oscillator.png')
plt.clf()



