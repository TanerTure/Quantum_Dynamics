import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('./')
import potentials

x_vals = np.linspace(-10,10,1000)
y_vals= potentials.abs(x_vals)

fig = (plt.figure(figsize=(6,3)))
plt.plot(x_vals,y_vals)
plt.show()
fig.savefig('abs_potential.png')
plt.clf()

y_vals_2=potentials.step(x_vals)

plt.plot(x_vals,y_vals_2)
plt.show()
fig.savefig('step_potential.png')
   

