import json
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


nside = 256
data = "./interp_maps/pix_resi.json"

delta_az = 83.1912

p_0 = [0, 0]
p_2 = [90, 6.80880000000001]
p_4 = [270, 6.80880000000001]

with open(data) as f:
    tilt_resi =  json.load(f)

# pixel of minimum residual power
phase_cen = {}

for tile in tilt_resi.keys():
    if "rf1" and "_0" in tile:
    #  if "S06" in tile:
        min_pix = np.argmax(tilt_resi[tile])
        #  if min_pix <= 600:
            #  phase_cen[f"{tile}"] = np.degrees(np.array(hp.pix2ang(nside, min_pix)))
        phase_cen[f"{tile}"] = np.array(hp.pix2ang(nside, min_pix))

theta = []
phi = []
for i in list(phase_cen.values()):
    theta.append(i[0])
    phi.append(i[1])


plt.style.use("seaborn")
figure = plt.figure(figsize=(6, 6))
ax = figure.add_subplot(111, polar=True)

#  ax.set_ylim(0, 12)
#  ax.set_rgrids([0, 2, 4, 6, 8, 10], angle=22)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

#  plt.scatter(np.radians([p_0[0], p_2[0], p_4[0]]), [p_0[1], p_2[1], p_4[1]] , color="#529471")
plt.scatter(phi, np.degrees(theta), color="#529471")
plt.tight_layout()
plt.show()

for i in range(len(phi)):
    print(np.degrees(phi[i]), np.degrees(theta[i]))
