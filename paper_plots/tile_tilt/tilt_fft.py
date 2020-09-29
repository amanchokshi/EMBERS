import numpy as np
from matplotlib import pyplot as plt

plt.style.use("seaborn")


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))


ang = np.radians(-10)
rot_mat = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])


x = np.linspace(-10, 10, 512)

gsn = gaussian(x, 10, 2)

c = np.dot(rot_mat, np.array([x, gsn]))

#  Mitch-isms
y = np.fft.fftshift(np.fft.fft(np.fft.fftshift(gsn))).imag
z = np.fft.fftshift(np.fft.fft(np.fft.fftshift(c[1]))).imag

plt.plot(x, gsn, label="Input Gaussian")
#  plt.plot(c[0], c[1], label="Rotated Gaussian")
plt.plot(x, y, label="Mitchism Gaussian")
#  plt.plot(x, z, label="Mitchism Rot")
ax = plt.gca()
#  ax.set_aspect(1)
plt.legend()
plt.tight_layout()
plt.show()
