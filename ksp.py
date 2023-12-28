import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
data = pd.read_csv("data1.csv")
parameters = ['Time',
              'Velocity',
              'GForce',
              'Acceleration',
              'Thrust',
              'TWR', 'Mass',
              'AltitudeFromTerrain',
              'AltitudeFromSea',
              'DownrangeDistance',
              'Latitude',
              'Longitude',
              'Apoapsis',
              'Periapsis',
              'Inclination',
              'OrbitalVelocity','TargetDistance',	'TargetVelocity',	'Stagevelocity_t',	'Vesselvelocity_t',	'Pressure']
for col in data.columns:
    if col not in parameters:
        data.drop(col, axis=1, inplace=True)
data.head()

V, h, g, F, m, alpha1, p = 0, 0, 9.8 / 5, 3749, 214, 0, 100
t = np.arange(0.0, 120, 0.1)
Vy, Vx, velocity_t, height_t = [0] * len(t), [0] * len(t), [0] * len(t), [0] * len(t)
mass_t, acceleration_t = [m] + ([0] * (len(t) - 1)), [5.54] + ([0] * (len(t) - 1))
force = [0] * len(t)
count = 0
alpha = [0] * len(t)
pressure = [0] * len(t)

for dot in t[1:]:
    dot_indx = int(dot * 10)
    rad_alpha = math.radians(alpha1 * math.pi)
    a_x = (F * math.cos(rad_alpha) - m * g * math.cos(rad_alpha) - 0.5 * 0.055 * p * 20 * (V ** 2) * math.cos(rad_alpha) * 1) / m
    a_y = (F * math.sin(rad_alpha) - m * g * math.sin(rad_alpha) - 0.5 * 0.055 * p * 20 * (V ** 2) * math.sin(rad_alpha) * 1) / m
    p = 100 * math.exp(-h / 5600)
    g = (6.67 * 10 ** -11) * 5.3 * 10 ** 22 / (600000 + h) ** 2
    h = (velocity_t[dot_indx - 1] ** 2) / (2 * acceleration_t[dot_indx - 1]) * (1.21 + 2.01 * math.cos(math.pi * count / (240 * 10)) * math.sin(math.pi * count / (240 * 10)))
    alpha[dot_indx] = alpha1
    alpha1 = math.sin(math.pi * count / (240 * 10))
    F += 1.3
    m -= 0.014
    force[dot_indx] = F
    pressure[dot_indx] = p
    count += 1
    Vx[dot_indx] = Vx[dot_indx - 1] + a_x * 0.1
    Vy[dot_indx] = Vy[dot_indx - 1] + a_y * 0.1
    velocity_t[dot_indx] = (Vy[dot_indx] ** 2 + Vx[dot_indx] ** 2) ** 0.5
    mass_t[dot_indx] = m - 0.126 * dot_indx
    
    acceleration_t[dot_indx] = (a_y ** 2 + a_x ** 2) ** 0.5
    height_t[dot_indx] = h
    print(alpha1)

velocity_t[0] = 0
acceleration_t[0] = 0
height_t[0] = 0
mass_t[0] = 214


plt.figure(figsize=(14, 12))
plt.subplot(221)
plt.xlabel(r'$t, с$')
plt.ylabel(r'$v, м/с$')
plt.plot(t, velocity_t)
plt.plot(data['Time'], data['Velocity'])
plt.grid(True)

plt.subplot(222)
plt.xlabel(r'$t, с$')
plt.ylabel(r'$a, м/с^2$')
plt.plot(t, acceleration_t)
plt.plot(data['Time'], data['Acceleration'])
plt.grid(True)

plt.subplot(223)
plt.xlabel(r'$t, с$')
plt.ylabel(r'$m, кг$')
plt.plot(t, mass_t)
plt.plot(data['Time'], data['Mass'])
plt.grid(True)

plt.subplot(224)
plt.xlabel(r'$t, с$')
plt.ylabel(r'$h, м$')
plt.plot(t, height_t)
plt.plot(data["Time"], list(map(lambda x: -1 * int(x) * 100, data["AltitudeFromSea"])))
plt.grid(True)
plt.show()