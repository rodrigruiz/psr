import numpy as np
import matplotlib.pyplot as plt

x = np.arange(58000, 59000, 0.1)

def y(x, a, b, c):
    return a * np.exp(- (x + b) * c)

fig, ax = plt.subplots()
plt.plot(x, y(x, 7 * 10**2, 58000, np.log(2)/500) - 7 * 10**2 * 0.96)
plt.show()
fig.savefig('test.png')