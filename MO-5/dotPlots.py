import matplotlib.pyplot as pyplot
import numpy as np

def DotsPrint(omegas, deltas, LambdaMas):
    # pyplot.scatter(omegas[0], deltas[0], color='black')
    pyplot.scatter(omegas[1], deltas[1], color='red')
    pyplot.scatter(omegas[2], deltas[2], color='orange')
    pyplot.scatter(omegas[3], deltas[3], color='yellow')
    pyplot.scatter(omegas[4], deltas[4], color='lime')
    pyplot.scatter(omegas[5], deltas[5], color='green')
    pyplot.scatter(omegas[6], deltas[6], color='skyblue')
    pyplot.scatter(omegas[7], deltas[7], color='blue')
    pyplot.scatter(omegas[8], deltas[8], color='purple')
    pyplot.scatter(omegas[9], deltas[9], color='pink')
    pyplot.scatter(omegas[10], deltas[10], color='brown')
    pyplot.legend([
        # f'lambda {LambdaMas[0]:>0.1f} ',
        f'lambda {LambdaMas[1]:>0.1f} ',
        f'lambda {LambdaMas[2]:>0.1f} ',
        f'lambda {LambdaMas[3]:>0.1f} ',
        f'lambda {LambdaMas[4]:>0.1f} ',
        f'lambda {LambdaMas[5]:>0.1f} ',
        f'lambda {LambdaMas[6]:>0.1f} ',
        f'lambda {LambdaMas[7]:>0.1f} ',
        f'lambda {LambdaMas[8]:>0.1f} ',
        f'lambda {LambdaMas[9]:>0.1f} ',
        f'lambda {LambdaMas[10]:>0.1f} '])
    pyplot.ylabel("delta")
    pyplot.xlabel("omega")
    pyplot.grid()
    pyplot.show()

l = np.array([l_elem / 10 for l_elem in range(0, 11)])

o1 = [13.29322, 5.99565, 5.99545, 5.99580, 5.99547, 5.99577, 5.99545, 5.99581, 5.99548, 5.99581, 5.99566]
d1 = [0.56380, 0.57609, 0.57611, 0.57616, 0.57610, 0.57616, 0.57611, 0.57616, 0.57610, 0.57616, 0.57614]

o2 = [9.18947, 3.57327, 3.56647, 3.57713, 3.57058, 3.56728, 3.58303, 3.59489, 3.59914, 3.59032, 3.58878]
d2 = [0.52912, 0.53525, 0.53545, 0.53533, 0.53558, 0.53532, 0.53567, 0.53523, 0.53496, 0.53508, 0.53584]


DotsPrint(o1, d1, l)
DotsPrint(o2, d2, l)
