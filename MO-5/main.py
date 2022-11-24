import random
import math
import numpy as np
import matplotlib.pyplot as pyplot


# Вариант 23
# Среднее гармоническое
# Манхэттен

def MovingAvgFind(K, M, FwithLine, alpha):
    MovingAvg = []
    for k in range(M, K - M + 1):
        sum = 0
        for j in range(k - M, k + M + 1):
            sum += alpha[j + M + 1 - k - 1] / FwithLine[j - 1]
        MovingAvg.append(sum)
    return MovingAvg


def AlphaMasGeneration(r, M):
    alphaMas = []
    for i in range(r):
        alphaMas.append(0)
    center = random.random()
    alphaMas[M] = center

    for m in range(M, r - 2):
        sum = 0
        for s in range(m, r - m):
            sum += alphaMas[s]
        center = random.uniform(0, 1 - sum) / 2
        alphaMas[m - 1] = center
        alphaMas[r - m] = center

    sum = 0
    for s in range(1, r - 1):
        sum += alphaMas[s]
    center = 0.5 * (1 - sum)
    alphaMas[0] = center
    alphaMas[r - 1] = center
    return alphaMas


def FindOmega(FwithLine, K, M):
    sum = 0
    for k in range(1, K + 1 - 2 * M):
        sum += abs(FwithLine[k] - FwithLine[k - 1])
    return sum


def FindDelta(FwithLine, FwithWave, K, M):
    sum = 0
    for k in range(0, K + 1 - 2 * M):
        sum += abs(FwithLine[k] - FwithWave[k])
    return sum / K


class Filtering(object):

    def __init__(self, signal, noise, lambdaMas, r):
        self.optimalAlpha = []
        self.omegas = []
        self.deltas = []
        self.Jmin = []
        self.distant = []
        self.r = r
        self.Signal = signal
        self.Noise = noise
        self.LambdaMas = lambdaMas

    def GarmonicMeanMethod(self, noise_arr, lambdaMas):
        r = self.r
        global K
        M = int((r - 1) / 2)
        N = int((math.log(1 - 0.95)) / math.log(1 - (0.01 / math.pi)))
        for i in range(len(lambdaMas)):
            alphaMas = AlphaMasGeneration(r, M)
            MovingAvg = MovingAvgFind(K, M, noise_arr, alphaMas)
            omega = FindOmega(MovingAvg, K, M)
            delta = FindDelta(MovingAvg, self.Noise, K, M)
            Jmin = lambdaMas[i] * omega + (1 - lambdaMas[i]) * delta
            for P in range(1, N):
                alphaMas = AlphaMasGeneration(r, M)
                MovingAvg = MovingAvgFind(K, M, noise_arr, alphaMas)
                omega = FindOmega(MovingAvg, K, M)
                delta = FindDelta(MovingAvg, self.Noise, K, M)
                J = lambdaMas[i] * omega + (1 - lambdaMas[i]) * delta
                if J < Jmin:
                    Jmin = J
                    bestAlphaMas = alphaMas.copy()
                    bestOmega = omega
                    bestDelta = delta
                MovingAvg.clear()
            self.Jmin.append(Jmin)
            self.optimalAlpha.append(bestAlphaMas)
            self.omegas.append(bestOmega)
            self.deltas.append(bestDelta)

    def dist(self):
        for index in range(0, len(self.omegas)):
            self.distant.append(abs(self.omegas[index]) + abs(self.deltas[index]))

    def OptimalSignSearch(self):
        r = self.r
        global K
        M = int((r - 1) / 2)
        minDist = min(self.distant)
        index = self.distant.index(minDist)
        optimalOmega = self.omegas[index]
        optimalDelta = self.deltas[index]
        optimalAlphas = self.optimalAlpha[index]
        optimalLambda = 0.1 * index
        optimalJ = optimalLambda * optimalOmega + (1 - optimalLambda) * optimalDelta

        filtering = []
        for k in range(M, K - M + 1):
            sum = 0
            for j in range(k - M, k + M + 1):
                sum += Noise[j - 1] * optimalAlphas[j + M + 1 - k - 1]
            filtering.append(sum)

        return optimalLambda, minDist, optimalAlphas, optimalOmega, optimalDelta, filtering, optimalJ

    def TextPrint(self, Lambda, omega, delta, J):
        r = self.r
        if r == 3:
            print("+-------+----------+----------------------+---------+---------+")
            print("|   h   |    dis   |         alpha        |    w    |    d    |")
            print("+-------+----------+----------------------+---------+---------+")
            for i in range(0, len(self.LambdaMas)):
                print(f"|  {self.LambdaMas[i]}", end="")
                print(f"  |{round(self.distant[i], 5):>9f}", end="")
                print(
                    f" | {self.optimalAlpha[i][0]:>0.4f} {self.optimalAlpha[i][1]:>0.4f} {self.optimalAlpha[i][2]:>0.4f} ",
                    end="")
                print(f"| {self.omegas[i]:>0.5f} |", end="")
                print(f" {self.deltas[i]:>0.5f} |")
            print("+-------+----------+----------------------+---------+---------+")
            print("+-------+---------+---------+---------+")
            print("|   h*  |    J    |    w    |     d   |")
            print("+-------+---------+---------+---------+")
            print(f"|  {Lambda:>0.1f}", end="")
            print(f"  | {J:>0.5f}", end="")
            print(f" | {omega:>0.5f} |", end="")
            print(f" {delta:>0.5f} |")
            print("+-------+---------+---------+---------+")
            print()

        if r == 5:
            print("+-------+----------+------------------------------------+---------+---------+")
            print("|   h   |    dis   |                alpha               |    w    |    d    |")
            print("+-------+----------+------------------------------------+---------+---------+")
            for i in range(0, len(self.LambdaMas)):
                print(f"|  {self.LambdaMas[i]}", end="")
                print(f"  |{round(self.distant[i], 5):>9f}", end="")
                print(
                    f" | {self.optimalAlpha[i][0]:>0.4f} {self.optimalAlpha[i][1]:>0.4f} {self.optimalAlpha[i][2]:>0.4f}"
                    f" {self.optimalAlpha[i][3]:>0.4f} {self.optimalAlpha[i][4]:>0.4f} ", end="")
                print(f"| {self.omegas[i]:>0.5f} |", end="")
                print(f" {self.deltas[i]:>0.5f} |")
            print("+-------+----------+------------------------------------+---------+---------+")
            print("+-------+---------+---------+---------+")
            print("|   h*  |    J    |    w    |     d   |")
            print("+-------+---------+---------+---------+")
            print(f"|  {Lambda:>0.1f}", end="")
            print(f"  | {J:>0.5f}", end="")
            print(f" | {omega:>0.5f} |", end="")
            print(f" {delta:>0.5f} |")
            print("+-------+---------+---------+---------+")
            print()

    def DotsPrint(self):

        # pyplot.scatter(self.omegas[0], self.deltas[0], color='black')
        pyplot.scatter(self.omegas[1], self.deltas[1], color='red')
        pyplot.scatter(self.omegas[2], self.deltas[2], color='orange')
        pyplot.scatter(self.omegas[3], self.deltas[3], color='yellow')
        pyplot.scatter(self.omegas[4], self.deltas[4], color='lime')
        pyplot.scatter(self.omegas[5], self.deltas[5], color='green')
        pyplot.scatter(self.omegas[6], self.deltas[6], color='skyblue')
        pyplot.scatter(self.omegas[7], self.deltas[7], color='blue')
        pyplot.scatter(self.omegas[8], self.deltas[8], color='purple')
        pyplot.scatter(self.omegas[9], self.deltas[9], color='pink')
        pyplot.scatter(self.omegas[10], self.deltas[10], color='brown')
        pyplot.legend([
            # f'{self.omegas[0]:>0.5f} {self.deltas[0]:>0.5f}',
            f'{self.omegas[1]:>0.5f} {self.deltas[1]:>0.5f}',
            f'{self.omegas[2]:>0.5f} {self.deltas[2]:>0.5f}',
            f'{self.omegas[3]:>0.5f} {self.deltas[3]:>0.5f}',
            f'{self.omegas[4]:>0.5f} {self.deltas[4]:>0.5f}',
            f'{self.omegas[5]:>0.5f} {self.deltas[5]:>0.5f}',
            f'{self.omegas[6]:>0.5f} {self.deltas[6]:>0.5f}',
            f'{self.omegas[7]:>0.5f} {self.deltas[7]:>0.5f}',
            f'{self.omegas[8]:>0.5f} {self.deltas[8]:>0.5f}',
            f'{self.omegas[9]:>0.5f} {self.deltas[9]:>0.5f}',
            f'{self.omegas[10]:>0.5f} {self.deltas[10]:>0.5f}'])
        pyplot.ylabel("delta")
        pyplot.xlabel("omega")
        pyplot.grid()
        pyplot.show()

    def PlotPrint(self, k, filteringSignal):
        global K
        global kFilter
        pyplot.plot(k, self.Signal, k, self.Noise, kFilter, filteringSignal)
        pyplot.title('Functions (r = ' + str(self.r) + ')')
        pyplot.ylabel('f(x)')
        pyplot.xlabel('x')
        pyplot.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
        pyplot.grid()
        pyplot.show()


K = 100
R = [3, 5]
k = np.array([float((k * math.pi) / K) for k in range(0, K + 1)])
Signal = np.array([(np.sin(elem) + 0.5) for elem in k])
Noise = np.array([(np.sin(elem) + 0.5) + (random.random() / 2 - 0.25) for elem in k])
LambaMas = np.array([l_elem / 10 for l_elem in range(0, 11)])

for r in R:
    filter = Filtering(Signal, Noise, LambaMas, r)
    filter.GarmonicMeanMethod(Noise, LambaMas)
    filter.dist()
    Lambda, dist, alphas, omega, delta, filteringSignal, J = filter.OptimalSignSearch()
    kFilter = np.array([float((k * math.pi) / K) for k in range(1, K + 3 - r)])

    filter.TextPrint(Lambda, omega, delta, J)
    filter.PlotPrint(k, filteringSignal)
    filter.DotsPrint()
