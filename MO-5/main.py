import random
import math
import numpy as np
import matplotlib.pyplot as plt


def FindOmega(FwithLine, M):
    sum = 0
    for k in range(1, 101 - 2 * M):
        sum += abs(FwithLine[k] - FwithLine[k - 1])
    return sum


def FindAlpha(FwithLine, FwithWave, M):
    sum = 0
    for k in range(0, 101 - 2 * M):
        sum += abs(FwithLine[k] - FwithWave[k])
    return sum / 100


def AlphaMasGeneration(r, M):
    if r == 3:
        alphaMas = np.array([0, 0, 0], float)
    if r == 5:
        alphaMas = np.array([0, 0, 0, 0, 0], float)
    center = random.random()
    alphaMas[M] = center

    for m in range(M, r - 2):
        summary = 0
        for s in range(m, r - m):
            summary += alphaMas[s]
        center = random.uniform(0, 1 - summary) / 2
        alphaMas[m - 1] = center
        alphaMas[r - m] = center

    summary = 0
    for s in range(1, r - 1):
        summary += alphaMas[s]
    center = 0.5 * (1 - summary)
    alphaMas[0] = center
    alphaMas[r - 1] = center
    return alphaMas


def Filtering(M, FwithLine, alpha):
    FwithLineSred = []
    for k in range(M, 100 - M + 1):
        summary = 0
        for j in range(k - M, k + M + 1):
            summary += alpha[j + M + 1 - k - 1] / FwithLine[j - 1]
        FwithLineSred.append(summary)
    return FwithLineSred


def signalGen(num):
    return np.sin(num) + 0.5


class solvingClass(object):

    def __init__(self, _signal, _noise, _lambda_arr, r):
        self.optimal_a = []
        self.omegas = []
        self.deltas = []
        self.min_J = []
        self.distant = []
        self.r = r
        self.signal = _signal
        self.noise = _noise
        self.lambdas = _lambda_arr

    def arithmeticMeanMethod(self, noise_arr, lambda_arr):
        r = self.r
        moving_aver = []
        M = int((r - 1) / 2)
        N = int((math.log(1 - 0.95)) / math.log(1 - (0.01 / math.pi)))
        for i in range(len(lambda_arr)):
            best_a = np.array([])
            best_omega = 0
            best_delta = 0
            min_J = 10000
            for P in range(0, N):
                a_arr = AlphaMasGeneration(r, M)
                moving_aver = Filtering(M, noise_arr, a_arr)
                omega = FindOmega(moving_aver, M)
                delta = FindAlpha(moving_aver, self.noise, M)
                J = lambda_arr[i] * omega + (1 - lambda_arr[i]) * delta
                if J < min_J:
                    min_J = J
                    best_a = a_arr.copy()
                    best_omega = omega
                    best_delta = delta
                moving_aver.clear()
            self.min_J.append(min_J)
            self.optimal_a.append(best_a)
            self.omegas.append(best_omega)
            self.deltas.append(best_delta)
        print('aboba')

    def distantSearch(self):
        for index in range(0, len(self.omegas)):
            self.distant.append(abs(self.omegas[index]) + abs(self.deltas[index]))

    def optimalSignSearch(self):
        r = self.r
        M = int((r - 1) / 2)
        min_dist = min(self.distant)
        index = self.distant.index(min_dist)
        optimal_omega = self.omegas[index]
        optimal_delta = self.deltas[index]
        optimal_alphas = self.optimal_a[index]
        optimal_lambda = 0.1 * index
        optimal_J = optimal_lambda * optimal_omega + (1 - optimal_lambda) * optimal_delta

        filtering = []
        for k_elem in range(M, 100 - M + 1):
            summary = 0
            for j in range(k_elem - M, k_elem + M + 1):
                summary += noiseUnit[j - 1] * optimal_alphas[j + M + 1 - k_elem - 1]
            filtering.append(summary)

        return optimal_lambda, min_dist, optimal_alphas, optimal_omega, optimal_delta, filtering, optimal_J

    def print(self, o_lambda, omega, delta, J):
        r = self.r
        if r == 3:
            print("+-------+----------+----------------------+---------+---------+")
            print("|   h   |    dis   |         alpha        |    w    |    d    |")
            print("+-------+----------+----------------------+---------+---------+")
            for i in range(0, len(self.lambdas)):
                print(f"|  {self.lambdas[i]}", end="")
                print(f"  |{round(self.distant[i], 5):>9f}", end="")
                print(f" | {self.optimal_a[i][0]:>0.4f} {self.optimal_a[i][1]:>0.4f} {self.optimal_a[i][2]:>0.4f} ",
                      end="")
                print(f"| {self.omegas[i]:>0.5f} |", end="")
                print(f" {self.deltas[i]:>0.5f} |")
                print("+-------+----------+----------------------+---------+---------+")

        if r == 5:
            print("+-------+----------+------------------------------------+---------+---------+")
            print("|   h   |    dis   |                alpha               |    w    |    d    |")
            print("+-------+----------+------------------------------------+---------+---------+")
            for i in range(0, len(self.lambdas)):
                print(f"|  {self.lambdas[i]}", end="")
                print(f"  |{round(self.distant[i], 5):>9f}", end="")
                print(f" | {self.optimal_a[i][0]:>0.4f} {self.optimal_a[i][1]:>0.4f} {self.optimal_a[i][2]:>0.4f}"
                      f" {self.optimal_a[i][3]:>0.4f} {self.optimal_a[i][4]:>0.4f} ", end="")
                print(f"| {self.omegas[i]:>0.5f} |", end="")
                print(f" {self.deltas[i]:>0.5f} |")
                print("+-------+----------+------------------------------------+---------+---------+")
            print("+-------+---------+---------+---------+")
            print("|   h*  |    J    |    w    |     d   |")
            print("+-------+---------+---------+---------+")
            print(f"|  {o_lambda:>0.1f}", end="")
            print(f"  | {J:>0.5f}", end="")
            print(f" | {omega:>0.5f} |", end="")
            print(f" {delta:>0.5f} |")
            print("+-------+---------+---------+---------+")
            print()

    def dotsPlot(self):

        # plt.scatter(0, 0, color='purple')
        plt.scatter(self.omegas[0], self.deltas[0], color='blue')
        plt.scatter(self.omegas[1], self.deltas[1], color='green')
        plt.scatter(self.omegas[2], self.deltas[2], color='grey')
        plt.scatter(self.omegas[3], self.deltas[3], color='black')
        plt.scatter(self.omegas[4], self.deltas[4], color='pink')
        plt.scatter(self.omegas[5], self.deltas[5], color='yellow')
        plt.scatter(self.omegas[6], self.deltas[6], color='lime')
        plt.scatter(self.omegas[7], self.deltas[7], color='#ad09a3')
        plt.scatter(self.omegas[8], self.deltas[8], color='orange')
        plt.scatter(self.omegas[9], self.deltas[9], color='red')
        plt.scatter(self.omegas[10], self.deltas[10], color='skyblue')
        plt.legend([
            f'{self.omegas[0]:>0.5f} {self.deltas[0]:>0.5f}', f'{self.omegas[1]:>0.5f} {self.deltas[1]:>0.5f}',
            f'{self.omegas[2]:>0.5f} {self.deltas[2]:>0.5f}', f'{self.omegas[3]:>0.5f} {self.deltas[3]:>0.5f}',
            f'{self.omegas[4]:>0.5f} {self.deltas[4]:>0.5f}', f'{self.omegas[5]:>0.5f} {self.deltas[5]:>0.5f}',
            f'{self.omegas[6]:>0.5f} {self.deltas[6]:>0.5f}', f'{self.omegas[7]:>0.5f} {self.deltas[7]:>0.5f}',
            f'{self.omegas[8]:>0.5f} {self.deltas[8]:>0.5f}', f'{self.omegas[9]:>0.5f} {self.deltas[9]:>0.5f}',
            f'{self.omegas[10]:>0.5f} {self.deltas[10]:>0.5f}'])
        plt.ylabel("delta")
        plt.xlabel("omega")
        plt.grid()
        plt.show()
        """        
        plt.scatter(self.omegas, self.deltas)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid()
        plt.show()
        """


k = np.array([float((k * math.pi) / 100) for k in range(0, 101)])
signalUnit = np.array([signalGen(elem) for elem in k])
noiseUnit = np.array([signalGen(elem) + (random.random() / 2 - 0.25) for elem in k])
lambda_l = np.array([l_elem / 10 for l_elem in range(0, 11)])

r = 3
mov_wind3 = solvingClass(signalUnit, noiseUnit, lambda_l, r)
mov_wind3.arithmeticMeanMethod(noiseUnit, lambda_l)
mov_wind3.distantSearch()
o3_lambda, o3_dist, o3_alphas, o3_omega, o3_delta, o3_filtering_signal, o3_J = mov_wind3.optimalSignSearch()
mov_wind3.print(o3_lambda, o3_omega, o3_delta, o3_J)

filter_k = np.array([float((k * math.pi) / 100) for k in range(1, 100)])
plt.plot(k, signalUnit, k, noiseUnit, filter_k, o3_filtering_signal)
plt.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Размер скользящего окна r = 3')
plt.grid()
plt.show()

mov_wind3.dotsPlot()

r = 5
mov_wind5 = solvingClass(signalUnit, noiseUnit, lambda_l, 5)
mov_wind5.arithmeticMeanMethod(noiseUnit, lambda_l)
mov_wind5.distantSearch()
o5_lambda, o5_dist, o5_alphas, o5_omega, o5_delta, o5_filtering_signal, o5_J = mov_wind5.optimalSignSearch()
mov_wind5.print(o5_lambda, o5_omega, o5_delta, o5_J)

filter_k5 = np.array([float((k * math.pi) / 100) for k in range(1, 98)])
plt.plot(k, signalUnit, k, noiseUnit, filter_k5, o5_filtering_signal)
plt.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Размер скользящего окна r = 5')
plt.grid()
plt.show()

mov_wind5.dotsPlot()
