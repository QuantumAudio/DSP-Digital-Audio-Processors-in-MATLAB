import wavio
import numpy as np

# Import input signal
inFile = wavio.read('original.wav')
s = inFile.data         # Input signal
fs = inFile.rate        # Sample rate
sw = inFile.sampwidth   # Sample width
N = s.size              # Number of samples
# Initial parameters
D1, D2 = int(0.25*fs), int(0.125*fs)
b0 = b1 = b2 = 1
a1, a2 = 0.2, 0.4
# Internal delay buffer 1 for s(n)
w1 = [0]*(D1 + 1)
# Internal delay buffer 2 for w1(n)
w2 = [0]*(D2 + 1)
# Delay buffer 1 index variable
q1 = 0
# Delay buffer 2 index variable
q2 = 0
# Delay buffer taps
tap1 = D1
tap2 = D2
# Empty list for output
y = np.array([])
# Loop through input signal
for n in range(0, N):
    s1 = w1[tap1]
    s2 = w2[tap2]
    y = np.append(y, b0*s[n] + b1*s1 + b2*s2)
    w2[q2] = s1 + a2*s2
    w1[q1] = s[n] + a1*s1
    q1 = q1 - 1             # Backshift index 1
    if q1 < 0:              # Circulate index 1
        q1 = D1
    tap1 = tap1 - 1         # Backshift tap1
    if tap1 < 0:            # Circulate tap1
        tap1 = D1
    q2 = q2 - 1             # Backshift index 2
    if q2 < 0:              # Circulate index 2
        q2 = D2
    tap2 = tap2 - 1         # Backshift tap2
    if tap2 < 0:            # Circulate tap2
        tap2 = D2

# Normalize results
y = y/y.max()
m = max(y)      # Maximum value of y
for i in range(len(y)):
    y[i] = y[i]/m

# Output results
wavio.write('output.wav', y, fs, sampwidth = sw)

