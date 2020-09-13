import wave
# Create Wave_read object
wav = wave.open('original.wav', mode='rb')
# Get number of frames, N
N = wav.getnframes()
# Get sample rate
fs = wav.getframerate() # sample rate is 16 kHz
# Read in audio file with echo
s = wav.readframes(N)

# Initial parameters
D = int(0.25*fs)
a = 0.45
# Internal delay buffer for s(n)
w = [0]*(3*D + 1)
# Delay buffer index variable
q = 0
# Delay buffer taps
tap1 = D
tap2 = 2*D
tap3 = 3*D
# Create empty list for output
y = []
# Loop through input signal
for n in range(0, N):
    # Read input into w
    w[q] = s[n]
    # y(n) = s(n) + as(n-D) + a^2s(n-2D) + a^3s(n-3D)
    y.append(w[q] + a*w[tap1] + a**2*w[tap2] + a**3*w[tap3])
    q = q - 1           # Backshift index
    if q < 0:           # Circulate index
        q = 3*D
    tap1 = tap1 - 1     # Backshift tap1
    if tap1 < 0:        # Circulate tap1
        tap1 = 3*D
    tap2 = tap2 - 1     # Backshift tap2
    if tap2 < 0:        # Circulate tap2
        tap2 = 3*D
    tap2 = tap3 - 1     # Backshift tap3
    if tap3 < 0:        # Circulate tap3
        tap3 = 3*D

# Normalize results
m = max(y)      # Maximum value of y
o = []          # Empty list for output
for i in y:
    o.append(i/m)

# Output results to file
out = wave.open('out.wav', mode='wb')
out.setsampwidth((wav.getsampwidth()))
out.setnchannels(wav.getnchannels())
out.setframerate(fs)
out.setnframes(N)
out.writeframes(o)

# Close files
wav.close()
out.close()
