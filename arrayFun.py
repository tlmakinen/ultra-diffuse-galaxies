# arrayFun
import numpy as np
import math

a = np.zeros((16,64))

def fib_formula(n):
	ratio = (1 + math.sqrt(5)) / 2
	value = (ratio**n - (1 - ratio)**n) / math.sqrt(5)
	return int(round(value))

b = np.zeros(16)

b = np.vstack(b)

for indx in (0,16):
	b[indx] = fib_formula(indx)

print(b)
