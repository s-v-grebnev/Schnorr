"""
How to run:
$sage fac2.sage N n prec trials
    :param N: the number to factor
    :param n: dimension of lattice
    :param prec: precision (integer to correct rounding errors)
    :return: prints resulting nearest-to-short vector t, aux vector c
"""

from fpylll import IntegerMatrix, SVP, LLL
import sys

def svp(B):
	A = IntegerMatrix.from_matrix(B)
	return SVP.shortest_vector(A, pruning=False)

def cvp(B, t):
	A = IntegerMatrix.from_matrix(B)
	A = LLL.reduction(A, delta=.75)		# LLL штатный
#	print(A)
#Babai
	t = SVP.shortest_vector(A, pruning=False)	# 'target'
#	print(t)
	b = vector(ZZ, t)
	A = matrix(ZZ, A)
	tt = []
	def my_rnd(a):
		if a > 0:
			return 1
		elif a < 0:
			return -1
		return 0	
	for j in range((A.nrows()-1), 0, -1):
#		c = (1. * b.dot_product(A.row(j)) / A.row(j).dot_product(A.row(j)))
		c = my_rnd(1. * b.dot_product(A.row(j)) / A.row(j).dot_product(A.row(j)))
		b = b - c * A[j] 	# check 
		tt.append(c)
	b = vector(ZZ, t) - b
	return b, vector(ZZ, tt)	

def first_primes(n):
	p = 1
	P = []
	while len(P) < n:
		p = next_prime(p)
		P += [p]
	return P

def is_smooth(x, P):
	y = x
	for p in P:
		while p.divides(y):
			y /= p
	return abs(y) == 1

# Test if a factoring relation was indeed found.
def test_Schnorr(P, N, n, prec=10):
	shuffle(f)

	# Scale up and round
	def sr(x):
		return round(x * 2^prec)

	diag = [sr(N*f[i]) for i in range(n)] + [sr(N*ln(N))]
	B = diagonal_matrix(diag, sparse=False)
	for i in range(n):
		B[i, n] = sr(N*ln(P[i]))


	print(B)
# Find a "short" vector
	b = svp(B)
# Reduce the vector by Babai's algorithm		
	a = cvp(B, b)
	print(a)

"""
	e = [b[i] / sr(N*f[i]) for i in range(n)]
	u = 1
	v = 1
	for i in range(n):
		assert e[i] in ZZ
		if e[i] > 0:
			u *= P[i]^e[i]
		if e[i] < 0:
			v *= P[i]^(-e[i])
	return is_smooth(u - v*N, P) 
"""

try:
	N = int(sys.argv[1])
except:
	N = 48567227

try:
	n = int(sys.argv[2])
except:
	n = 47

try:
	pprec = int(sys.argv[3])
except:
	pprec = 4

try:
	trials = int(sys.argv[4])
except:
	trials = 1


N = 48567227
P = first_primes(n)
f = list(range(1, n+1))
for i in range(trials):
	test_Schnorr(P, N, n, prec = pprec)

