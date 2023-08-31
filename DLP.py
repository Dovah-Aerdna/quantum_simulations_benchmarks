#https://rawcdn.githack.com/rmhowe425/Qiskit/31620d163a12256be7c400e439bdefe09c53f9be/Final%20Product.html

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer
from qiskit.circuit.library import QFT
from qiskit.extensions import UnitaryGate
from numpy import log, sqrt, floor, pi, arange
import numpy as np

# Calculates the multiplicative inverse of x mod N
# (the number y such that xy = 1 (mod N)) using
# the extended Euclidean algorithm.
def invert(x, N):
    q = [0, 0]
    r = [N, x]
    t = [0, 1]

    while r[-1] != 0:
        q.append(r[-2]//r[-1])
        r.append(r[-2] - (q[-1]*r[-1]))
        t.append(t[-2] - (q[-1]*t[-1]))
    
    return int(t[-2] % N)

# Returns a unitary matrix which has the effect of multiplying each
# input |x> by a in mod N, resulting in the state |ax>.
def create_controlled_unitary(a, N):
    dim = 2**int(np.ceil(np.log(N)/np.log(2)) + 1)
    U = np.zeros((dim, dim))
    # Generate a permutation of the multiplicative group of Z_N.
    for i in range(int(dim/2)):
        U[i,i] = 1
    for i in range(N):
        U[int(dim/2) + i, ((a*i) % N)+int(dim/2)] = 1
    # The remaining states are irrelevant.
    for i in range(N, int(dim/2)):
        U[int(dim/2) + i, int(dim/2) + i] = 1
    return U

def create_unitary(a, N):
    a = int(np.round(a) % N)
    dim = 2**int(np.ceil(np.log(N)/np.log(2)))
    U = np.zeros((dim, dim))
    # Generate a permutation of the multiplicative group of Z_N.
    for i in range(N):
        U[i, ((a*i) % N)] = 1
    # The remaining states are irrelevant.
    for i in range(N, dim):
        U[i, i] = 1
    return U

# b is some power of a, and the oracle outputs the half-bit of m,
# where b = a^m (mod N), with >50% probability.
def oracle(a, b, N, verbose=False, draw=False):
    # Calculate the order of a
    r = 1
    product = a
    while product != 1:
        product = (product * a) % N
        r += 1

    # Find number of bits needed to store a value from 0 to N-1 (n)
    # and initialize 2 quantum registers of size n.
    n = int(np.ceil(np.log(N)/np.log(2)))
    qr1, qr2 = QuantumRegister(n), QuantumRegister(n)
    cr1, cr2 = ClassicalRegister(n), ClassicalRegister(1)
    qc = QuantumCircuit(qr1, qr2, cr1, cr2)
    
    #Change second register to state |00...01>
    qc.x(qr2[n-1])
    
    #Add H gate to first register
    for i in range(n):
        qc.h(qr1[i])
    
    # We need log_2(n) different matrices U_(a^(2^x))
    for i in range(n):
        U = create_controlled_unitary(a**(2**(n-i)) % N, N)
        qubits = [qr1[i]] + [qr2[j] for j in range(n)]
        qc.iso(U, qubits, [])

    qc.append(QFT(n), [qr1[i] for i in range(n)])

    for i in range(n):
        qc.measure(qr1[i], cr1[i])
    
    # Reuse the first quantum register, because we don't need it anymore.
    for i in range(2**(n-1), 2**n):
        qc.x(qr1[0]).c_if(cr1, i)

    qc.h(qr1[0])

    qc.barrier()

    # Now cr1 is in state y. We define k to be the closest integer to y*r / 2**n.
    # I don't think there's any way to get the result of the measurement mid-circuit
    # in qiskit. So this is a stopgap measure for now.
    for y in range(2**n):
        k = int(np.round(y*r/(2**n))) % r
        kInv = bin(invert(k, r))[2:]

        # Pad kInv with initial zeros
        while len(kInv) < n:
            kInv = '0' + kInv

        if '1' in kInv:
            for i in range(len(kInv)):
                bit = int(kInv[i])
                if bit == 1:
                    # Apply U operation only if the value of cr1 is y.
                    U = create_unitary(b**(2**i) % N, N)
                    qc.iso(U, [qr2[i] for i in range(n)], []).c_if(cr1, y)
    
    qc.barrier()
    qc.rz(-np.pi/2 , qr1[0])
    qc.h(qr1[0])
    qc.measure(qr1[0], cr2[0])
    
    if draw:
        print(qc.draw(output='text'))

    # Execute the circuit on a simulator.
    backend_name = 'qasm_simulator'
    backend = Aer.get_backend(backend_name)
    #backend.set_options(device='GPU')
    if verbose:
        print("Running circuit on", backend_name, "...")
    result = execute(qc, backend, shots=int(np.ceil(8*np.pi))).result().get_counts(qc)

    if verbose:
        print(result)
    
    # See whether the half-bit was measured more often as 0 or as 1, and
    # return that value as the half-bit.
    zeros_count = 0
    ones_count = 0

    for k in result.keys():
        half_bit = k[0]
        if half_bit == '0':
            zeros_count += result[k]
        else:
            ones_count += result[k]

    if verbose:
        print('Zeros:', zeros_count, '\tOnes:', ones_count)
    
    if zeros_count > ones_count:
        return 0
    else:
        return 1
    


from random import randint, random

e = 1/pi   # Maximum epsilon from the Magic Box paper.

def calculate_order(a, N):
    power = 1
    curr = a
    while (curr != 1):
        curr = (curr * a) % N
        power += 1
    return power

'''
   Determines the most significant bit of an integer c
   @param c: Base
   @param n: Modulus
   @return: 1 or -1
'''
def n1(c, n):
    val = -1

    if c % n >= (n / 2):
        val = 1
        
    return val

'''
   Decides which 'guess' best fits with the output of the oracle using cross correlation.
   @param G: Base of the logarithm.
   @param X: Power of the logarithm.
   @param n: Modulus.
   @param l: Lag of cross correlation - refers to how far the series are offset.
   @param d: The probability that the estimation will be incorrect.
'''
def EstimateCrossCorrelation(G, X, n, l, d, period):
    total = 0  # Called 'sum' in algorithm

    # Compute the number of trials
    m = int(round(2 / (sqrt(d) * e)))

    # Compute estimate
    for trial in range(m):
        t = randint(1, n)
        output = oracle(G, X, n)  # Oracle(base, LHS, mod)

        if output == 0:
            output = -1
        
        total = total + (output*(X * (G ** period)) * n1(t + l, n))

    return (total / m)

'''
   Main algorithm for computing the discrete logarithm.
   @param G: Base of a the logarithm.
   @param X: Power of the logarithm.
   @param n: Modulus
'''
def Logarithm(G, X, n):
    step = (n * e)
    l = int(log(n)/log(2))  # Number of iterations
    d = round((l / 4), 10)  # Limit on probability error

    period = calculate_order(G, n)
    itw=0
    itf=0
    # Repeat until logarithm is found.
    while True:
        print("iter wh=",itw)
        itw=itw+1
        # Make initial guess. 
        for c in arange(0, n-1+step, step):
            
            # Refine guess.
            for i in range(l-1, -1, -1):
                print("iter for=",itf)
                itf=itf+1
                
                Xi = (X ** (2**i)) % n
                if (EstimateCrossCorrelation(G, Xi, n, c/2, d, period)
                    > EstimateCrossCorrelation(G, Xi, n, c/2 + n/2, d, period)):
                    print("if")
                    c = (c/2) % period
                else:
                    print("else")
                    c = (c/2 + n/2) % period
            
            potential = int(floor(c))
            if ((G ** potential) % n) == X:
                print("ret potential")
                return potential
            
def is_primitive_root(g, p):
    tracker = []
    for i in range(1, p):
        eq = pow(g, i) % p
        if (eq in tracker):
            #print(f"{g} is NOT a primitive root of {p}.")
            return False
            break
        else:
            tracker.append(eq)
        if (i == p-1):
            #print(f"{g} is a primitive root of {p}.")
            return True
        


import time
import sympy 
import random

'''
start = time.time()
print(Logarithm(7, 13, 15)) #7^x=13 mod 15
end = time.time()

print("Execution took", int(round(end - start)), "seconds.")
'''
#############################


list_p=sympy.primerange(10,12)

print("p\tbitlen(p)\tg\tpk\ttime")

for p in list_p:
    for g in range(2,p):
        if is_primitive_root(g, p):
            
            pk=random.randint(2,p-2)
            print("g=",g," p=",p," pk=",pk)
            start = time.time()
            x=Logarithm(g, pk, p) #g^x=pk mod p
            end = time.time()
            print(p,"\t",p.bit_length(),"\t",g,"\t",pk,"\t",(end-start))
