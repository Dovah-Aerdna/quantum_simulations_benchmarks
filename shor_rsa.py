import numpy as np
import qiskit 
from math import sqrt,log,gcd
import random
from random import randint
import rsa
import sympy
import time
from qiskit.algorithms import Shor #use qiskit 0.36

def mod_inverse(a, m):
    for x in range(1, m):
        if (a * x) % m == 1:
            return x
    return -1


def generate_my_keys(list_n,verbose=False):
    
    keys=[] #((e, n), (d, n))

    for n, bitl, lam in list_n:

        e=lam+1
        for k in range(int(np.floor(np.log(lam)))):
            if sympy.prevprime(e)>3:
                e=sympy.prevprime(e)
            else:
                break
        
        while(gcd(e,lam)!=1):
            e=sympy.prevprime(e)
            
        d=modular_inverse(e,lam)
        if verbose:
            print("n=",n," e=",e, " d=",d, " lam=",lam)
        keys.append(((e, n), (d, n)))

    return keys

def encrypt(m, package):
    e, n = package
    c = pow(m, e, n) 
    return c

def decrypt(c, package):
    d, n = package
    m = pow(c, d, n) 
    return (m)

def modular_inverse(a,m): 
    a = a % m; 
    for x in range(1, m) : 
        if ((a * x) % m == 1) : 
            return x 
    return 1


def period(a,N):
    
    available_qubits = N.bit_length()
    r=-1   
    
    if N >= 2**available_qubits:
        print(str(N)+' is too big for IBMQX')
    
    qr = qiskit.QuantumRegister(available_qubits)   
    cr = qiskit.ClassicalRegister(available_qubits)
    qc = qiskit.QuantumCircuit(qr,cr)
    x0 = randint(1, N-1) 
    x_binary = np.zeros(available_qubits, dtype=bool)
    
    for i in range(1, available_qubits + 1):
        bit_state = (N%(2**i)!=0)
        if bit_state:
            N -= 2**(i-1)
        x_binary[available_qubits-i] = bit_state    
    
    for i in range(0,available_qubits):
        if x_binary[available_qubits-i-1]:
            qc.x(qr[i])
    x = x0
    
    while np.logical_or(x != x0, r <= 0):
        r+=1
        qc.measure(qr, cr) 
        for i in range(0,3): 
            qc.x(qr[i])
        qc.cx(qr[2],qr[1])
        qc.cx(qr[1],qr[2])
        qc.cx(qr[2],qr[1])
        qc.cx(qr[1],qr[0])
        qc.cx(qr[0],qr[1])
        qc.cx(qr[1],qr[0])
        qc.cx(qr[3],qr[0])
        qc.cx(qr[0],qr[1])
        qc.cx(qr[1],qr[0])
        
        result = qiskit.execute(qc,backend = qasm_sim, shots=1024).result()
        counts = result.get_counts()
        #print("qubits=",qc.num_qubits,"depth:",qc.depth())
        
        results = [[],[]]
        for key,value in counts.items(): 
            results[0].append(key)
            results[1].append(int(value))
        s = results[0][np.argmax(np.array(results[1]))]
    return r

def shors_breaker(N):
    N = int(N)
    while True:
        a=randint(0,N-1)
        g=gcd(a,N)
        if g!=1 or N==1:
            #print("not quantum!")
            return g,N//g
        else:
            #print("calling quantum subroutine...")
            r=period(a,N) 
            if r % 2 != 0:
                continue
            elif pow(a,r//2,N)==-1:
                continue
            else:
                p=gcd(pow(a,r//2)+1,N)
                q=gcd(pow(a,r//2)-1,N)
                if p==N or q==N:
                    continue
                return p,q
            


def shor_factorization(N):

    obj = Shor(quantum_instance = qasm_sim)
    res = obj.factor(N)
    #result = Shor_gen.run()
    print(res)
    p, q = res.factors[0]
    return (p, q)
############################################################

# list_n = (n=pq, bit_len(n), lambda=lcm(p-1,q-1)=ab/gcd(a,b))
maxprime=100
primes=[]
l=0
for prime_n in sympy.primerange(3, maxprime):
    primes.append(prime_n)

for p in primes:
    print(p)

list_n=[(primes[a]*primes[b], (primes[a]*primes[b]).bit_length(),int((primes[a]-1)*(primes[b]-1)/gcd(primes[a]-1,primes[b]-1))) for a in range(len(primes)) for b in range(a+1, len(primes))]

for n in list_n:
    print(n)

list_n.sort(key = lambda x: x[0])

# e>2 -> n>=10
#list_n=list_n[3:]

keys=generate_my_keys(list_n)

############################################################
file = open("qiskit_shor.dat", "w")
file.write("bitlen(n)\tn\tres\ttime\tqubits\tdepth\n")


candidate_N=[2*k+1 for k in range(1,60) if not sympy.isprime(2*k+1)]

print(candidate_N)
qasm_sim = qiskit.Aer.get_backend('qasm_simulator')

for N in candidate_N:

    start = time.time()
    obj = Shor(quantum_instance = qasm_sim)
    res = obj.factor(N)
    end = time.time()

    circ = Shor().construct_circuit(N)

    line=str(N.bit_length())+"\t"+str(N)+"\t"+str(res)+"\t"+str(end-start)+"\t"+str(circ.num_qubits)+"\t"+str(circ.depth())
    file.write(line+"\n")

file.close() 