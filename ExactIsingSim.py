# Import the Qiskit
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, QiskitError
from qiskit import execute, IBMQ, BasicAer
from qiskit.providers.ibmq import least_busy

#import Qconfig  
# useful additional packages
import matplotlib.pyplot as plt

import numpy as np
from scipy import linalg as la

# Authenticate for access to remote backends
try:
    IBMQ.load_accounts()
except:
    print("""WARNING: There's no connection with the API for remote backends.
             Have you initialized a file with your personal token?
             For now, there's only access to local simulator backends...""")
def had(qc,q):
    qc.h(q)
    
def digit_sum(n):
    num_str = str(n)
    sum = 0
    for i in range(0, len(num_str)):
        sum += int(num_str[i])
    return sum

# CZ (Controlled-Z)
# control qubit: q0
# target qubit: q1
def CZ(qp,q0,q1):
    qp.h(q1)
    qp.cx(q0,q1)
    qp.h(q1)
# f-SWAP
# taking into account the one-directionality of CNOT gates in the available devices
def fSWAP(qp,q0,q1):
    qp.cx(q0,q1)
    qp.h(q0)
    qp.h(q1)
    qp.cx(q0,q1)
    qp.h(q0)
    qp.h(q1)
    qp.cx(q0,q1)
    CZ(qp,q0,q1)
    
# CH (Controlled-Haddamard)
# control qubit: q1
# target qubit: q0
def CH2(qp,q0,q1):
    qp.sdg(q0)
    qp.h(q0)
    qp.tdg(q0)
    qp.h(q0)
    qp.h(q1)
    qp.cx(q0,q1)
    qp.h(q0)
    qp.h(q1)
    qp.t(q0)
    qp.h(q0)
    qp.s(q0)
# Fourier transform gates
def F2(qp,q0,q1):
    qp.cx(q0,q1)
    CH2(qp,q0,q1)
    qp.cx(q0,q1)
    CZ(qp,q0,q1) 
def F0(qp,q0,q1):
    F2(qp,q0,q1)   
def F1(qp,q0,q1):
    F2(qp,q0,q1)
    qp.sdg(q0)
    
from math import pi

# ROTATIONAL GATES
def RZ(qp,th,q0):
    qp.u1(-th,q0)
def RY(qp,th,q0):
    qp.u3(th,0.,0.,q0)
def RX(qp,th,q0):
    qp.u3(th,0.,pi,q0)

# CRX (Controlled-RX)
# control qubit: q0
# target qubit: q1
def CRX(qp,th,q0,q1):
    RZ(qp,pi/2.0,q1)
    RY(qp,th/2.0,q1)
    qp.cx(q0,q1)
    RY(qp,-th/2.0,q1)
    qp.cx(q0,q1)
    RZ(qp,-pi/2.0,q1)
# Bogoliubov B_1
def B(qp,thk,q0,q1):
    qp.x(q1)
    qp.cx(q1,q0)
    CRX(qp,thk,q0,q1)
    qp.cx(q1,q0)
    qp.x(q1)
    
# This circuit can be implemented in ibmqx5 using qubits (q0,q1,q2,q3)=(6,7,11,10)
# It can also be implemented between other qubits or in ibqmx2 and ibqmx4 using fermionic SWAPS
# For instance, the lines commented correspond to the implementations:
# ibmqx2 (q0,q1,q2,q3)=(4,2,0,1)
# ibmqx4 (q0,q1,q2,q3)=(3,2,1,0)
def Udisg(qc,lam,q0,q1,q2,q3):
    k=1
    n=4
    th1=-np.arccos((lam-np.cos(2*pi*k/n))/np.sqrt((lam-np.cos(2*pi*k/n))**2+np.sin(2*pi*k/n)**2))
    B(Udis,th1,q0,q1)
    F1(Udis,q0,q1)
    F0(Udis,q2,q3)
    #fSWAP(Udis,q2,q1) # for ibmqx2
    #fSWAP(Udis,q1,q2) # for ibmqx4
    F0(Udis,q0,q2)
    F0(Udis,q1,q3)
    #fSWAP(Udis,q2,q1) # for ibmqx2
    #fSWAP(Udis,q1,q2) # for ibmqx4


def Initial(qc,lam,q0,q1,q2,q3):
    if lam <1:
        qc.x(q3)
def Ising(ini,udis,mes,lam,q0,q1,q2,q3,c0,c1,c2,c3):
    Initial(ini,lam,q0,q1,q2,q3)
    Udisg(udis,lam,q0,q1,q2,q3)
    mes.measure(q0,c0)
    mes.measure(q1,c1)
    mes.measure(q2,c2)
    mes.measure(q3,c3)
    return ini+udis+mes

try:
    IBMQ.load_accounts()
except:
    print("""WARNING: There's no connection with the API for remote backends.
             Have you initialized a file with your personal token?
             For now, there's only access to local simulator backends...""")
print("IBMQ backends: ", IBMQ.backends())
print("IBMQ backends: ", BasicAer.backends())
backend_sim = BasicAer.get_backend('qasm_simulator')
print("Connected to" + str(backend_sim))
shots = 1024
coupling_map = None
mag_sim = []
for i in range(8):
    q = QuantumRegister(4,"q")
    c = ClassicalRegister(4,"c")
    Udis = QuantumCircuit(q,c)
    ini = QuantumCircuit(q,c)
    mes = QuantumCircuit(q,c)

    lam=i*0.25
    sumter = Ising(ini,Udis,mes,lam,q[0],q[1],q[2],q[3],c[0],c[1],c[2],c[3])

    #Isex.set_api(Qconfig.APItoken, Qconfig.config["url"]) # set the APIToken and API url   

    result = execute(sumter, backend=backend_sim,
                    coupling_map=coupling_map, shots=shots).result()
    res=result.get_counts(sumter)
    r1=list(res.keys())
    r2=list(res.values())
    M=0
    for j in range(0,len(r1)):
        M=M+(4-2*digit_sum(r1[j]))*r2[j]/shots
    #print("$\lambda$: ",lam,", $<\sigma_{z}>$: ",M/4)
    mag_sim.append(M/4)


def exact(lam):
    if lam <1:
        return lam/(2*np.sqrt(1+lam**2))
    if lam >1:
        return 1/2+lam/(2*np.sqrt(1+lam**2))
    return None

vexact = np.vectorize(exact)
l=np.arange(0.0,2.0,0.01)
l1=np.arange(0.0,2.0,0.25)
plt.figure(figsize=(9,5))
#plt.plot(l,vexact(l),'k',label='exact')
plt.plot(l1, mag_sim, 'bo')
plt.xlabel('$g$')
plt.ylabel('$<\sigma_{x}>$')
plt.legend()
plt.title('Transverse magnetization of the ground state of n=4 Ising spin chain')
plt.show()
