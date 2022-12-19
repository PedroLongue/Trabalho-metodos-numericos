import numpy as np
import matplotlib.pyplot as plot
from random import random

def createMatrix(constants, Node):
    DeltaX = constants[0]/(Node[1]+1)
    s = DeltaX**2 * K/D
    aux = np.arange(DeltaX,L-DeltaX,DeltaX)
    internalNodeCalculation = [0 for i in range(len(aux)+3)]

    internalNodeCalculation[0]=0
    internalNodeCalculation[(len(internalNodeCalculation)-1)]=internalNodeCalculation[len(internalNodeCalculation)-5]+DeltaX
    for i in range(0,len(aux)):
        internalNodeCalculation[i+1]=(aux[i])

    A=[[0 for aux in range(Node[1])]for aux2 in range(Node[1])]
    for i in range(Node[1]):
        for j in range(Node[1]):
            if(i==j and j==0):
                A[i][j] = 2+s
                A[i][j+1]= -1
                continue

            elif(i==j and j==(Node[1]-1)):
                A[i][j] = 2+s
                A[i][j-1]= -1
                continue

            elif(i==j):
                A[i][j] = 2+s
                A[i][j+1]= -1
                A[i][j-1]= -1
        
    B=[0 for aux in range(Node[1])]
    B[0] = Node[0]
    B[Node[1]-1] = Node[2]
                        
    initialSolution = [0.5 for i in range(len(A))]                        

    for i in range(0, 3000):            
        concentration = matrixCalculation(A, B, initialSolution)

    x=[]
    for i in range(Node[1]+2):
        x.append(0)
    x[0]=Node[0]
    x[Node[1]+1]=Node[2]
    for i in range(0,Node[1]):
        x[i+1]=(concentration[i])
       
    for i in range(len(x)):
        print('concentação '+str(i+1)+" = ",x[i])
    print("\n")

    return internalNodeCalculation, aux, concentration, x 

def matrixCalculation(a, b, initialSolution):
    n = len(a)                   
    
    for j in range(0, n):        
        aux = b[j]                  
          
        for i in range(0, n):     
            if(j != i):
                aux = aux - (a[j][i] * initialSolution[i])      
        initialSolution[j] = aux / a[j][j]      

    return initialSolution   

def renderGraph(constants, Node):
    hex = '#%06X' % round(random() * 0xffffff)
    concentration=createMatrix(constants, Node)
    plot.plot(concentration[3], concentration[0], color=hex)
    plot.xlabel("delta X")
    plot.ylabel("concentração")
    plot.show()

def Refinamento(Node, constants, i):
    for k in range(i[0], i[1]+1):
        hex = '#%06X' % round(random() * 0xffffff)
        system = createMatrix([constants[0], constants[1], constants[2]], [Node[0], k, Node[1]])
        plot.plot(system[3], color = hex, label = 'Nos Iternos=' + str(k))
        plot.legend(title='Numero de Nos Internos')
        plot.xlabel("Delta X")
        plot.ylabel("Concentração")

    plot.show()

def Sensibilidade(Node, constants, passoK, passoD):
    d = [D-2*passoD, D-passoD, D, D+passoD, D+2*passoD]
    k = [K-2*passoK, K-passoK, K, K+passoK, K+2*passoK]
    for i in range(0, len(d)):
        hex = '#%06X' % round(random() * 0xffffff)
        SolverK = createMatrix([constants[0], k[i], constants[2]], [Node[0], Node[1], Node[2]])
        plot.plot(SolverK[1], SolverK[3], color=hex, label="K="+str(k[i]))
        plot.legend(title="K")
    plot.show()
    for i in range(0, len(k)-1):
        hex = '#%06X' % round(random() * 0xffffff)
        SolverD = createMatrix([constants[0], constants[1], d[i]], [Node[0], Node[1], Node[2]])
        plot.plot(SolverD[1], SolverD[3], color=hex, label="D="+str(d[i]))
        plot.legend(title="D")
    plot.show()

firstNode = 0.1
lastNode = 0
internalNode = 20

convention = 10**(-6)

L= 4
K = 8*convention
D = 2*convention

Ii = 3
If = 6
deltaK = 1*convention
deltaD = 0.5*convention 

renderGraph([L, K, D], [firstNode, internalNode, lastNode])
Refinamento([firstNode, lastNode], [L, K, D], [Ii, If])
Sensibilidade([firstNode, internalNode, lastNode], [L, K, D], deltaK, deltaD)