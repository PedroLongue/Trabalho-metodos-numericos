import numpy as np
import matplotlib.pyplot as plot
from random import random

#Calcula matrix por Gauss Seidel
def matrixCalculation(A, B, initialSolution):
    Alength = len(A)                   
    
    for k in range(0, Alength):        
        aux = B[k]                  
          
        for i in range(0, Alength):     
            if(k != i):
                aux = aux - (A[k][i] * initialSolution[i])      
        initialSolution[k] = aux / A[k][k]      

    return initialSolution   

#Construção da Matriz
def createMatrix(constants, Node):
    DeltaX = constants[0]/(Node[1]+1)
    aux = np.arange(DeltaX,L,DeltaX)
    internalNodeCalculation = [0 for i in range(len(aux)+2)]

    internalNodeCalculation[0]=0
    internalNodeCalculation[len(aux)+1]=constants[0]
    for i in range(0,len(aux)):
        internalNodeCalculation[i+1]=(aux[i])

    A=[[0 for aux in range(Node[1])]for aux2 in range(Node[1])]
    for i in range(Node[1]):
        for j in range(Node[1]):
            if(i==j and j==0):
                A[i][j] = 2+ (DeltaX**2 * K/D)
                A[i][j+1]= -1
                continue

            elif(i==j and j==(Node[1]-1)):
                A[i][j] = 2+(DeltaX**2 * K/D)
                A[i][j-1]= -1
                continue

            elif(i==j):
                A[i][j] = 2+ (DeltaX**2 * K/D)
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

    return internalNodeCalculation, aux, concentration, x 

#Gráfico delta X x concentração
def renderGraph(constants, Node):
    hex = '#%06X' % round(random() * 0xffffff)
    concentration=createMatrix(constants, Node)
    plot.plot(concentration[3], concentration[0], color=hex)
    plot.xlabel("delta X")
    plot.ylabel("concentração")
    plot.show()

#Gráfico da sensibilidade
def sensitivity(Node, constants, steps):
    sensitivityD = [constants[2]-2*steps[1], constants[2]-steps[1], constants[2], constants[2]+steps[1], constants[2]+2*steps[1]]
    sensitivityK = [constants[1]-2*steps[0], constants[1]-steps[0], constants[1], constants[1]+steps[0], constants[1]+2*steps[0]]

    for i in range(0, len(sensitivityK)):
        hex = '#%06X' % round(random() * 0xffffff)
        descriptionD = createMatrix([constants[0], constants[1], sensitivityD[i]], [Node[0], Node[1], Node[2]])
        plot.plot(descriptionD[0], descriptionD[3], color=hex, label="sensibilidade D: "+str(sensitivityD[i]))
        plot.legend(title="D")
    plot.show()

    for i in range(0, len(sensitivityD)):
        hex = '#%06X' % round(random() * 0xffffff)
        descriptionK = createMatrix([constants[0], sensitivityK[i], constants[2]], [Node[0], Node[1], Node[2]])
        plot.plot(descriptionK[0], descriptionK[3], color=hex, label="sensibilidade K: "+str(sensitivityK[i]))
        plot.legend(title="K")
    plot.show()

#Gráfico do refinamento
def refinement(Node, constants, i):
    for k in range(i[0], i[1]+1):
        hex = '#%06X' % round(random() * 0xffffff)
        system = createMatrix([constants[0], constants[1], constants[2]], [Node[0], k, Node[1]])
        plot.plot(system[3], color = hex, label = 'Nos interiores: ' + str(k))
        plot.xlabel("delta X")
        plot.ylabel("Concentração")
        plot.legend(title="Refinamento")
    plot.show()

#Node
firstNode = 0.1
lastNode = 0
internalNode = 80

#convenção adotada
convention = 10**(-6)

#constants
L= 8
D = 4*convention #alfa
K = 8*convention #beta

Ii = 2
If = 6
deltaK = 0.5*convention
deltaD = 2*convention 

renderGraph([L, K, D], [firstNode, internalNode, lastNode])
refinement([firstNode, lastNode], [L, K, D], [Ii, If])
sensitivity([firstNode, internalNode, lastNode], [L, K, D], [deltaK, deltaD])

