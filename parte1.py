import numpy as np
import matplotlib.pyplot as plot

from random import random

#Realiza o calculo de RK3 com sistemas
def rungeKuttaSystem(equations, t, tMax, constant, deltaT):
    increment = np.arange(t, tMax, deltaT) #incrementa do t inical ate o t final usando o delta t

    c_a=[len(increment)]
    c_b=[len(increment)]
    c_c=[len(increment)]

    c_a[0]= equations[0]
    c_b[0]= equations[1]
    c_c[0]= equations[2]

    #Função de CA
    def derivateC_A(t,c_a,c_b,c_c):
        return ((-constant[0]*c_a*c_c)+c_b)

    #Função de CB
    def derivateC_B(t,c_a,c_b,c_c):
        return ((constant[1]*c_a*c_c)-c_b)

    #Função de CC
    def derivateC_C(t,c_a,c_b,c_c):
        return ((-constant[2]*c_a*c_c+c_b)-2*c_c)

    for j in range(t, len(increment)-1):
        k1Ca = derivateC_A(increment[j], c_a[j], c_b[j], c_c[j])
        k1Cb = derivateC_B(increment[j], c_a[j], c_b[j], c_c[j])
        k1Cc = derivateC_C(increment[j], c_a[j], c_b[j], c_c[j])

        k2Ca= derivateC_A(increment[j]+ deltaT/2, c_a[j] + k1Ca*deltaT/2, c_b[j] + k1Cb*deltaT/2, c_c[j] + k1Cc*deltaT/2 )
        k2Cb= derivateC_B(increment[j]+ deltaT/2, c_a[j] + k1Ca*deltaT/2, c_b[j] + k1Cb*deltaT/2, c_c[j] + k1Cc*deltaT/2)
        k2Cc= derivateC_C(increment[j]+ deltaT/2, c_a[j] + k1Ca*deltaT/2, c_b[j] + k1Cb*deltaT/2, c_c[j] + k1Cc*deltaT/2 )

        k3Ca= derivateC_A(increment[j]+ deltaT, c_a[j]-k1Ca*deltaT + 2*k2Ca*deltaT, c_b[j]-k1Cb*deltaT + 2*k2Cb*deltaT, c_c[j]-k1Cc*deltaT + 2*k2Cc*deltaT)
        k3Cb= derivateC_B(increment[j]+ deltaT, c_a[j]-k1Ca*deltaT + 2*k2Ca*deltaT, c_b[j]-k1Cb*deltaT + 2*k2Cb*deltaT, c_c[j]-k1Cc*deltaT + 2*k2Cc*deltaT)
        k3Cc= derivateC_C(increment[j]+ deltaT, c_a[j]-k1Ca*deltaT + 2*k2Ca*deltaT, c_b[j]-k1Cb*deltaT + 2*k2Cb*deltaT, c_c[j]-k1Cc*deltaT + 2*k2Cc*deltaT)


        c_a.append(c_a[j] + (1/6)*(k1Ca+4*k2Ca+k3Ca)*deltaT)
        c_b.append(c_b[j] + (1/6)*(k1Cb+4*k2Cb+k3Cb)*deltaT)
        c_c.append(c_c[j] + (1/6)*(k1Cc+4*k2Cc+k3Cc)*deltaT)

    # print('k1Ca: ', k1Ca)
    # print('k1Cb: ', k1Cb)
    # print('k1Cc: ', k1Cc)
    # print('\n')
    # print('k2Ca: ', k2Ca)
    # print('k2Cb: ', k2Cb)
    # print('k2Cc: ', k2Cc)
    # print('\n')
    # print('k3Ca: ', k3Ca)
    # print('k3Cb: ', k3Cb)
    # print('k3Cc: ', k3Cc)

    return c_a, c_b, c_c, increment
 
#Gráfico tempo x C
def renderGraph(equations, t, tMax, constant, deltaT):
    hex = '#%06X' % round(random() * 0xffffff)
    hex1 = '#%06X' % round(random() * 0xffffff)
    hex2 = '#%06X' % round(random() * 0xffffff)
    resultRk = rungeKuttaSystem(equations,t,tMax,constant,deltaT)

    plot.plot(resultRk[3],resultRk[0],color=hex, label='ca')
    plot.plot(resultRk[3],resultRk[1],color=hex1, label='cb')
    plot.plot(resultRk[3],resultRk[2],color=hex2, label='cc')
    plot.legend(title="C")
    plot.xlabel("t(s)")
    plot.ylabel("C")
    plot.show()

    print("O aproximado de Ca=", resultRk[0][len(resultRk[0])-1])
    print("O aproximado de Cb=", resultRk[1][len(resultRk[1])-1])
    print("O aproximado de Cc=", resultRk[2][len(resultRk[2])-1])

#Gráfico refinamento
def refinamento(equations, t, tMax, constant, deltaT0, deltaT1, deltaT):
    DeltaT = np.arange(deltaT0, deltaT1, 0.1)
    x, (ca, cb, cc)= plot.subplots(1,3)
    for i in range(0,len(DeltaT)):
        hex = '#%06X' % round(random() * 0xffffff)
        rungeKuttaSystem1 = rungeKuttaSystem(equations, t, tMax, constant, DeltaT[i])
        print(rungeKuttaSystem1[0],rungeKuttaSystem1[3])
        ca.plot(rungeKuttaSystem1[3],rungeKuttaSystem1[0],color=hex, label='%.2f'%(float(DeltaT[i])))
        ca.legend(title="delta T")
        ca.set_xlabel("t(s)")
        ca.set_ylabel("ca")
        cb.plot(rungeKuttaSystem1[3],rungeKuttaSystem1[1],color=hex, label='%.2f'%(float(DeltaT[i])))
        cb.legend(title="delta T")
        cb.set_xlabel("t(s)")
        cb.set_ylabel("cb")
        cc.plot(rungeKuttaSystem1[3],rungeKuttaSystem1[2],color=hex, label='%.2f'%(float(DeltaT[i])))
        cc.legend(title="delta T")
        cc.set_xlabel("t(s)")
        cc.set_ylabel("cc") 
    plot.show()  

#Gráfico sensibilidade
def sensibilidade(equations, t, tMax, constant, deltaT):
    sensitivityAlpha = [constant[0]-0.1*constant[0],constant[0]-0.05*constant[0], constant[0], constant[0]+0.05*constant[0], constant[0]+0.10*constant[0]]
    sensitivityBeta = [constant[1]-0.1*constant[1],constant[1]-0.05*constant[1], constant[1], constant[1]+0.05*constant[1], constant[1]+0.10*constant[1]]
    sensitivityGamma = [constant[2]-0.1*constant[2],constant[2]-0.05*constant[2], constant[2], constant[2]+0.05*constant[2], constant[2]+0.10*constant[2]]
    x, (alphaGraph, betaGraph, gammaGraph)= plot.subplots(1,3)

    for i in range(0, len(sensitivityAlpha)):
        hex = '#%06X' % round(random() * 0xffffff)
        rungeKuttaAlpha = rungeKuttaSystem(equations, t, tMax, constant, deltaT)
        rungeKuttaBeta = rungeKuttaSystem(equations, t, tMax, constant, deltaT)
        rungeKuttaGamma = rungeKuttaSystem(equations, t, tMax, constant, deltaT)

        alphaGraph.plot(rungeKuttaAlpha[3],rungeKuttaAlpha[0],color=hex, label='%.2f'%(float(sensitivityAlpha[i])))
        alphaGraph.legend(title="alfa")
        alphaGraph.set_xlabel("t(s)")
        alphaGraph.set_ylabel("ca")
        betaGraph.plot(rungeKuttaBeta[3],rungeKuttaBeta[1],color=hex, label='%.2f'%(float(sensitivityBeta[i])))
        betaGraph.legend(title="beta")
        betaGraph.set_xlabel("t(s)")
        betaGraph.set_ylabel("cb")
        gammaGraph.plot(rungeKuttaGamma[3],rungeKuttaGamma[2],color=hex, label='%.2f'%(float(sensitivityGamma[i])))
        gammaGraph.legend(title="gamma")
        gammaGraph.set_xlabel("t(s)")
        gammaGraph.set_ylabel("cc") 
    plot.show()

ca, cb, cc = 10, 0, 15
t0 = 0
tMax = 3
deltaT = 1.1
alfa, beta, gamma = 10, 15, 17

#--------------REFINAMENTO--------------#
deltaT0 = 0.0001
deltaT1 = 0.0002

rungeKuttaSystem([ca, cb, cc],t0,tMax,[alfa, beta, gamma],deltaT)

# print("Ca: ", rungeKuttaSystem([ca, cb, cc],t0,tMax,[alfa, beta, gamma],deltaT)[0])
# print("Cb: ", rungeKuttaSystem([ca, cb, cc],t0,tMax,[alfa, beta, gamma],deltaT)[1])
# print("Cc: ", rungeKuttaSystem([ca, cb, cc],t0,tMax,[alfa, beta, gamma],deltaT)[2])

renderGraph([ca, cb, cc],t0,tMax,[alfa, beta, gamma],deltaT)
refinamento([ca,cb,cc],t0, tMax, [alfa,beta, gamma],deltaT0, deltaT1,deltaT)
sensibilidade([ca,cb,cc], t0, tMax, [alfa,beta, gamma], deltaT)

