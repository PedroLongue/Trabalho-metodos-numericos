import numpy as np
import matplotlib.pyplot as plt

def rungeKutta3System(equations, t, tMax, constant, deltaT):
    increment = np.arange(t, tMax, deltaT) #incrementa do t inical ate o t final usando o delta t

    c_a=[len(increment)]
    c_b=[len(increment)]
    c_c=[len(increment)]

    #condições iniciais
    c_a[0]= equations[0]
    c_b[0]= equations[1]
    c_c[0]= equations[2]

    def derivateC_A(t,c_a,c_b,c_c):
        return ((-constant[0]*c_a*c_c)+c_b)

    def derivateC_B(t,c_a,c_b,c_c):
        return ((constant[1]*c_a*c_c)-c_b)

    def derivateC_C(t,c_a,c_b,c_c):
        return ((-constant[2]*c_a*c_c+c_b)-2*c_c)

    for i in range(len(equations)):
        print(f'f{i+1}({t}) = {equations[i]},  {t=},  {deltaT=}')

    print('\n')

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
 
def renderGraph(equations, t, tMax, constant, deltaT):
    # print("aaa", len(resultRk[3]))
    resultRk = rungeKutta3System(equations,t,tMax,constant,deltaT)
    plt.plot(resultRk[3],resultRk[0],color='purple', label='ca')
    plt.plot(resultRk[3],resultRk[1],color='black', label='cb')
    plt.plot(resultRk[3],resultRk[2],color='blue', label='cc')
    plt.legend(title="C",loc=0)
    plt.xlabel("t(s)")
    plt.ylabel("C")
    plt.show()
    print("a=", resultRk[0][len(resultRk[0])-1])
    print("b=", resultRk[1][len(resultRk[1])-1])
    print("c=", resultRk[2][len(resultRk[2])-1])

ca, cb, cc = 1, 0, 1
t0 = 0
tMax = 1.5
deltaT = 0.5
alfa, beta, gama = 1, 1, 1

# rungeKutta3System([ca, cb, cc],t0,tMax,[alfa, beta, gama],deltaT)

# print("Ca: ", rungeKutta3System([ca, cb, cc],t0,tMax,[alfa, beta, gama],deltaT)[0])
# print("Cb: ", rungeKutta3System([ca, cb, cc],t0,tMax,[alfa, beta, gama],deltaT)[1])
# print("Cc: ", rungeKutta3System([ca, cb, cc],t0,tMax,[alfa, beta, gama],deltaT)[2])

renderGraph([ca, cb, cc],t0,tMax,[alfa, beta, gama],deltaT)