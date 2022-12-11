import matplotlib.pyplot as plt
import scipy
from scipy.integrate import quad
import numpy as np

def L_prac(L,di,n,d,f):
    if(n/d == 1/2):
        Lprac = 0.48*((L/di)/(1+L/di))*(3*10**8/f)
    elif(n/d == 1):
        Lprac = 0.96*((L/di)/(1+L/di))*(3*10**8/f)
    elif(n/d == 3/2):
        Lprac = 1.44*((L/di)/(1+L/di))*(3*10**8/f)
    elif(n/d == 2):
        Lprac = 1.92*((L/di)/(1+L/di))*(3*10**8/f)
    return Lprac
def E(x,Im,r,k,H):
    n = 120*np.pi
    return ((n*Im/(2*np.pi*r))*((np.cos(k*H*np.cos(x))-np.cos(k*H))/np.sin(x)))**2*r**2*np.sin(x)/n
def I2(x):
    return 1

def potencia_radiada(Im,r,f,L):        #L=2H
    k = 2*np.pi/(3*10**8/f)            #k = numero de onda
    H=L/2
    P1,err1 = quad(E, 0, np.pi, args=(Im,r,k,H))  #x = tetha
    P2,err2 = quad(I2, 0, 2*np.pi)                #integral de phi de 0-2pi
    print("\t\tEn watts")
    print("\t\t{0:.5f}".format(P2*P1),"W")          #se restringe a una presición de 5 decimales
    if(P2*P1 < 1):
        print("\t\tEn miliwatts")
        print("\t\t{0:.5f}".format(P2*P1*1000),"mW")
    return P1*P2
def R_r(P,Im):
    Rr = 2*P/Im
    print("\t\t{0:.5f}".format(Rr),"ohms")
    return Rr
def R_p(f,L, a, o):	#L: longitud de la antena, a: radio de la antena, o: conductancia del material
    uo = 4*np.pi*10**(-7) 		#uo permeabilidad magnetica
    Rs = np.sqrt((np.pi*f*uo)/o)		
    Rp = Rs*L/(2*np.pi*a)
    print("\t\t{0:.5f}".format(Rp), "ohms")
    return Rp

def reactancia(L,a,f):
    #H = L/2
    #Zo = 120*(np.log(2*H/a)-1)
    #k = 2*np.pi/(3*10**8/f)            #k = numero de onda
    #X = -Zo*(1/np.tan(H*k))
    X = 0
    print("\t\t{0:.1f}".format(X),"j ohms")
    return X

def e_cd(Rr, Rp):
    ecd = Rr/(Rr+Rp)
    print("\t\t{0:.5f}".format(ecd))
    return ecd

def max_E(n,Im,r,k,H):
    maxim = 0.0
    aux = 0.0
    #como python no deja incrementar de a 0.01 el for, se incrementa de a 1 y la x se divide en 100 en los calculos, 2pi=6.3 y se multiplica por 100 para la equivalencia
    for x in range(1,630,1):
        x1 = x/100 #Variable auxiliar por que python arroja error al dividir dentro del coseno
        aux = (np.cos(k*H*np.cos(x1))-np.cos(k*H))/np.sin(x1)
        if(maxim < aux):
            maxim = aux
    maxim = maxim*(n*Im/(2*np.pi*r))
    return maxim

def En(x,Im,r,k,H):
    n = 120*np.pi
    Emax = max_E(n,Im,r,k,H)
    return (((n*Im/(2*np.pi*r))*((np.cos(k*H*np.cos(x))-np.cos(k*H))/np.sin(x)))/Emax)*np.sin(x)/n

def ancho_haz(Im,r,f,L):
    k = 2*np.pi/(3*10**8/f)            #k = numero de onda
    H=L/2
    A1,err1 = quad(En, 0, np.pi, args=(Im,r,k,H))  #x = tetha
    A2,err2 = quad(I2, 0, 2*np.pi)                #integral de phi de 0-2pi
    ancho = np.rad2deg(abs(A2*A1*180/np.pi))
    print("\t\t{0:.5f}".format(ancho))
    return ancho

def directividad(ancho):
    D = 4*np.pi/ancho
    print("\t\t{0:.5f}".format(D))
    return D

def apertura_Efec(D, lambd):
    Ap_Efec = (D*lambd**2)/(4*np.pi)
    print("\t\t{0:.5f}".format(Ap_Efec))
    return Ap_Efec

def gain_antena(e_cd, D):
    G = e_cd*D
    print("\t\t{0:.5f}".format(G),"veces")
    print("\t\t{0:.5f}".format(10*np.log10(G)),"dBi")
    return G

def dipolo():
    print("\n\n════════════════════════════════════════════════════════════════════════════")
    print("Este trabajo calcula los parametros de la antena dipolo de λ/2, λ, 3λ/2 o 2λ")
    print("dados los siguientes parametros:\n")
    print("════════════════════════════════════════════════════════════════════════════")
    print("\n>>>Ingrese la frecuencia de trabajo [Hz]:")
    fo = float(input())
    while(fo < 0):
        print("\n═════════════════════════════════")
        print("Dato incorrecto, ingrese de nuevo")
        print("═════════════════════════════════")
        print(">>>Ingrese la frecuencia de trabajo [Hz]:")
        fo = float(input())
        
    print("\n>>>Ingrese la Corriente máxima [A]:")
    Im = float(input())
    while(Im < 0):
        print("\n═════════════════════════════════")
        print("Dato incorrecto, ingrese de nuevo")
        print("═════════════════════════════════")
        print(">>>Ingrese la Corriente máxima [A]:")
        Im = float(input())
    print("\n>>>Ingrese la longitud de la antena, en terminos de (n/d)*lambda:")
    n = int(input("n: "))
    d = int(input("d: "))
    if d == 0:
	    L = 0
    else:
	    L = n/d
    #se verifica que no sea una antena pequeña, que el denominador no sea cero y que este en las longitudes establecidas
    while ((L <= 1/10) or (d == 0) or ((n/d != 1/2) and (n/d != 3/2) and (n/d != 1) and (n/d != 2))):
            print("\n═════════════════════════════════")
            print("Dato incorrecto, ingrese de nuevo")
            print("═════════════════════════════════")
            n = int(input("n: "))
            d = int(input("d: "))
            if (d != 0):
                L = n/d
    L =(3*10**8/fo)*L
    ra = 1
    print("\n>>>Ingrese el diámetro de la antena:")
    di = float(input())
    while(di < 0):
        print("\n═════════════════════════════════")
        print("Dato incorrecto, ingrese de nuevo")
        print("═════════════════════════════════")
        print(">>>Ingrese el diámetro de la antena:")
        di = float(input())
    print("\n>>>Ingrese la conductividad del material de la antena:")
    o = float(input())
    while(di < 0):
        print("\n═════════════════════════════════")
        print("Dato incorrecto, ingrese de nuevo")
        print("═════════════════════════════════")
        print(">>>Ingrese la conductividad del material de la antena:")
        o = float(input())
    #calculo de parámetros
    print("\n\n═════════════════════════════════════")
    print(" Parametros de antena dipolo de", n, end = "")
    print("λ/", end = "")
    print(d)
    print("═════════════════════════════════════")
    print("\n\t>>>Lambda (λ):")
    print("\t\t",3*10**8/fo, "m")
    print("\n\t>>>L nominal:")
    print("\t\t",L, "m")
    print("\n\t>>>L práctica:")
    print("\t\t{0:.5f}".format(L_prac(L,di,n,d,fo)), "m")
    print("\n\t>>>Potencia Radiada")
    P = potencia_radiada(Im,ra,fo,L)
    print("\n\t>>>Resistencia Radiada")
    Rr = R_r(P,Im)
    print("\n\t>>>Resistencia de Pérdidas")
    a = di/2
    Rp = R_p(fo,L,a,o)
    print("\n\t>>>Reactancia")
    X = reactancia(L,a,fo)
    print("\n\t>>>Impedancia")
    if X >= 0: 
	    print("\t\t{0:.5f}".format(Rp+Rr),"+","{0:.1f}".format(X), "j ohms")
    else:
	    print("\t\t{0:.5f}".format(Rp+Rr),"-","{0:.1f}".format(-X), "j ohms")
    print("\n\t>>>Eficiencia")
    ecd = e_cd(Rr, Rp)
    print("\n\t>>>Ancho del haz")
    ancho = ancho_haz(Im,ra,fo,L)
    print("\n\t>>>Directividad")
    D = directividad(ancho)
    print("\n\t>>>Apertura efectiva")
    apertura_Efec(D, 3*10**8/fo)
    print("\n\t>>>Ganancia")
    G = gain_antena(ecd, D)
    print("\n\t>>>Patrón de radiación")
    n = 120*np.pi
    k = 2*np.pi/(3*10**8/fo)            #k = numero de onda
    H=L/2
    Emax = max_E(n,Im,ra,k,H)
    x = np.linspace(0,2*np.pi,1000)  #x = tetha
    r = 20*np.log10(abs((n*Im/(2*np.pi*ra))*((np.cos(k*H*np.cos(x))-np.cos(k*H))/(np.sin(x)*Emax))))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    ax.plot(x,r)
    #x = np.arange(0.001,2*np.pi-0.001,0.1)
    #y = 10*np.log10((n*Im/(2*np.pi*r))*((np.cos(k*H*np.cos(x))-np.cos(k*H))/(np.sin(x)*Emax)))
    #plt.plot(x,y)
    #plt.xlabel('Degrees [°]')
    #plt.ylabel('P [W]')
    #plt.title('Patrón de Radiación')
    plt.show()

def guia_onda():
    print("\n\n═════════════════════════════════════════════════════════════════════════════")
    print("Para la guía de onda rectangular parcialmente cargada que se muestra en la")
    print("figura adjunta, con ß=0, encontrar los parámetros de la guía en el modo TE10,")
    print("para Z>0 y Z<0. Suponga que a = 2,286 cm, b = a/2, Er = 2,25 y la frecuencia")
    print("de operación es de 10% más a la de corte.")
    print("═════════════════════════════════════════════════════════════════════════════")
    print("\t\t_____________________________________________")
    print("\t\t                     |")
    print("\t\t                     |")
    print("\t\t                     |")
    print("\t\t          E1         |           Eo          ")
    print("\t\t                     |")
    print("\t\t                     |")
    print("\t\t                     |")
    print("\t\t¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯")
    print("\t\t                    z=0")
    print("\nSe tiene que: Er1=2.25, Er2=0\n")
    Er=2.25
    a = 2.286*10**-2
    b = a/2
    m = 1
    n = 0
    print("\t>>>Frecuencia de corte:")
    fc1 = ((3*10**8)/(2*np.sqrt(Er)))*np.sqrt((m/a)**2+(n/b)**2)
    print("\t\t{0:.5f}".format(fc1/1000000000), "GHz            Z<0")
    fc2 = ((3*10**8)/2)*np.sqrt((m/a)**2+(n/b)**2)
    print("\t\t{0:.5f}".format(fc2/1000000000), "GHz            Z>0")
    print("\t>>>Longitud de onda de corte")
    lam1 = 3*10**8/fc1
    lam2 = 3*10**8/fc2
    print("\t\t{0:.5f}".format(lam1), "m              Z<0")
    print("\t\t{0:.5f}".format(lam2), "m              Z>0")
    print("\t>>>Impedancia modo TE10")
    n1 = np.sqrt((4*np.pi*10**-7)/(2.25*8.85*10**-12))
    n2 = np.sqrt((4*np.pi*10**-7)/(8.85*10**-12))
    z1 = n1/np.sqrt(1-(1/1.1)**2)
    z2 = n2/np.sqrt(1-(1/1.1)**2)
    print("\t\t{0:.5f}".format(z1), "ohms         Z<0")
    print("\t\t{0:.5f}".format(z2), "ohms         Z>0")
    print("\t>>>Velocidad de fase")
    vfas1 = (3*10**8/(1.1*fc1))*(1.1*fc1/np.sqrt((1.1*fc1)**2-fc1**2))
    vfas2 = (3*10**8/(1.1*fc2))*(1.1*fc2/np.sqrt((1.1*fc2)**2-fc2**2))
    print("\t\t{0:.5f}".format(vfas1), "m/s            Z<0")
    print("\t\t{0:.5f}".format(vfas2), "m/s            Z>0")
    print("\t>>>Velocidad de grupo")
    vgrup = 3*10**8*np.sqrt(1-(1/1.1)**2)
    #vgrup2 = 3*10**8*np.sqrt(1-(1/1.1)**2) tienen la misma vel de grupo, no depende de Er
    print("\t\t{0:.5f}".format(vgrup/1000), "Km/s      Z<0")

def circuit_micro():
    print("\n\n═════════════════════════════════════════════════════════════════════════════")
    print("\t\t           _______________________           ")
    print("\t\t          │->a1               b3->│          ")
    print("\t\t     P1═══│                       │═══P3     ")
    print("\t\t          │<-b1               a3<-│          ")
    print("\t\t          │                       │          ")
    print("\t\t          │->a2               b4->│          ")
    print("\t\t     P2═══│                       │═══P4     ")
    print("\t\t          │<-b2               a4<-│          ")
    print("\t\t           ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯           ")
    print("Una red de cuatro puertos tiene la matriz de dispersión que se")
    print("muestra a continuación:")
    magni = [[0.178,0.6,0.4,0],[0.6,0,0,0.3],[0.4,0,0,0.5],[0,0.3,0.5,0]]
    fase = [[90,45,45,0],[45,0,0,-45],[45,0,0,-45],[0,-45,-45,0]]
    print("\t\t┌─                                         ─┐")
    print("\t\t│  0.178└90°  0.6└45°   0.4└45°     0       │")
    print("\t\t│   0.6└45°      0         0     0.3└-45°   │")
    print("\t\t│   0.4└45°      0         0     0.5└-45°   │")
    print("\t\t│      0      0.3└-45°  0.5└-45°    0       │")
    print("\t\t└─                                         ─┘")
    print("a) ¿Es esta red sin pérdidas?")
    print("   Parq no tener pérdidas, las matriz de dispersión [S] debe ser unitaria,")
    print("   por ende, de la primera fila se debe cumplir:")
    sumat = 0
    for i in range(0,4):
        sumat = sumat + magni[0][i]**2
    print("   |s11|^2+|s12|^2+|s21|^2+|s22|^2=(0.178)^2+(0.6)^2+(0.4)^2+(0)^2 ={0:.3f}".format(sumat))
    if sumat == 0:
        print("   Por ende al ser igual a cero si tiene pérdidas")
    else:
        print("   Por ende al ser diferente de cero si tiene pérdidas\n")
    print("b) ¿Esta red es recíproca?")
    cont = 0
    for i in range(0,4):
        for j in range(0,4):
            if ((magni[i][j]==magni[j][i])and(fase[i][j]==fase[j][i])):
                cont+=1
    if cont==16:     #igual a 16 ya que son 16 terminos de la matriz
        print("   Al ser la matriz simétrica, se tiene que la red es recíproca\n")
    else:
        print("   Al ser no simétrica, se tiene que la red no es recíproca\n")
        
    print("c) ¿Cuál es la pérdida de retorno en el puerto 1 cuando todos los demás")
    print("   puertos terminan con cargas coincidentes?")
    print("   Como el coeficiente de reflección de cada puerto es cada s de la diagonal")
    print("   se tiene que:")
    RL = -20*np.log10(magni[0][0])
    print("   Pérdidas RL=-20log|s11|={0:.3f}".format(RL),"dB\n")
    
    print("d) ¿Cuál es la pérdida de inserción entre los puertos 2 y 4 cuando todos")
    print("   los demás puertos terminan con cargas coincidentes?")
    print("   Cuando los puertos 1 y 3 terminan en Zo, se tiene:")
    print("   V1+=0, V2+=0, por lo que V4-=s42*V2+")
    print("   Como los demás términos de la matriz diferentes a la diagonal son")
    print("   coeficientes de transmisión, las perdidas por inserción son:")
    IL = -20*np.log10(magni[3][1])
    print("   Pérdidas IL=-20log|s42|={0:.3f}".format(IL),"dB\n")

def main():
    flag = 0
    print("\n\n═════════════════════════════════════════════════════════════════════════════")
    print("\t\tMateo Felipe Umaña Gordillo, METX 2022-2")
    print("\tParámetros de antena dipolo, guía de onda y circuito de microondas")
    print("═════════════════════════════════════════════════════════════════════════════")
    while flag == 0:
        print("\nEscoja alguna de las siguientes opciones.\n")
        print("\t1)Diseño de antena de longitud n lambda.")
        print("\t2)Guía de onda.")
        print("\t3)Circuitos de microondas.\n")
        flag = int(input("Ingrese su opción: "))
        while(flag != 1 and flag != 2 and flag != 3):
            print("\n═════════════════════════════════")
            print("Dato incorrecto, ingrese de nuevo")
            print("═════════════════════════════════")
            flag = int(input("Ingrese su opción: "))
        if(flag == 1):
            dipolo()
        elif(flag == 2):
            guia_onda()
        elif(flag == 3):
            circuit_micro()
        print("\n\n════════════════════════════════════")
        print("¿Desea escoger otra opción del menú?")
        print("════════════════════════════════════")
        print("\t    Si: 4, No: 5")
        flag = int(input())
        while(flag!= 4 and flag != 5):
            print("\n═════════════════════════════════")
            print("Dato incorrecto, ingrese de nuevo")
            print("═════════════════════════════════")
            print("\t    Si: 4, No: 5")
            flag = int(input())
        print("\n\n")
        if flag==4:
            flag = 0    #Se hace flag=0 para entrar de nuevo al menú
        elif flag==5:
            flag = 1    #Flag=1 para salir del blucle
            print("Gracias por su atención.")
            print("Hasta la próxima.")
main()
