# -*- coding: utf-8 -*-
#########################################################################
  #  Aplicação do modelo de Hodgkin-Huxley pelo método de Runge_Kutta
  #
  #         @Autor: Anderson Ferreira Sepulveda
  #         Universidade de São Paulo - Instituto de Física 
#########################################################################

import math
import matplotlib.pyplot as plt

def Volt(V_i,C_m,i_m,I_e,A,delta):
    '''
    Potencial de membrana
    '''    
    k1 = delta*((1./C_m)*(-i_m + I_e/A))
    k2 = delta*((1./C_m)*(-i_m + I_e/A + 0.5*k1))
    k3 = delta*((1./C_m)*(-i_m + I_e/A + 0.5*k2))
    k4 = delta*((1./C_m)*(-i_m + I_e/A + k3))
    V = V_i + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)    
    return V

def Alfa_n(V):
    '''
    Os valores de alfa_n
    '''     
    a_n = (0.01*(V + 55.))/(1. - math.e**(-(V + 55.)/10.))
    return a_n

def Beta_n(V):
    '''
    Os valores de beta_n
    '''    
    b_n = 0.125*math.exp(-(V + 65)/80.)
    return b_n

def Alfa_m(V):
    '''
    Os valores de alfa_m
    '''
    a_m =  (0.1*(V + 40.))/(1. - math.e**(-(V + 40.)/10.))
    return a_m
    
def Beta_m(V):
    '''
    Os valores de beta_m
    '''    
    b_m = 4.*math.e**(-(V + 65.)/18.)
    return b_m

def Alfa_h(V):
    '''
    Os valores de alfa_h
    ''' 
    a_h = 0.07*math.e**(-(V + 65.)/20.)
    return a_h

def Beta_h(V):
    '''
    Os valores de beta_h
    '''    
    b_h = 1./(1. + math.e**(-(V + 35.)/10.))
    return b_h
    
    
def N(delta,V,n_i):
    '''
    Calcula a probabilidade de n
    '''
    k1 = delta*(Alfa_n(V)*(1. - n_i) - Beta_n(V)*n_i)
    k2 = delta*(Alfa_n(V)*(1. - (n_i + 0.5*k1)) - Beta_n(V)*(n_i + 0.5*k1))
    k3 = delta*(Alfa_n(V)*(1. - (n_i + 0.5*k2)) - Beta_n(V)*(n_i + 0.5*k2))
    k4 = delta*(Alfa_n(V)*(1. - (n_i + k3)) - Beta_n(V)*(n_i + k3))
    n = n_i + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
    return n
    
def M(delta,V,m_i):
    '''
    Calcula a probabilidade de m
    '''
    k1 = delta*(Alfa_m(V)*(1. - m_i) - Beta_m(V)*m_i)
    k2 = delta*(Alfa_m(V)*(1. - (m_i + 0.5*k1)) - Beta_m(V)*(m_i + 0.5*k1))
    k3 = delta*(Alfa_m(V)*(1. - (m_i + 0.5*k2)) - Beta_m(V)*(m_i + 0.5*k2))
    k4 = delta*(Alfa_m(V)*(1. - (m_i + k3)) - Beta_m(V)*(m_i + k3))
    m = m_i + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
    return m
    
def H(delta,V,h_i):
    '''
    Calcula a probabilidade h
    '''
    k1 = delta*(Alfa_h(V)*(1. - h_i) - Beta_h(V)*h_i)
    k2 = delta*(Alfa_h(V)*(1. - (h_i + 0.5*k1)) - Beta_h(V)*(h_i + 0.5*k1))
    k3 = delta*(Alfa_h(V)*(1. - (h_i + 0.5*k2)) - Beta_h(V)*(h_i + 0.5*k2))
    k4 = delta*(Alfa_h(V)*(1. - (h_i + k3)) - Beta_h(V)*(h_i + k3))
    h = h_i + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
    return h

                
def main():    
    
    '''
    Vetores
    '''    
    I = []
    
    
    V_11 = []  
    
    V_e = []
        
    
    T_11 =[]    
        
    
    MUm_11 = []   
           
    
    MUh_11 =[]   
           
    
    MUn_11 =[]   
        
    N_e = []
    M_e = []
    H_e = []
    
    N_inf = []
    M_inf = []
    H_inf = []
    
    tau_N = []
    tau_M = []
    tau_H = []
    
    '''
    Parâmetros
    '''
    F = 96.5                    #Constante de Faraday (coulombs/mmol)
    R = 8.314                   #COnstante dos gases (J/K)
    T_celsius = 6.3             #Temperatura a qual é submetido o axônio (em celsius) 
    T = T_celsius + 273         #Temperatura em K
    
    t_i = 0.      #Instante inicial (ms)
    t_f = 10     #Instante final  (ms)
    V_rep = -65.   #Potencial de repouso (em mV)    
    #R_m = 10000.    #Resistência da membrana
    C_m = 0.1        #Capacitância da membrana (uF/cm^2)
    A = 0.001      #Área da superfície celular (cm^2)
    I_e = 0       #Corrente injetada na célula    
    
    G_L = 0.3     #Condutância (mS/cm^2)
    G_K = 36.      #Condutância a ions potassio
    G_Na = 120.    #Condutancia a ions sodio
    
    #tau = 10.      #Constante de tempo (ms)
    
    Na_o = 491.    #Concentração extracelular de sodio (mmol/l)
    Na_i = 50.     #Concentração intracelular de sódio
    K_o = 20.      #Concentração extracelular de potassio
    K_i = 400.     #Concentração intracelular de potassio
    
    E_L = -49.    #Força eletromotriz
    E_Na = R*T*math.log(Na_o/Na_i)/F
    E_K = R*T*math.log(K_o/K_i)/F
    
    n_i = Alfa_n(V_rep)/(Alfa_n(V_rep) + Beta_n(V_rep))
    m_i = Alfa_m(V_rep)/(Alfa_n(V_rep) + Beta_m(V_rep))
    h_i = Alfa_h(V_rep)/(Alfa_h(V_rep) + Beta_h(V_rep))
    
    '''
    Inicializadores
    '''
    
    V_11.append(V_rep)    
    
    V_e.append(V_rep)
    
      
    MUn_11.append(n_i)
    
    
    MUm_11.append(m_i)
    
    
    MUh_11.append(h_i)
    
    N_inf.append(Alfa_n(V_rep)/(Alfa_n(V_rep) + Beta_n(V_rep)))
    M_inf.append(Alfa_m(V_rep)/(Alfa_m(V_rep) + Beta_m(V_rep)))
    H_inf.append(Alfa_h(V_rep)/(Alfa_h(V_rep) + Beta_h(V_rep)))
    
    tau_N.append(1./(Alfa_n(V_rep) + Beta_n(V_rep)))
    tau_M.append(1./(Alfa_m(V_rep) + Beta_n(V_rep)))
    tau_H.append(1./(Alfa_h(V_rep) + Beta_h(V_rep)))
       
    N_e.append(n_i)
    M_e.append(m_i)
    H_e.append(h_i)
       
    
    T_11.append(t_i)
        
    '''
    Aplicação do modelo de H-H quando mu = 11 (que mais aproxima do valor exato
    '''    
    mu = 11    
    I_Na = G_Na*(MUm_11[0]**3)*MUh_11[0]*(V_11[0] - E_Na)
    I_K = G_K*(MUn_11[0]**4)*(V_11[0] - E_K)
    I_L = G_L*(V_11[0] - E_L)
    I_m = I_Na + I_K + I_L
    I.append(I_m)
    for j in range(4**mu):
        DeltaT = (t_f - t_i)/(4.**mu)
        I_Na = G_Na*(MUm_11[j]**3)*MUh_11[j]*(V_11[j] - E_Na)
        I_K = G_K*(MUn_11[j]**4)*(V_11[j] - E_K)
        I_L = G_L*(V_11[j] - E_L)
        I_m = I_Na + I_K + I_L
        I.append(I_m)
        v = Volt(V_11[j],C_m,I_m,I_e,A,DeltaT)
        V_11.append(v)        
        n =  N(DeltaT,V_11[j],MUn_11[j])
        MUn_11.append(n)
        m = M(DeltaT,V_11[j],MUm_11[j])
        MUm_11.append(m)
        h = H(DeltaT,V_11[j],MUh_11[j])
        MUh_11.append(h)
        t = T_11[j] + DeltaT
        T_11.append(t)
    
    '''
    Valores exatos
        
    mu = 11    
    for j in range(4**mu):                    
        n_inf = Alfa_n(V_11[j])/(Alfa_n(V_11[j]) + Beta_n(V_11[j]))
        N_inf.append(n_inf)
        tau_n = 1./(Alfa_n(V_11[j]) + Beta_n(V_11[j]))
        tau_N.append(tau_n)
        
        
        m_inf = Alfa_m(V_11[j])/(Alfa_m(V_11[j]) + Beta_m(V_11[j]))
        M_inf.append(m_inf)
        tau_m = 1./(Alfa_m(V_11[j]) + Beta_m(V_11[j]))
        tau_M.append(tau_m)
        
        
        h_inf = Alfa_h(V_11[j])/(Alfa_h(V_11[j]) + Beta_h(V_11[j]))
        H_inf.append(h_inf)
        tau_h = 1./(Alfa_h(V_11[j]) + Beta_h(V_11[j]))
        tau_H.append(tau_h)
        
    '''       
       
    '''
    plot
    '''    
    '''
    plt.figure()
    plt.plot(T_11,MUn_11,'k--',label = 'n')
    plt.plot(T_11,MUm_11,'k:',label = 'm')
    plt.plot(T_11,MUh_11,'k-.',label = 'h')
    plt.title('Probabilidade de abertura dos canais com o tempo')
    plt.xlabel('Tempo (ms)')
    plt.ylabel('Probabilidade')
    plt.legend(loc = 'best')
    plt.show()
    
    
    plt.figure()
    plt.plot(V_11,N_inf,'k--',label = 'n_inf')
    plt.plot(V_11,M_inf,'k:',label = 'm_inf')
    plt.plot(V_11,H_inf,'k-.',label = 'h_inf')
    plt.xlabel('Voltagem (mV)')
    plt.legend(loc = 'best')
    plt.show()
    
    
    plt.figure()
    plt.plot(V_11,tau_N,'k--',label = '$tau_n$')
    plt.plot(V_11,tau_M,'k:',label = '$tau_m$')
    plt.plot(V_11,tau_H,'k-.',label = '$tau_h$')
    plt.xlabel('Voltagem (mV)')
    plt.ylabel('tau (ms)')
    plt.legend(loc = 'best')
    plt.show()
    
    plt.figure()
    plt.plot(T_11,I,'k--')
    plt.title('Corrente')
    plt.xlabel('Tempo (ms)')
    plt.ylabel('i ($\mu\,A/mm^2$)')
    plt.show()
    '''
    plt.figure()
    plt.subplot(5,1,1)
    plt.plot(T_11,V_11,'k:')            
    plt.ylabel('V (mV)')
    plt.subplot(5,1,2)    
    plt.plot(T_11,MUn_11,'k:')            
    plt.ylabel('n')    
    plt.subplot(5,1,3)
    plt.plot(T_11,MUm_11,'k:')                
    plt.ylabel('m')    
    plt.subplot(5,1,4)
    plt.plot(T_11,MUh_11,'k:')        
    plt.ylabel('h')    
    plt.subplot(5,1,5)
    plt.plot(T_11,I,'k:')    
    plt.xlabel('Tempo (ms)')    
    plt.ylabel('i ($\mu\,A/mm^2$)')
    plt.show()
                                      
main()
