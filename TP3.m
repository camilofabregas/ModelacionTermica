clear all
close all

#----------------------Inicializacion de variables------------------------------

To = 40; #Temperatura inicial del cuerpo.
Tinf = 20; #Temperatura del fluido (fija).
h = 0.1; #Paso.
tau = 1; #Constante.


#--------------------------Funciones auxiliares---------------------------------

function graf = graficar(x, y, titulo, nom_x, nom_y)
  figure();
  plot(x,y);
  title(titulo);
  xlabel(nom_x);
  ylabel(nom_y);
endfunction


function sol = solExacta(Ti, Te, t, q)
  sol = Te + (Ti - Te) * exp(-(t/q));
endfunction


#-------------------------------Ejercicio A-------------------------------------

function [Ti, t] = Euler_A(Ti, Te, h, q)
  t = h; #Tiempo, inicializado en h.
  i = 1; #Iterador.
  T(i) = t; #Vector donde se almacenan los tiempos t de c/ iteracion.
  Tt(i) = Ti; #Vector donde se almacenan las temps. Ti de c/ iteracion.
  while (abs(Te-Ti) > 0.1)
    Ti = Ti + h*((Te - Ti) / q); #Calculo la nueva Ti.
    t += h;
    i += 1;
    Tt(i) = Ti;
    T(i) = t;
  endwhile
  #Resultado.
  disp(['[EULER] La diferencia para el paso ',num2str(h),'seg es menor a 0.1 en t=']);
  disp(t);
  disp('Euler, Orden 1=');
  disp(Ti);
  #Grafico.
  #graficar(T, Tt, 'Comportamiento de la Temp (Euler)', 'Tiempo', 'Temperatura');
endfunction


function [Ti, t] = CrankNicolson_A(Ti, Te, h, q)
  t = h; #Tiempo, inicializado en h.
  i = 1; #Iterador.
  T(i) = t; #Vector donde se almacenan los tiempos t de c/ iteracion.
  Tt(i) = Ti; #Vector donde se almacenan las temps. Ti de c/ iteracion.
  while (abs(Te-Ti) > 0.1)
    Ti_Euler = Ti + h*((Te - Ti) / q); #Calculo de Ti en adelanto con Euler.
    Ti = Ti + (h / 2)*(((Te - Ti) / q) + ((Te - Ti_Euler) / q));
    t += h;
    i += 1;
    Tt(i) = Ti;
    T(i) = t;
  endwhile
  disp(['[C-N] La diferencia para el paso ',num2str(h),'seg es menor a 0.1 en t=']);
  disp(t);
  disp('C-N, Orden 2=');
  disp(Ti);
  #Grafico.
  #graficar(T, Tt, 'Comportamiento de la Temp (C-N)', 'Tiempo', 'Temperatura');
endfunction


function [Ti, t] = RungeKutta4_A(Ti, Te, h, q)
  t = h; #Tiempo, inicializado en h.
  i = 1; #Iterador.
  T(i) = t; #Vector donde se almacenan los tiempos t de c/ iteracion.
  Tt(i) = Ti; #Vector donde se almacenan las temps. Ti de c/ iteracion.
  while (abs(Te-Ti) > 0.1)
    k1 = h*((Te - Ti) / q);
    k2 = h*((Te - (Ti + k1/2)) / q);
    k3 = h*((Te - (Ti + k2/2)) / q);
    k4 = h*((Te - (Ti + k3)) / q);
    Ti = Ti + 1/6 *(k1 + 2*k2 + 2*k3 + k4); #Calculo la nueva Ti
    t += h;
    i += 1;
    Tt(i) = Ti;
    T(i) = t;
  endwhile
  #Resultado.
  disp(['[RK-4] La diferencia para el paso ',num2str(h),'seg es menor a 0.1 en t=']);
  disp(t);
  disp('Runge-Kutta, Orden 4=');
  disp(Ti);
  #Grafico.
  #graficar(T, Tt, 'Comportamiento de la Temp (RK4)', 'Tiempo', 'Temperatura');
endfunction


#[A_Ti_Euler, A_t_Euler] = Euler_A(To, Tinf, h, tau);
#[A_Ti_CN, A_t_CN] = CrankNicolson_A(To, Tinf, h, tau);
#[A_Ti_RK4, A_t_RK4] = RungeKutta4_A(To, Tinf, h, tau);


#-------------------------------Ejercicio B-------------------------------------

function [Ti,t] = Euler_B(Ti, h, q)
  Te = 10; #Temp inicial del fluido en t=0.
  t = h; #Tiempo, inicializado en h.
  i = 1; #Iterador.
  T(i) = t; #Vector donde se almacenan los tiempos t de c/ iteracion.
  Tt(i) = Ti; #Vector donde se almacenan las temps. Ti de c/ iteracion.
  while (abs(Te-Ti) > 0.1)
    Te = 10 + sin(10*t);
    Ti = Ti + h*((Te - Ti) / q);
    t += h;
    i += 1;
    Tt(i) = Ti;
    T(i) = t;
  endwhile
  #Resultado.
  disp(['[EULER] La diferencia para el paso ',num2str(h),'seg es menor a 0.1 en t=']);
  disp(t);
  disp('Euler, Orden 1=');
  disp(Ti);
  #Grafico.
  graficar(T, Tt, 'Comportamiento de la Temp (Euler)', 'Tiempo', 'Temperatura');
endfunction


function [Ti,t] = CrankNicolson_B(Ti, h, q)
  Te = 10; #Temp inicial del fluido en t=0.
  t = h; #Tiempo, inicializado en h.
  i = 1; #Iterador.
  T(i) = t; #Vector donde se almacenan los tiempos t de c/ iteracion.
  Tt(i) = Ti; #Vector donde se almacenan las temps. Ti de c/ iteracion.
  while (abs(Te-Ti) > 0.1)
    Te = 10 + sin(10*t);
    Ti_Euler = Ti + h*((Te - Ti) / q); #Calculo de Ti en adelanto con Euler.
    Ti = Ti + (h / 2)*(((Te - Ti) / q) + ((Te - Ti_Euler) / q));
    t += h;
    i += 1;
    Tt(i) = Ti;
    T(i) = t;
  endwhile
  disp(['[C-N] La diferencia para el paso ',num2str(h),'seg es menor a 0.1 en t=']);
  disp(t);
  disp('C-N, Orden 2=');
  disp(Ti);
  #Grafico.
  graficar(T, Tt, 'Comportamiento de la Temp (C-N)', 'Tiempo', 'Temperatura');
endfunction


function [Ti,t] = RungeKutta4_B(Ti, h, q)
  Te = 10; #Temp inicial del fluido en t=0.
  t = h; #Tiempo, inicializado en h.
  i = 1; #Iterador.
  T(i) = t; #Vector donde se almacenan los tiempos t de c/ iteracion.
  Tt(i) = Ti; #Vector donde se almacenan las temps. Ti de c/ iteracion.
  while (abs(Te-Ti) > 0.1)
    Te = 10 + sin(10*t);
    k1 = h*(Te/q- Ti/q);
    k2 = h*((Te - (Ti + k1/2)) / q);
    k3 = h*((Te - (Ti + k2/2)) / q);
    k4 = h*((Te - (Ti + k3)) / q);
    Ti = Ti + 1/6 *(k1 + 2*k2 + 2*k3 + k4);
    t += h;
    i += 1;
    Tt(i) = Ti;
    T(i) = t;
  endwhile
  #Resultado.
  disp(['[RK4] La diferencia para el paso ',num2str(h),'seg es menor a 0.1 en t=']);
  disp(t);
  disp('Runge-Kutta, Orden 4=');
  disp(Ti);
  #Grafico.
  graficar(T, Tt, 'Comportamiento de la Temp (RK4)', 'Tiempo', 'Temperatura');
endfunction


#[B_Ti_Euler, B_t_Euler] = Euler_B(To, h, tau);
#[B_Ti_CN, B_t_CN] = CrankNicolson_B(To, h, tau);
#[B_Ti_RK4, B_t_RK4] = RungeKutta4_B(To, h, tau);


#-------------------------------Ejercicio C-------------------------------------

#-----EULER-----EULER-----EULER-----EULER-----EULER-----EULER-----EULER-----#

[C_Ti_Euler_005, C_t_Euler_005] = Euler_A(To, Tinf, 0.05, tau);
C_Ti_Euler_005_Exacta = solExacta(To, Tinf, C_t_Euler_005, tau)
C_Ti_Euler_005_ErrTrunc = abs(C_Ti_Euler_005_Exacta - C_Ti_Euler_005)

[C_Ti_Euler_001, C_t_Euler_001] = Euler_A(To, Tinf, 0.01, tau);
C_Ti_Euler_001_Exacta = solExacta(To, Tinf, C_t_Euler_001, tau)
C_Ti_Euler_001_ErrTrunc = abs(C_Ti_Euler_001_Exacta - C_Ti_Euler_001)

[C_Ti_Euler_0005, C_t_Euler_0005] = Euler_A(To, Tinf, 0.005, tau);
C_Ti_Euler_0005_Exacta = solExacta(To, Tinf, C_t_Euler_0005, tau)
C_Ti_Euler_0005_ErrTrunc = abs(C_Ti_Euler_0005_Exacta - C_Ti_Euler_0005)

[C_Ti_Euler_0001, C_t_Euler_0001] = Euler_A(To, Tinf, 0.001, tau);
C_Ti_Euler_0001_Exacta = solExacta(To, Tinf, C_t_Euler_0001, tau)
C_Ti_Euler_0001_ErrTrunc = abs(C_Ti_Euler_0001_Exacta - C_Ti_Euler_0001)


#-----CN-----CN-----CN-----CN-----CN-----CN-----CN-----#

[C_Ti_CN_005, C_t_CN_005] = CrankNicolson_A(To, Tinf, 0.05, tau);
C_Ti_CN_005_Exacta = solExacta(To, Tinf, C_t_CN_005, tau)
C_Ti_CN_005_ErrTrunc = abs(C_Ti_CN_005_Exacta - C_Ti_CN_005)

[C_Ti_CN_001, C_t_CN_001] = CrankNicolson_A(To, Tinf, 0.01, tau);
C_Ti_CN_001_Exacta = solExacta(To, Tinf, C_t_CN_001, tau)
C_Ti_CN_001_ErrTrunc = abs(C_Ti_CN_001_Exacta - C_Ti_CN_001)

[C_Ti_CN_0005, C_t_CN_0005] = CrankNicolson_A(To, Tinf, 0.005, tau);
C_Ti_CN_0005_Exacta = solExacta(To, Tinf, C_t_CN_0005, tau)
C_Ti_CN_0005_ErrTrunc = abs(C_Ti_CN_0005_Exacta - C_Ti_CN_0005)

[C_Ti_CN_0001, C_t_CN_0001] = CrankNicolson_A(To, Tinf, 0.001, tau);
C_Ti_CN_0001_Exacta = solExacta(To, Tinf, C_t_CN_0001, tau)
C_Ti_CN_0001_ErrTrunc = abs(C_Ti_CN_0001_Exacta - C_Ti_CN_0001)


#-----RK4-----RK4-----RK4-----RK4-----RK4-----RK4-----RK4-----#

[C_Ti_RK4_005, C_t_RK4_005] = RungeKutta4_A(To, Tinf, 0.05, tau);
C_Ti_RK4_005_Exacta = solExacta(To, Tinf, C_t_RK4_005, tau)
C_Ti_RK4_005_ErrTrunc = abs(C_Ti_RK4_005_Exacta - C_Ti_RK4_005)

[C_Ti_RK4_001, C_t_RK4_001] = RungeKutta4_A(To, Tinf, 0.01, tau);
C_Ti_RK4_001_Exacta = solExacta(To, Tinf, C_t_RK4_001, tau)
C_Ti_RK4_001_ErrTrunc = abs(C_Ti_RK4_001_Exacta - C_Ti_RK4_001)

[C_Ti_RK4_0005, C_t_RK4_0005] = RungeKutta4_A(To, Tinf, 0.005, tau);
C_Ti_RK4_0005_Exacta = solExacta(To, Tinf, C_t_RK4_0005, tau)
C_Ti_RK4_0005_ErrTrunc = abs(C_Ti_RK4_0005_Exacta - C_Ti_RK4_0005)

[C_Ti_RK4_0001, C_t_RK4_0001] = RungeKutta4_A(To, Tinf, 0.001, tau);
C_Ti_RK4_0001_Exacta = solExacta(To, Tinf, C_t_RK4_0001, tau)
C_Ti_RK4_0001_ErrTrunc = abs(C_Ti_RK4_0001_Exacta - C_Ti_RK4_0001)


#X = categorical({'Euler','Crank-Nicolson','Runge-Kutta 4'});
#y = [C_Ti_Euler_005_ErrTrunc C_Ti_Euler_001_ErrTrunc C_Ti_Euler_0005_ErrTrunc C_Ti_Euler_0001_ErrTrunc; C_Ti_CN_005_ErrTrunc C_Ti_CN_001_ErrTrunc C_Ti_CN_0005_ErrTrunc C_Ti_CN_0001_ErrTrunc; C_Ti_RK4_005_ErrTrunc C_Ti_RK4_001_ErrTrunc C_Ti_RK4_0005_ErrTrunc C_Ti_RK4_0001_ErrTrunc]
#bar(y)
#set(gca,'xticklabel',{'Euler','Crank-Nicolson','Runge-Kutta 4'});
#title('Error de truncamiento. Pasos 0.05, 0.01, 0.005, 0.001');

#-------------------------------Ejercicio D-------------------------------------

#[D_Ti_Euler_1, D_t_Euler_1] = Euler_A(To, Tinf, 1, tau);
#[D_Ti_CN_1, D_t_CN_1] = CrankNicolson_A(To, Tinf, 1, tau);
#[D_Ti_RK4_1, D_t_RK4_1] = RungeKutta4_A(To, Tinf, 1, tau);

#[D_Ti_Euler_15, D_t_Euler_15] = Euler_A(To, Tinf, 1.5, tau);
#[D_Ti_CN_15, D_t_CN_15] = CrankNicolson_A(To, Tinf, 1.5, tau);
#[D_Ti_RK4_15, D_t_RK4_15] = RungeKutta4_A(To, Tinf, 1.5, tau);

#[D_Ti_Euler_19, D_t_Euler_19] = Euler_A(To, Tinf, 1.9, tau);
#[D_Ti_CN_19, D_t_CN_19] = CrankNicolson_A(To, Tinf, 1.9, tau);
#[D_Ti_RK4_19, D_t_RK4_19] = RungeKutta4_A(To, Tinf, 1.9, tau);

#[D_Ti_Euler_2, D_t_Euler_2] = Euler_A(To, Tinf, 2, tau);        #ROMPE EN 2seg
#[D_Ti_CN_2, D_t_CN_2] = CrankNicolson_A(To, Tinf, 2, tau);      #ROMPE EN 2seg
#D_Ti_RK4_2, D_t_RK4_2] = RungeKutta4_A(To, Tinf, 2, tau);



#[D_Ti_Euler_2, D_t_Euler_2] = Euler_A(To, Tinf, 2.8, tau);        #ROMPE EN 2seg
#[D_Ti_CN_2, D_t_CN_2] = CrankNicolson_A(To, Tinf, 2.8, tau);      #ROMPE EN 2seg
#[D_Ti_RK4_28, D_t_RK4_28] = RungeKutta4_A(To, Tinf, 2.8, tau); #ROMPE EN 2.8seg

#-------------------------------Ejercicio E-------------------------------------

#---------------Cambio de Tau para RK4---------------#
#[E_Ti_RK4_tau5, E_t_RK4_tau5] = RungeKutta4_A(To, Tinf, h, 5);
#[E_Ti_RK4_tau100, E_t_RK4_tau100] = RungeKutta4_A(To, Tinf, h, 100);
#[E_Ti_RK4_tau500, E_t_RK4_tau500] = RungeKutta4_A(To, Tinf, h, 500);

#---------------Cambio de To para RK4---------------#
#[E_Ti_RK4_To150, E_t_RK4_To150] = RungeKutta4_A(150, Tinf, h, tau);
#[E_Ti_RK4_To500, E_t_RK4_To500] = RungeKutta4_A(500, Tinf, h, tau);
#[E_Ti_RK4_ToNeg, E_t_RK4_ToNeg] = RungeKutta4_A(-20, Tinf, h, tau);

#---------------Cambio de Tinf para RK4---------------#
#[E_Ti_RK4_Tinf150, E_t_RK4_Tinf150] = RungeKutta4_A(To, 150, h, tau);
#[E_Ti_RK4_Tinf500, E_t_RK4_Tinf500] = RungeKutta4_A(To, 500, h, tau);
#[E_Ti_RK4_TinfNeg, E_t_RK4_TinfNeg] = RungeKutta4_A(To, -20, h, tau);


#-------------------------------Ejercicio F-------------------------------------

function derivExactaCentrada = derivadaExactaCentrada(To, Tinf, h, t, tau)
  fSuperior = solExacta(To, Tinf, t+h, tau);
  fInferior = solExacta(To, Tinf, t-h, tau);
  derivExactaCentrada = (fSuperior - fInferior) / (2*h);
endfunction

#derivExactaCentradaT1 = derivadaExactaCentrada(To, Tinf, h, 1, tau)
#derivExactaCentradaT2 = derivadaExactaCentrada(To, Tinf, h, 2, tau)
#derivExactaCentradaT3 = derivadaExactaCentrada(To, Tinf, h, 3, tau)


function derivExactaAdelanto = derivadaExactaAdelanto(To, Tinf, h, t, tau)
  f3h = solExacta(To, Tinf, t+3*h, tau);
  f2h = solExacta(To, Tinf, t+2*h, tau);
  f1h = solExacta(To, Tinf, t+h, tau);
  ft = solExacta(To, Tinf, t, tau);
  derivExactaAdelanto = ((f3h / 3) - (f2h * (3/2)) + (f1h * 3) - (ft * (11/6))) / h;
endfunction

#derivExactaAdelantoT1 = derivadaExactaAdelanto(To, Tinf, h, 1, tau)
#derivExactaAdelantoT2 = derivadaExactaAdelanto(To, Tinf, h, 2, tau)
#derivExactaAdelantoT3 = derivadaExactaAdelanto(To, Tinf, h, 3, tau)


function derivExactaAtraso = derivadaExactaAtraso(To, Tinf, h, t, tau)
  f3h = solExacta(To, Tinf, t-3*h, tau);
  f2h = solExacta(To, Tinf, t-2*h, tau);
  f1h = solExacta(To, Tinf, t-h, tau);
  ft = solExacta(To, Tinf, t, tau);
  derivExactaAtraso = ((ft * (11/6)) - (f1h * 3) + (f2h * (3/2)) - (f3h / 3)) / h;
endfunction

#derivExactaAtrasoT1 = derivadaExactaAtraso(To, Tinf, h, 1, tau)
#derivExactaAtrasoT2 = derivadaExactaAtraso(To, Tinf, h, 2, tau)
#derivExactaAtrasoT3 = derivadaExactaAtraso(To, Tinf, h, 3, tau)