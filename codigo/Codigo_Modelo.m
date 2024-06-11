# ANÁLISIS DINÁMICO DE LA SUSPENSIÓN DE UN VEHÍCULO CON UN RECIPIENTE CON LÍQUIDO EN EL INTERIOR

clc; clear; close;

###################################
### FUNCION DE PICOS ALEATORIOS ###
###################################
% Funcion de amplitud erratica para camino aleatorio
function amplitud = amplitud_aleatoria(t)
    amplitud_base = 0.02*exp(-0.1*t).*sin(2*pi*t/5);
    impulsos = abs(rand(size(t)))>0.9; % Genera impulsos aleatorios con probabilidad del 90%
    amplitud = amplitud_base.*(1+5*randn(size(t))).*(1+impulsos*80000); % Combina la amplitud base, el ruido y los impulsos
end

##################
### PARAMETROS ###
##################

# Automovil
H = 1.46; % Altura del automovil [m]
L = 4.49; % Longitud del chasis [m]
r = 0.5; % Distancia al centro de masa geometrico [m]
l_eje = 2.74; % Distancia entre ejes [m]
l2 = 1.1; % Distancia del centro de masa a rueda delantera [m]
l1 = l_eje-l2; % Distancia del centro de masa a rueda trasera [m]

# Masa
m2 = 25; % Masa rueda frontal [kg]
m1 = 25; % Masa rueda trasera [kg]
m3 = 985; % Masa del chasis [kg] 
I = 1/12*m3*(L^2+H^2) + m3*r^2; % Inercia y Steiner [kg.m2]

# Rigidez
k1 = 25530; % Suspension trasera [N/m]
k2 = 35350; % Suspension frontal [N/m]
k3 = 154850; % Rigidez rueda trasera [N/m]
k4 = 154850; % Rigidez rueda frontal [N/m]

% Amortiguamiento
c1 = 1100; % Amortiguador trasero [N.s/m]
c2 = 1200; % Amortiguador frontal [N.s/m]

# Desplazamiento del Vehiculo
long_onda = 5; % [m]
amplitud_onda = 0.05; % [m]
velocidad = 20; % [m/s]
bache = 0.2; % Tamaño del bache [m]
largo_bache = 0.5; % Largo/diametro del bache [m]
duracion = largo_bache/velocidad; % Duración del impulso [s]

# Vector de tiempo
max_t = 10; % Cuantos segundos quiero observar [s]
vect_t = linspace(0, max_t, 3000); % Vector de tiempo

# Fuerza externa aplicada como impulso
t_pulso_1 = 5; % Tiempo en el que se aplica el impulso a la rueda delantera [s]
t_pulso_2 = t_pulso_1 + (l2+l1)/velocidad; % Tiempo en el que se aplica el impulso a la rueda trasera [s]

###################################
### RECIPIENTE DE CARGA LIQUIDA ###
###################################
d = 0.06; % Diametro del recipiente [m]
hL0 = 0.08; % Altura del liquido en reposo [m]
hr = 0.1; % Altura del recipiente [m]
h_taza = 1.2; % Posicion/altura del recipiente en el auto [m]
g = 9.81; % Aceleracion de la gravedad [m/s2]

#############################
### CONDICIONES INICIALES ###
#############################
% Condiciones iniciales
x0D = [0; 0; 0; 0]; % Desplazamientos
x0V = [0; 0; 0; 0]; % Velocidades

################
### MATRICES ###
################
# Matriz de Masa
mat_M = diag([m2, m1, m3, I]);

# Matriz de Rigidez
mat_K = [(k4+k2),      0,           -k2,           -k2*l2;
       0,       (k3+k1),        -k1,            k1*l1; 
      -k2,        -k1,        (k2+k1),      (k2*l2-k1*l1); 
    -k2*l2,      k1*l1,    (k2*l2-k1*l1), (k2*l2^2+k1*l1^2)];

# Matriz de Amortiguamiento
mat_C = [c2,          0,          -c2,           -c2*l2;
     0,           c1,         -c1,            c1*l1; 
    -c2,         -c1,        (c2+c1),      (c2*l2-c1*l1); 
   -c2*l2,      c1*l1,    (c2*l2-c1*l1), (c2*l2^2+c1*l1^2)]; 


########################
### FUERZAS EXTERNAS ###
########################

# Fuerza Externa aplicada en los Neumaticos
F1 = @(t) [bache*k4*(t>=t_pulso_1 & t<t_pulso_1+duracion); bache*k3*(t>=t_pulso_2 & t<t_pulso_2+duracion); 0; 0];

# Fuerza externa aplicada por el camino
Ff = amplitud_onda*k4; % Fuerza sobre rueda delantera
Fr = amplitud_onda*k3; % Fuerza sobre rueda trasera

w_f = 2*pi*velocidad/long_onda; % Frecuencia forzada
phi = (l_eje)*2*pi/long_onda; % Desfasaje de la rueda trasera

# Camino con relieve senoidal
F2 = @(t) [Ff*sin(w_f*t); Fr*sin(w_f*t+phi); 0; 0];

# Fuerza externa aleatoria simulando un camino sinuoso
F3 = @(t) [amplitud_aleatoria(t); amplitud_aleatoria(t); 0; 0];

# Fuerza combinada
F = @(t) F1(t)+F2(t)+F3(t);

##################
### AUXILIARES ###
##################
% Calculo de autovalores y autovectores
[Autovectores, Autovalores] = eig(mat_K, mat_M); % V: autovectores (modos), D: autovalores

% Calculo de frecuencias naturales
vect_w_n = sqrt(diag(Autovalores)); # Frecuencia Natual del Sistema
vect_zita = diag(Autovectores'*mat_C*Autovectores)./(2*vect_w_n); % Vector de Coeficientes de amortiguacion modal


############################
### DESCOMPOSICION MODAL ###
############################
% Transformacion a coordenadas modales de las condiciones iniciales
q0D = Autovectores'*mat_M*x0D; # Condiciones Iniciales de desplazamiento modal
q0V = Autovectores'*mat_M*x0V; # Condiciones Iniciales de velocidad modal

% Inicializacion de matriz modal
q_modal = zeros(length(vect_w_n), length(vect_t));

% Precomputo de fuerza externa a lo largo del tiempo (para resolver luego la integral)
F_t = cell2mat(arrayfun(F, vect_t, 'UniformOutput', false));


#######################################
### SOLUCION EN COORDENADAS MODALES ###
#######################################
for i = 1:length(vect_w_n)
    w_d = vect_w_n(i)*sqrt(1-vect_zita(i)^2); % Frecuencia amortiguada

    # Auxiliar 1
    A = q0D(i);

    # Auxiliar 2
    B = (q0V(i)+vect_zita(i)*vect_w_n(i)*q0D(i))/w_d;
    
    # Inicializacion de vector para solucion permanente
    solucion_permanente = zeros(1, length(vect_t));
    
    # Resolucion de integral para solucion permanente
    for j = 1:length(vect_t)
        tau = vect_t(1:j);
        Q_modal = Autovectores(:,i)'*F_t(:, 1:j); % Transformacion de la fuerza en el tiempo a espacio modal
        
        integrand = exp(-vect_zita(i)*vect_w_n(i)*(vect_t(j)-tau)).*sin(w_d*(vect_t(j)-tau)).*Q_modal;
        solucion_permanente (j) = trapz(tau, integrand, 2); % Metodo del trapecio para resolver
    end
    
    % Resolucion de la parte transitoria + permanente
    q_modal(i, :) = exp(-vect_zita(i)*vect_w_n(i)*vect_t) .* (A*cos(w_d*vect_t)+B*sin(w_d*vect_t)) + (1/(w_d*(1-vect_zita(i)^2)))*solucion_permanente;
end


###############################
### COORDENADAS GEOMETRICAS ###
###############################
% Transformacion a coordenadas geometricas
x_total = Autovectores*q_modal;
vect_ang_tita = x_total(4,:)*180/pi; % Desplazamiento angular del chasis en grados


###############################
### VELOCIDAD Y ACELERACION ###
###############################

# Velocidades
aprox_Velocidad = diff(x_total, 1, 2) ./ diff(vect_t); % Aproximacion de la velocidad
t_velocidad = vect_t(1:end-1); % El vector de tiempo para velocidades es un punto mas corto

% Aceleraciones
aprox_Aceleracion = diff(aprox_Velocidad, 1, 2) ./ diff(t_velocidad); % Aproximacion de la aceleracion
t_aceleracion = t_velocidad(1:end-1); % El vector de tiempo para aceleraciones es dos puntos mas corto

# Aceleracion Traslacional
a_lineal = aprox_Aceleracion(4,:)*h_taza; % Convertir de Movimiento rotacional a traslacional del chasis


##########################
#### ESTADO DEL FLUIDO ###
##########################
# Altura del Fluido
vect_h_fluido = hL0 + abs(a_lineal)*d/2./(g+abs(aprox_Aceleracion(3,:))); 

disp("Estado del Sistema")
if max(vect_h_fluido)>hr % Comparacion con la altura
    disp('-- El liquido se derramo. --');
else
    disp('-- El liquido NO se derramo. --');
end

vect_h_recipiente = ones(1,length(t_aceleracion)) * hr; % Vector para graficar altura limite del liquido



###################
#### RESPUESTAS ###
###################
disp("\nModos de Vibracion del Sistema")
disp(Autovectores)

disp("\nFrecuencias Natural por cada GLD")
disp(vect_w_n)

disp("\nMaxima altura  alcanzada por el liquido[m]")
disp(max(vect_h_fluido))

#################
#### GRAFICOS ###
#################

% Figura 1: Fuerzas Externas
figure(1)
subplot(2,1,1)
plot(vect_t, F_t(1,:), 'r', LineWidth=1)
title('Fuerza sobre rueda delantera')
xlabel('Tiempo [s]')
ylabel('Fuerza [N]')
grid on

subplot(2,1,2)
plot(vect_t, F_t(2,:), 'b', LineWidth=1)
title('Fuerza sobre rueda trasera')
xlabel('Tiempo (s)')
ylabel('Fuerza (N)')
grid on

% Figura 2: Comportamiento del sistema
figure(2)
subplot(3,1,1)
plot(vect_t, x_total(1,:), 'b', vect_t, x_total(2,:), 'm', LineWidth=1)
title('Desplazamientos');
legend('Rueda frontal', 'Rueda trasera');
xlabel('Tiempo [s]')
ylabel('Desplazamiento [m]')
grid on

subplot(3,1,2)
plot( vect_t, x_total(3,:), 'b', LineWidth=1)
title('Desplazamiento del chasis')
xlabel('Tiempo [s]')
ylabel('Desplazamiento [m]')
grid on

subplot(3,1,3)
plot(vect_t, vect_ang_tita, 'k', LineWidth=1)
title('Desplazamiento angular en grados (cabeceo)')
xlabel('Tiempo [s]')
ylabel('Angulo [°]')
grid on

% Figura 3:Comportamiento del Fluido
figure(3)
plot(t_aceleracion, vect_h_fluido, 'b', t_aceleracion, vect_h_recipiente, 'r', LineWidth=1)
legend('Maxima Altura del Fluido', 'Borde del Recipiente');
title('Movimiento del líquido')
axis([0, max_t, 0.07, 0.11]);
xlabel('Tiempo [s]')
ylabel('Desplazamiento [m]')
grid on