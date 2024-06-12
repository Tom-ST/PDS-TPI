pkg load signal
clc; clear all; close all;
##variables (se pueden cambiar)
duracion = 0.1; %cantidad de minutos a leer del archivo

t_ini = 10; #Segundo que se comienza a analizar la cancion
t_fin = 70; #Termina el analisis

## Codigo
canciones = cargar_canciones();
seleccion = menu("Elige una cancion", canciones);

if seleccion > 0
  nombre_archivo = canciones{seleccion};
  disp(["Has seleccionado: ", nombre_archivo]);
else
  disp("No has seleccionado ninguna canción.");
  return;
endif

archivo_audio = strcat("../samples/",nombre_archivo);


% Cargar la cantidad de minutos especificados
info = audioinfo(archivo_audio);

#duracion_segundos = duracion * 60;
duracion_segundos = t_fin - t_ini;

##n_muestras = fix(duracion_segundos * info.SampleRate);
muestra_ini = fix(t_ini * info.SampleRate);
muestra_fin = fix(t_fin * info.SampleRate);

##if info.TotalSamples < n_muestras
##  error("El archivo es mas corto que la duracion especificada");
##endif



#[y, Fs] = audioread(archivo_audio, [1 n_muestras]);
[y, Fs] = audioread(archivo_audio, [muestra_ini muestra_fin]);
##t = (0:n_muestras-1) / Fs;
t = t_ini:1/Fs:t_fin;

%Convierto audio stereo a mono
#Audio_read nos da dos columnas, una para el sonido q sale por el parlante izquierdo
#y otra que sale por el parlante derecho, hacemos un promedio
if size(y ,2) == 2
  y = (y(:,1) + y(:,2)) / 2;
endif


##figure(1);
##plot(t,y);
##xlabel("Tiempo (s)");
##ylabel("Amplitud");
##title("Audio en funcion del tiempo")

%ventana para transformada de fourier de corto tiempo: 10ms (por que 10? se podria calcular este valor?)
%analizo la transformada de fourier en ventanas de 10 ms para detectar
% cambios en las frecuencias, solo utilizo la magnitud de la fft
% utilizo una funcion de compresion C(x) para que aquellos componentes de
% frecuencia alta (como platillos) no sean "masked(?)" por componentes
% de baja frecuencia pero alta amplitud.
% uso C(x) = x^(1/2), o arcsin(x). no uso log porque no se comporta bien cerca del cero

tam_ventana = 0.01;% en segundos (10 ms)
tam_ventana_m = fix(tam_ventana * Fs); %tamaño de la ventana en muestras
desplazamiento = fix(tam_ventana_m / 2); %desplazamiento de la ventana (50% overlap)

## Funcion de compresion
C = @(x) sqrt(x);

num_ventanas = fix(length(y) / tam_ventana_m);

fragmentos = zeros(num_ventanas, tam_ventana_m);

for i = 1:num_ventanas
  fragmentos(i,:) = y((i-1)*tam_ventana_m + 1 : i*tam_ventana_m);
endfor


fft_resultados = zeros(num_ventanas, tam_ventana_m);

for i = 1:num_ventanas
  fft_resultados(i, :) = C(abs(fft(fragmentos(i,:)))); #las filas son 10ms, hace la fft cada 10ms
endfor

##figure(2)
##stem(fft_resultados(100,:))
##title("Transformada de fourier para una cierta fila")


%----- Calculo del flujo de la energia ------
% Encontrar los índices correspondientes a 100 Hz y 10000 Hz
frecuencia_100Hz = 100;
frecuencia_10000Hz = 10000;

indice_100Hz = round(frecuencia_100Hz * tam_ventana_m / Fs);
indice_10000Hz = round(frecuencia_10000Hz * tam_ventana_m / Fs);

E_hat = zeros(num_ventanas - 1, 1); % Se resta 1 porque se calcula la diferencia entre frames sucesivos

j=(indice_100Hz + 1):indice_10000Hz;

for i = 2:num_ventanas
  E_hat(i-1) = sum(fft_resultados(i,j) - fft_resultados(i-1, j));
endfor

% Rectificacion media onda
E = max(E_hat, 0);

f_m = info.SampleRate; %frecuencia de muestreo

##figure(3);
####plot(0.01:0.01:t(end),E);
##plot((t_ini+0.01):0.01:(t_fin-0.01),E);
##
##xlabel("Tiempo");
##ylabel("Flujo de energía");
##title("Flujo de energía en función del tiempo");









% =========================== METODO DEL PAPER 2 ===============================
[peaks, peak_locs] = findpeaks(E);

% Ploteamos los picos en la señal de flujo de energía
##figure(4);
##plot((t_ini+0.01):0.01:(t_fin-0.01),E);
##hold on;
##plot((t_ini+0.01)+peak_locs*0.01, peaks, 'ro'); % Convertimos los índices a tiempo
##hold off;
##xlabel("Tiempo");
##ylabel("Flujo de energía");
##title("Picos de flujo de energía en función del tiempo");


% Establecer un umbral como un porcentaje del máximo pico encontrado
umbral = 0.5; % Por ejemplo, seleccionamos picos que estén por encima del 50% del máximo pico

% Calculamos el máximo pico
max_peak = max(peaks);

% Filtramos los picos significativos que superan el umbral
picos_significativos = peak_locs(peaks >= umbral * max_peak);
valores_significativos = peaks(peaks >= umbral * max_peak);

% Ploteamos los picos significativos en la señal de flujo de energía
figure(5);
plot((t_ini+0.01):0.01:(t_fin-0.01),E);
hold on;
plot((t_ini+0.01)+picos_significativos*0.01, valores_significativos, 'ro'); % Convertimos los índices a tiempo
hold off;

xlabel("Tiempo");
ylabel("Flujo de energía");
title("Picos de flujo de energía significativos en función del tiempo");

% Ploteamos la señal de audio en función del tiempo
##figure(6);
##plot(t, y);
##hold on;
##
##% Ploteamos líneas verticales en las ubicaciones de los picos
tiempo_picos = (t_ini + 0.01) + picos_significativos * 0.01; % Convertimos los índices de picos a tiempo
##for i = 1:length(tiempo_picos)
##    line([tiempo_picos(i), tiempo_picos(i)], ylim, 'Color', 'r', 'LineStyle', '--'); % Graficamos una línea vertical en la posición del pico
##end
##
##hold off;
##
##xlabel("Tiempo (s)");
##ylabel("Amplitud");
##title("Picos identificados en la señal de audio");



%============================ Funcion de probabilidad ==================
% La idea es probar en la funcion L_t distintas combinaciones de T, S ,y b1 hasta que se obtenga
% la mayor likelihood

% ---Avances--
% La funcion p_t ya esta terminada, te da la distribucion de probabilidad de la ubicacion de los beats
% la ubicacion de los beats 1 a 4 ya estan calculados

% Tempos, Swings y beats1 son distintos valores que se van a tener que probar para maximizar la likelihood
% por ahora conviene probar con un solo valor de esos
% no cambiar los rangos que esos estan asi por el paper

% ---Falta hacer--
% 1.
% en pt(end) el valor es 0 (hacer zoom en el grafico de las posiciones probables de los beats
% no deberia bajar asi derecho, sino tomar el valor que corresponde en la campana de gauss. Corregir

% 2.
% Hay que terminar la funcion L_t que recibe (T, S, b1, t_picos) y da un valor numerico, no se por que me da -inf
%CAMBIANDO P_T AGREGANDO UN EPSILON PARA QUE NO EVALUE EN 0 SE SOLUCIONA EL ERROR

% una vez que de un numero, hay que probar todas las combinaciones posibles de T, b1 y S y el mayor L_t es el correcto

Tempos = [70:140];
Swings = [0:0.1:0.4];%El swing (S) es un porcentaje que determina que tanto se atrasa el segundo y cuarto cuarto-beat (es una propiedad de algunos generos musicales como rock y jazz)


tiempo_seleccionado = t_fin-t_ini;
T = Tempos(1);
S = Swings(1);
periodo_beat = 60 / T;
beats1 = linspace(0, periodo_beat, 32)+t_ini;
b1 = beats1(1);
b2 = b1+periodo_beat;
b3 = b2+periodo_beat;
b4 = b3+periodo_beat;

b1_4 = [b1 b2+periodo_beat*S b3 b4+periodo_beat*S];

b1_n = b1;  % Inicializamos el primer beat
beat = b1 + periodo_beat;  % Inicializamos el siguiente beat

% Generamos los beats adicionales hasta alcanzar o superar el tiempo t_fin
while beat <= t_fin
    b1_n = [b1_n, beat];
    beat = beat + periodo_beat;  % Calculamos el siguiente beat
end

% Ajustamos los beats en las segundas y cuartas posiciones
for i = 2:2:length(b1_n)
    b1_n(i) = b1_n(i) + periodo_beat * S;
end

for i = 4:4:length(b1_n)
    b1_n(i) = b1_n(i) + periodo_beat * S;
end

##t=linspace(t_ini, t_fin,10000);
pt = p_t (t, T, b1_n);

#metromono 80beat
#Con epsilon
#T = 70
#S = 0
#b1 = 10
#retval = -1860.6

#Sin epsilon
#T = 70
#S = 0
#b1 = 10
#retval = inf

figure();
hold on
plot(t,pt)
title("Posicion probable donde se encuentran los cuarto-beats");
xlabel("Tiempo");

retval = L_t (T, S, b1, tiempo_picos)

%-------MAYOR LIKELIHOOD-----------------

% Discretización de las variables
##T_values = Tempos(1):10:Tempos(end); % Reducir el rango y aumentar los pasos de T
##S_values = linspace(0.5, 2, 10); % Reducir el rango y aumentar los pasos de S
##b0_values = linspace(0, 0.5, 10); % Reducir el rango y aumentar los pasos de b0s

T_values = 70:10:120;
S_values = 0;
b1_values=beats1;

max_likelihood = -Inf;
best_T = T_values(1);
best_S = S_values(1);
best_b1 = b1_values(1);

total_iterations = length(T_values) * length(S_values) * length(b1_values);
iteration = 0;
likelihood_values = [];
for i = 1:length(T_values)

    for j= 1:length(S_values)

      for k =1:length( b1_values)
            iteration++;
            likelihood = L_t(T_values(i), S_values(j), b1_values(k), tiempo_picos);
             likelihood_values = [likelihood_values, likelihood]; % Almacenar el valor de likelihood
            if likelihood > max_likelihood
                max_likelihood = likelihood;
                best_T = T_values(i);
                best_S = S_values(j);
                best_b1 = b1_values(k);
            endif
            if mod(iteration, 10) == 0
                disp(['Iteración ', num2str(iteration), ' de ', num2str(total_iterations)]);
            endif
        endfor
    endfor
endfor

disp(['Mejor T: ', num2str(best_T)]);
disp(['Mejor S: ', num2str(best_S)]);
disp(['Mejor b1: ', num2str(best_b1)]);
disp(['Max Likelihood: ', num2str(max_likelihood)]);



%====================== Fin




##% Calcular la función de verosimilitud para diferentes BPM
##max_bpm = 140;
##min_bpm = 70;
##bps = (min_bpm:max_bpm) / 60;
##L = zeros(length(bps), 1);
##
##for k = 1:length(bps)
##  bpm = bps(k);
##  interval = Fs / bpm; % Intervalo en muestras
##  score = 0;
##  for i = 1:length(picos_significativos)
##    % Para cada pico, calcular la distancia a los múltiplos del intervalo
##    distancias = abs(mod(picos_significativos(i) - picos_significativos(1), interval));
##    % Considerar también la distancia al múltiplo siguiente para evitar errores de fase
##    distancias = min(distancias, interval - distancias);
##    % Sumar las probabilidades
##    prob = exp(-distancias / (Fs / 4));
##    score = score + sum(prob);
##  endfor
##  L(k) = score;
##endfor
##
##[~, idx_max] = max(L);
##BPM_estimado = bps(idx_max) * 60;
##
##figure(7);
##plot((min_bpm:max_bpm), L);
##xlabel('BPM');
##ylabel('Likelihood');
##title('Likelihood de diferentes BPM');
##
##disp(["El BPM estimado es: ", num2str(BPM_estimado)]);


##
##%========================== CALCULO SIMPLE DE BPM UTILIZANDO PROMEDIO ==========================
##%========================== DESCOMENTAR PARA PROBAR ============================================
##% Calculamos los intervalos de tiempo entre los picos consecutivos
##intervalos_tiempo = diff(picos_significativos) * 0.01; % Convertimos los índices de picos a tiempo y calculamos los intervalos en segundos
##
##% Calculamos el promedio de los intervalos de tiempo
##promedio_intervalo_tiempo = mean(intervalos_tiempo);
##
##% Convertimos el intervalo de tiempo promedio a BPM
##bpm = round(60 / promedio_intervalo_tiempo);
##
##disp(['El tempo de la canción es aproximadamente ', num2str(bpm), ' BPM']);
##
##
##% Función p_t
##function pt_value = pt(t_i, T, S, b_0)
##  pt_value = exp(-((t_i - T) ^ 2) / (2 * S ^ 2)) / (S * sqrt(2 * pi));
##endfunction
##
##% Función L_t
##function log_likelihood = L_t(T, S, b_0, t_i_values)
##  log_likelihood = sum(log(arrayfun(@(t_i) pt(t_i, T, S, b_0), t_i_values)));
##endfunction
##
##% Discretización de las variables
##T_values = [70:140];  % Reducción a 20 puntos
##S_values = linspace(0.1, 5, 20); % Reducción a 20 puntos
##b0_values = linspace(0, 1, 20);  % Reducción a 20 puntos
##
##t_i_values = tiempo_picos - t_ini; % Ajustar los valores de t_i para que correspondan a los picos detectados
##
##max_likelihood = -Inf;
##best_T = T_values(1);
##best_S = S_values(1);
##best_b0 = b0_values(1);
##
##total_iterations = length(T_values) * length(S_values) * length(b0_values);
##iteration = 0;
##
##for T = T_values
##  for S = S_values
##    for b1 = b0_values
##      iteration++;
##      likelihood = L_t(T, S, b1, t_i_values);
##      if likelihood > max_likelihood
##        max_likelihood = likelihood;
##        best_T = T;
##        best_S = S;
##        best_b0 = b1;
##      endif
##      if mod(iteration, 100) == 0  % Imprimir progreso cada 100 iteraciones
##        disp(['Iteración ', num2str(iteration), ' de ', num2str(total_iterations)]);
##      endif
##    endfor
##  endfor
##endfor
##
##disp(['Mejor T: ', num2str(best_T)]);
##disp(['Mejor S: ', num2str(best_S)]);
##disp(['Mejor b1: ', num2str(best_b0)]);
##disp(['Max Likelihood: ', num2str(max_likelihood)]);

##% Calculo de p_e(e)
##t_values = linspace(-10, 10, 1000); % Asumiendo valores de tiempo para calcular p_t(t)
##T = best_T;
##S = best_S;
##b_0 = best_b0;
##
##% Calculo de p_t(t)
##p_t_values = arrayfun(@(t) pt(t, T, S, b_0), t_values);
##
##% Calculo de p_e(e) usando convolución
##p_e_values = conv(p_t_values, flip(p_t_values), 'same');
##
##% Graficar p_t(t) y p_e(e)
##figure(7);
##subplot(2, 1, 1);
##plot(t_values, p_t_values);
##xlabel('Tiempo');
##ylabel('p_t(t)');
##title('Función de densidad de probabilidad p_t(t)');
##
##subplot(2, 1, 2);
##plot(t_values, p_e_values);
##xlabel('e');
##ylabel('p_e(e)');
##title('Función de densidad de probabilidad p_e(e)');


