pkg load signal
clc; clear all; close all;
##variables (se pueden cambiar)
t_ini = 50; #Segundo que se comienza a analizar la cancion
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

duracion_segundos = t_fin - t_ini;

muestra_ini = fix(t_ini * info.SampleRate);
muestra_fin = fix(t_fin * info.SampleRate);

if or(info.TotalSamples < muestra_ini, info.TotalSamples < muestra_fin)
  disp(["El rango de tiempo seleccionado esta por fuera de la duracion del archivo: ", num2str(info.Duration), "s"]);
  return
endif


[y, Fs] = audioread(archivo_audio, [muestra_ini muestra_fin]);
t = t_ini:1/Fs:t_fin;

%Convierto audio stereo a mono
#Audio_read nos da dos columnas, una para el sonido q sale por el parlante izquierdo
#y otra que sale por el parlante derecho, hacemos un promedio
if size(y ,2) == 2
  y = (y(:,1) + y(:,2)) / 2;
endif


figure();
plot(t,y);
xlabel("Tiempo (s)");
ylabel("Amplitud");
title("Audio en funcion del tiempo")

% analizo la transformada de fourier en ventanas de 10 ms para detectar
% cambios en las frecuencias, solo utilizo la magnitud de la fft
% utilizo una funcion de compresion C(x) para que aquellos componentes de
% frecuencia alta (como platillos) no sean enmascarados por componentes
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

figure();
subplot(2,1,1);
hold on;
title("Transformada de fourier para una cierta ventana de 10ms");
stem(abs(fft(fragmentos(100,:))));

subplot(2,1,2)
hold on;
title("Transformada de fourier para la misma ventana de 10ms, con funcion de compresion");
stem(fft_resultados(100,:));

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

figure(3);
##plot(0.01:0.01:t(end),E);
plot((t_ini+0.01):0.01:(t_fin-0.01),E);

xlabel("Tiempo");
ylabel("Flujo de energía");
title("Flujo de energía en función del tiempo");



% =========================== Identificacion de Picos Significativos ===============================
[peaks, peak_locs] = findpeaks(E);

% Ploteamos los picos en la señal de flujo de energía
figure();
subplot(2,1,1)
plot((t_ini+0.01):0.01:(t_fin-0.01),E);
hold on;
plot((t_ini+0.01)+peak_locs*0.01, peaks, 'ro'); % Convertimos los índices a tiempo
hold off;
xlabel("Tiempo");
ylabel("Flujo de energía");
title("Picos de flujo de energía en función del tiempo");
legend('Flujo de energía', 'Picos');

% Establecer un umbral como un porcentaje del máximo pico encontrado
umbral = 0.5; % Por ejemplo, seleccionamos picos que estén por encima del 50% del máximo pico

% Calculamos el máximo pico
max_peak = max(peaks);

% Filtramos los picos significativos que superan el umbral
picos_significativos = peak_locs(peaks >= umbral * max_peak);
valores_significativos = peaks(peaks >= umbral * max_peak);

% Ploteamos los picos significativos en la señal de flujo de energía
subplot(2,1,2);
plot((t_ini+0.01):0.01:(t_fin-0.01),E);
hold on;
plot((t_ini+0.01)+picos_significativos*0.01, valores_significativos, 'ro'); % Convertimos los índices a tiempo
% Agregar la línea horizontal intermitente del umbral
y_val = umbral * max_peak;
plot([t_ini, t_fin], [y_val, y_val], 'r--'); % Línea horizontal intermitente en el valor del umbral
hold off;

xlabel("Tiempo");
ylabel("Flujo de energía");
title("Picos de flujo de energía significativos en función del tiempo");
legend('Flujo de energía', 'Picos significativos', 'Umbral');

% Ploteamos la señal de audio en función del tiempo
figure();
plot(t, y);
hold on;

% Ploteamos líneas verticales en las ubicaciones de los picos
tiempo_picos = (t_ini + 0.01) + picos_significativos * 0.01; % Convertimos los índices de picos a tiempo
for i = 1:length(tiempo_picos)
    line([tiempo_picos(i), tiempo_picos(i)], ylim, 'Color', 'r', 'LineWidth', 1.5); % Graficamos una línea vertical en la posición del pico
end

hold off;

xlabel("Tiempo (s)");
ylabel("Amplitud");
title("Picos identificados en la señal de audio");



%============================ Establecer rangos de T, S y b1 ==================
Tempos = [70:140];
Swings = [0:0.1:0.4]; %El swing (S) es un porcentaje que determina que tanto se atrasa el segundo y cuarto cuarto-beat (es una propiedad de algunos generos musicales como rock y jazz)

% ==== A continuacion codigo para debuggear ==
##periodo_beat = 60 / T;
##beats1 = linspace(0, periodo_beat, 32)+t_ini;

##T = 140;
##S = Swings(1);
##b1 = beats1(1);
##b1 = tiempo_picos(1);

##b1_n = b1;  % Inicializamos el primer beat
##beat = b1 + periodo_beat;  % Inicializamos el siguiente beat
##
##% Generamos los beats adicionales hasta alcanzar o superar el tiempo t_fin
##while beat <= t_fin
##    b1_n = [b1_n, beat];
##    beat = beat + periodo_beat;  % Calculamos el siguiente beat
##end
##
##% Ajustamos los beats en las segundas y cuartas posiciones
##for i = 2:2:length(b1_n)
##    b1_n(i) = b1_n(i) + periodo_beat * S;
##end
##
##for i = 4:4:length(b1_n)
##    b1_n(i) = b1_n(i) + periodo_beat * S;
##end
##
####t=linspace(t_ini, t_fin,10000);
##pt = p_t (t, T, b1_n);
##
##
####figure();
##hold on
##plot(t,pt)
##title("Posicion probable donde se encuentran los cuarto-beats");
##xlabel("Tiempo");

% ==== Fin codigo para debugear


%-------MAYOR LIKELIHOOD-----------------
[T,S,b1,likelihood] = mejor_T_S_b1 (Tempos, 0, t_ini, t_fin, tiempo_picos)
%====================== Fin

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
