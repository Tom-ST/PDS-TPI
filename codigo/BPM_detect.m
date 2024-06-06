clc; clear all;
##variables (se pueden cambiar)
duracion = 0.1; %cantidad de minutos a leer del archivo

t_ini = 10; #Segundo que se comienza a analizar la cancion
t_fin = 15; #Termina el analisis

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


figure(1);
plot(t,y);
xlabel("Tiempo (s)");
ylabel("Amplitud");
title("Audio en funcion del tiempo")

%ventana para transformada de fourier de corto tiempo: 10ms (por que 10? se podria calcular este valor?)
%analizo la transformada de fourier en ventanas de 10 ms para detectar
% cambios en las frecuencias, solo utilizo la magnitud de la fft
% utilizo una funcion de compresion G(x) para que aquellos componentes de
% frecuencia alta (como platillos) no sean "masked(?)" por componentes
% de baja frecuencia pero alta amplitud.
% uso G(x) = x^(1/2), o arcsin(x). no uso log porque no se comporta bien cerca del cero

tam_ventana = 0.01;% en segundos (10 ms)
tam_ventana_m = fix(tam_ventana * Fs); %tamaño de la ventana en muestras
desplazamiento = fix(tam_ventana_m / 2); %desplazamiento de la ventana (50% overlap)

## Funcion de compresion
G = @(x) sqrt(x);

num_ventanas = fix(length(y) / tam_ventana_m);

fragmentos = zeros(num_ventanas, tam_ventana_m);

for i = 1:num_ventanas
  fragmentos(i,:) = y((i-1)*tam_ventana_m + 1 : i*tam_ventana_m);
endfor


fft_resultados = zeros(num_ventanas, tam_ventana_m);

for i = 1:num_ventanas
  fft_resultados(i, :) = G(abs(fft(fragmentos(i,:)))); #las filas son 10ms, hace la fft cada 10ms
endfor

figure(2)
stem(fft_resultados(100,:))
title("Transformada de fourier para una cierta fila")


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

convolucion=conv(y,y);
vector_conv=1:length(convolucion);
vector_conv = vector_conv/Fs;


figure(4)
plot(vector_conv,convolucion)
title("Convolucion")

%----- Calculo de los beats ------
% Calculo del escalar K
% Suponiendo que E_R(i) es una referencia conocida o un patrón que se espera encontrar
M = length(E);  % Horizonte en el cual se minimiza la distancia
E_R = rand(M, 1); % Ejemplo: generar una referencia aleatoria (reemplazar con la referencia correcta)

numerador = 0;
denominador = 0;

for i = 1:M
  numerador += E(i) * E_R(i);
  denominador += E_R(i) ^ 2;
endfor

K = numerador / denominador;

% Calculo de la norma L2 minimizada
L2_norm = 0;
for i = 1:M
  L2_norm += (E(i) - K * E_R(i)) ^ 2;
endfor


##
##% Calculo de la correlacion cruzada
##N_R = 100; % Numero de valores de R (ejemplo: 100 valores de R entre 60 y 150 BPM)
##tempos = linspace(60, 150, N_R);
##N_D = 10; % Numero de tiempos discretizados
##
##C = zeros(N_R, N_D);
##for r = 1:N_R
##  R = tempos(r);
##  T_r = 1 / (R / 60); % Periodo en segundos
##  for t = 1:N_D
##    j = (t-1) * T_r * Fs + 1;
##    if j + M - 1 <= length(E_R)
##      E_R_shifted = E_R(floor(j):floor(j) + M - 1);
##      C(r, t) = sum(E .* E_R_shifted);
##    end
##  end
##end
##
##[max_corr, idx] = max(C(:));
##[r_idx, t_idx] = ind2sub(size(C), idx);
##best_R = tempos(r_idx);
##best_t = (t_idx - 1) * (1 / (best_R / 60));
##
##disp(['Mejor tempo (BPM): ', num2str(best_R)]);
##disp(['Mejor candidato a downbeat (s): ', num2str(best_t)]);
##
##figure(4);
##imagesc(C);
##xlabel("Tiempos discretizados");
##ylabel("Tempos (BPM)");
##title("Correlación cruzada");
##colorbar;
