clc; clear all;
##variables (se pueden cambiar)
duracion = 0.01; %cantidad de minutos a leer del archivo

##MUY BUENAS A TODOS

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
duracion_segundos = duracion * 60;
n_muestras = fix(duracion_segundos * info.SampleRate);
if info.TotalSamples < n_muestras
  error("El archivo es mas corto que la duracion especificada");
endif


[y, Fs] = audioread(archivo_audio, [1 n_muestras]);
t = (0:n_muestras-1) / Fs;

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

%ventana para transformada de fourier de corto tiempo: 10ms (por que 10? se podria calcular este valor?)
%analizo la transformada de fourier en ventanas de 10 ms para detectar
% cambios en las frecuencias, solo utilizo la magnitud de la fft
% utilizo una funcion de compresion G(x) para que aquellos componentes de
% frecuencia alta (como platillos) no sean "masked(?)" por componentes
% de baja frecuencia pero alta amplitud.
% uso G(x) = x^(1/2), o arcsin(x). no uso log porque no se comporta bien cerca del cero

tam_ventana = 0.001;% en segundos (10 ms)
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

figure
stem(fft_resultados(100,:))
title("Transformada de fourier para una cierta fila")


% Calculo del flujo de la energia
E_hat = zeros(num_ventanas - 1, 1); % Se resta 1 porque se calcula la diferencia entre frames sucesivos

for i = 2:num_ventanas
  E_hat(i-1) = sum(fft_resultados(i, :) - fft_resultados(i-1, :));
endfor

% Rectificacion media onda
E = max(E_hat, 0);

f_m = info.SampleRate #frecuencia de muestreo

figure;
plot(E);
xlabel("Frames");
ylabel("Flujo de energía");
title("Flujo de energía en función del tiempo");




##
##E_hat = zeros(1, num_ventanas - 1);
##for i = 2:num_ventanas
##  E_hat(i-1) = sum(fft_resultados(i, :) - fft_resultados(i-1, :));
##endfor
##
##E = max(E_hat, 0); % toma los valores positivos 2da ecuacion del paper
##
##% Graficar el flujo de energía
##tiempo_ventanas = (1:length(E)) * (desplazamiento / Fs);
##
##figure;
##plot(tiempo_ventanas, E);
##xlabel("Tiempo en segundos");
##ylabel("Flujo de Energía E(n)");
##
##


