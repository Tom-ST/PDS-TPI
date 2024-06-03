clc; clear all;
##variables (se pueden cambiar)
duracion = 0.5; %cantidad de minutos a leer del archivo



## Codigo
canciones = cargar_canciones();
seleccion = menu("Elige una cancion", canciones);

if seleccion > 0
  nombre_archivo = canciones{seleccion};
  disp(["Has seleccionado: ", nombre_archivo]);
else
  disp("No has seleccionado ninguna canción.");
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

%Convierto audio stereo a mono
if size(y ,2) == 2
  y = (y(:,1) + y(:,2)) / 2;
endif


plot(y);

%ventana para transformada de fourier de corto tiempo: 10ms (por que 10? se podria calcular este valor?)
%analizo la transformada de fourier en ventanas de 10 ms para detectar
% cambios en las frecuencias, solo utilizo la magnitud de la fft
% utilizo una funcion de compresion G(x) para que aquellos componentes de
% frecuencia alta (como platillos) no sean "masked(?)" por componentes
% de baja frecuencia pero alta amplitud.
% uso G(x) = x^(1/2), o arcsin(x). no uso log porque no se comporta bien cerca del cero

ventana = 0.001;% en segundos (10 ms)
G = @(x) arcsin(x);

