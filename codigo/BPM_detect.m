archivo_audio = '../samples/Canon_in_D.mp3';



##variables (se pueden cambiar)
archivo = 1; %numero de archivo, del 1 al 10
duracion = 1; %cantidad de minutos a leer del archivo







% Cargar la cantidad de minutos especificados
% falta chequear que pasa cuando el archivo es mas corto que la duracion requerida
info = audioinfo(archivo_audio);
duracion_segundos = duracion * 60;
n_muestras = duracion_segundos * info.SampleRate;

[y, Fs] = audioread(archivo_audio, [1 n_muestras]);

%Convierto audio stereo a mono
if size(y ,2) == 2
  y = (y(:,1) + y(:,2)) / 2;
endif

plot(y)
