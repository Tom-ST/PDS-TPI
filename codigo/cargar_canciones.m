function lista = cargar_canciones()
  archivos = dir("../samples/");
  for i=1:length(archivos)
    archivos(i).name;
  endfor
endfunction

