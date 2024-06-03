function lista = cargar_canciones()
  archivos = dir("../samples/");
  lista = {};

  for i=3:length(archivos)
    lista(end+1) = archivos(i).name;
  endfor
endfunction

