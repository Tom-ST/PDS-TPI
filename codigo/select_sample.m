function filename = select_sample(number)
  switch (number)
    case 1
      filename = '../samples/Canon_in_D.mp3';
	otherwise
		error('No existe ese archivo de audio');
  endswitch
  display(strcat('Archivo seleccionado: ',filename(12:end)));
endfunction
