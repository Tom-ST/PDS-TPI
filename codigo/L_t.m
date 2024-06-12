function retval = L_t (T, S, b1, t_picos,t_ini,t_fin)
  periodo_beat = 60 / T;

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


  sum = 0;
  for i=1:length(t_picos)
    sum += log(p_t(t_picos(i), T, b1_n));
  endfor
  retval = sum;

endfunction
