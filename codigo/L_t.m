function retval = L_t (T, S, b1, t_picos)
  periodo_beat = 60 / T;

  b2 = b1+periodo_beat;
  b3 = b2+periodo_beat;
  b4 = b3+periodo_beat;

  b1_4 = [b1 b2+periodo_beat*S b3 b4+periodo_beat*S];

  sum = 0;
  for i=1:length(t_picos)
    sum += log(p_t(t_picos(i), T, b1_4));
  endfor
  retval = sum;
endfunction
