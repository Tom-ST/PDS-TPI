function [T,S,b1,likelihood] = mejor_T_S_b1 (T_values, S_values, t_ini, t_fin, tiempo_picos, T_aprox)
  ajuste = 20; % cuantos beats por encima y debajo del T_aprox utilizar

  max_likelihood = -Inf;
  best_T = nan;
  best_S = nan;
  best_b1 = nan;


  T_min = max(70, T_aprox - ajuste);
  T_max = min(140, T_aprox + ajuste);
  T_values = T_min:T_max;

  total_iterations = length(T_values) * length(S_values) * 32; %beats1 contiene 32 elementos

  iteration = 0;
  for i = 1:length(T_values)
    periodo_beat = 60 / T_values(i);
    b1_values = linspace(0, periodo_beat, 32)+t_ini;
    for j= 1:length(S_values)
      for k =1:length( b1_values)
        iteration++;
        likelihood = L_t(T_values(i), S_values(j), b1_values(k), tiempo_picos, t_fin);
        if likelihood > max_likelihood
          max_likelihood = likelihood;
          best_T = T_values(i)
          best_S = S_values(j)
          best_b1 = b1_values(k)
        endif
        if mod(iteration, 10) == 0
          disp(['Iteración ', num2str(iteration), ' de ', num2str(total_iterations)]);
        endif
      endfor
    endfor
  endfor

  T = best_T;
  S = best_S;
  b1 = best_b1;
  likelihood = max_likelihood;
endfunction
