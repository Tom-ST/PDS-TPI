function ER = expected_energy_flux(E, tempo, downbeat_location)
    % Parámetros de la función
    beats_per_sec = tempo / 60; % Conversión de BPM a beats por segundo
    beat_period = 1 / beats_per_sec; % Período del beat en segundos

    % Definir los momentos de tiempo para los pulsos principales, secundarios y terciarios
    main_pulse_times = downbeat_location + (0:beat_period:(length(E)-1)*beat_period);
    secondary_pulse_times = downbeat_location + (beat_period / 2):(beat_period):(length(E)-1)*beat_period;
    tertiary_pulse_times = downbeat_location + (beat_period / 4):(beat_period):(length(E)-1)*beat_period;

    % Generar las señales de pulsos
    main_pulses = zeros(size(E));
    secondary_pulses = zeros(size(E));
    tertiary_pulses = zeros(size(E));

    main_pulses(round(main_pulse_times * length(E))) = 1;
    secondary_pulses(round(secondary_pulse_times * length(E))) = 0.5; % Amplitud reducida para pulsos secundarios
    tertiary_pulses(round(tertiary_pulse_times * length(E))) = 0.25; % Amplitud reducida aún más para pulsos terciarios

    % Sumar los pulsos para obtener la señal de flujo de energía esperada
    ER = main_pulses + secondary_pulses + tertiary_pulses;
end

