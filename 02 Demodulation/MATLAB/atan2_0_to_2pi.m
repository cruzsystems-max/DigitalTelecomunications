function angles = atan2_0_to_2pi(y, x)
    angles = atan2(y, x);  % Calcula los ángulos en el rango de -π a π
    negative_angles = angles < 0;
    angles(negative_angles) = angles(negative_angles) + 2*pi;  % Ajusta los ángulos negativos
end