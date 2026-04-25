%% Modulación M-PSK
% Script para generar señales M-PSK con visualización de formas de onda
% y constelación
% Autor: Cristian Cruz

clear all;
close all;
clc;

%% Parámetros
Num_intervalos = 5;

% Probabilidad de bits
p0 = 0.5;
p1 = 1 - p0;

% Parámetros de modulación M-PSK
M = 8;
n = log2(M);
Ts = 0.001;
fs = 1/Ts;

% Parámetros para la gráfica de señales
fmuestreo = 1000000;
fc = 10000;
Ac = 1;

% Número de bits de la secuencia a generar
N = Num_intervalos * n;

%% Generación de secuencia de bits
rng(42); % Semilla para reproducibilidad
secuencia = (rand(1, N) > p0);
fprintf('Secuencia generada: ');
fprintf('%d ', secuencia);
fprintf('\n');

%% Función para generar código Gray
function gray = gray_code(n)
    if n == 0
        gray = {''};
        return;
    end
    prev = gray_code(n-1);
    gray = cell(1, 2^n);
    for i = 1:length(prev)
        gray{i} = ['0' prev{i}];
    end
    for i = 1:length(prev)
        gray{length(prev) + i} = ['1' prev{length(prev) - i + 1}];
    end
end

%% Convertidor de Datos M-PSK
function [Ik, Qk, bk, tabla, gray_secuencia] = convertidor_datos_psk(bits, M, N_param)
    n = log2(M);

    % Asegurar múltiplo de n
    L = length(bits);
    if mod(L, n) ~= 0
        bits = bits(1:L - mod(L, n));
    end

    % Agrupar bits
    num_simbolos = length(bits) / n;
    grupos = reshape(bits, n, num_simbolos)';

    % Generar código Gray
    gray = gray_code(n);

    % Crear mapping
    mapping = containers.Map();
    for i = 1:M
        key = gray{i};
        mapping(key) = i - 1;
    end

    % Convertir grupos a índices b_k
    bk = zeros(1, num_simbolos);
    gray_secuencia = cell(1, num_simbolos);

    for i = 1:num_simbolos
        bits_str = '';
        for j = 1:n
            bits_str = [bits_str num2str(grupos(i,j))];
        end
        bk(i) = mapping(bits_str);
        gray_secuencia{i} = bits_str;
    end

    % Calcular fases: φ_k = π(2*b_k + N)/M
    phi_k = pi * (2*bk + N_param) / M;

    % Calcular componentes I_k y Q_k
    Ik = cos(phi_k);
    Qk = sin(phi_k);

    % Crear tabla de correspondencia
    tabla = struct();
    for i = 1:M
        code = gray{i};
        b = i - 1;
        phi = pi * (2*b + N_param) / M;
        I = cos(phi);
        Q = sin(phi);
        tabla(i).gray = code;
        tabla(i).bk = b;
        tabla(i).phi_k = phi;
        tabla(i).Ik = I;
        tabla(i).Qk = Q;
    end
end

% Convertir datos
[Ik, Qk, bk, tabla, gray_secuencia] = convertidor_datos_psk(secuencia, M, 1);

% Mostrar tabla de correspondencia
fprintf('\nTabla de correspondencia M-PSK:\n');
fprintf('%-10s %-5s %-15s %-15s %-15s\n', 'Gray', 'bk', 'phi_k', 'Ik=cos(phi_k)', 'Qk=sen(phi_k)');
fprintf('----------------------------------------------------------------------\n');
for i = 1:M
    fprintf('%-10s %-5d %-15.4f %-15.4f %-15.4f\n', ...
        tabla(i).gray, tabla(i).bk, tabla(i).phi_k, tabla(i).Ik, tabla(i).Qk);
end

%% Función pulso rectangular
function p = pulso_rectangular(t, Ts, fmuestreo)
    eps = 1/(2*fmuestreo);
    p = double((t >= -eps) & (t < Ts - eps));
end

%% Generar señales xi(t) y xq(t)
function [t, xi_t] = generar_xi(simbolos_I, Ts, fmuestreo)
    Ns = round(Ts * fmuestreo);
    t = (0:length(simbolos_I)*Ns-1) / fmuestreo;
    xi_t = zeros(size(t));

    for k = 1:length(simbolos_I)
        xi_t = xi_t + simbolos_I(k) * pulso_rectangular(t - (k-1)*Ts, Ts, fmuestreo);
    end
end

function [t, xq_t] = generar_xq(simbolos_Q, Ts, fmuestreo)
    Ns = round(Ts * fmuestreo);
    t = (0:length(simbolos_Q)*Ns-1) / fmuestreo;
    xq_t = zeros(size(t));

    for k = 1:length(simbolos_Q)
        xq_t = xq_t + simbolos_Q(k) * pulso_rectangular(t - (k-1)*Ts, Ts, fmuestreo);
    end
end

%% Generar señales para visualización
num_simbolos = 20;
[t, x_i] = generar_xi(Ik(1:min(num_simbolos, length(Ik))), Ts, fmuestreo);
[~, x_q] = generar_xq(Qk(1:min(num_simbolos, length(Qk))), Ts, fmuestreo);

% Pulso rectangular
pulso_rect = pulso_rectangular(t, Ts, fmuestreo);

% Portadoras
wc = 2 * pi * fc;
portadora_cos = cos(wc * t);
portadora_sin = -sin(wc * t);

% Productos de señales con portadoras
xi_cos = x_i .* portadora_cos;
xq_sin = -x_q .* sin(wc * t);

% Señal modulada
x_c = Ac * (x_i .* portadora_cos - x_q .* sin(wc * t));

fprintf('\nSeñales generadas exitosamente para %d símbolos\n', num_simbolos);
fprintf('Duración total: %.6f segundos\n', t(end));
fprintf('Número de muestras: %d\n', length(t));

%% Gráficas de las formas de onda
figure('Position', [100 50 1200 900]);

% 1. Pulso rectangular
subplot(8,1,1);
stairs(t, pulso_rect, 'LineWidth', 1.5, 'Color', [0 0.08 0.98]);
hold on;
area(t, pulso_rect, 'FaceColor', [0 0.08 0.98], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('Pulso Rectangular p_D(t) de duración Ts');
ylabel('Amplitud');
grid on;
xlim([0 t(end)]);

% 2. Señal en fase xi(t)
subplot(8,1,2);
stairs(t, x_i, 'LineWidth', 1.5, 'Color', [0.04 0.97 0]);
hold on;
area(t, x_i, 'FaceColor', [0.04 0.97 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('Señal en fase: xi(t)');
ylabel('Amplitud');
grid on;
xlim([0 t(end)]);
% Añadir códigos Gray
for k = 1:min(length(gray_secuencia), num_simbolos)
    if k <= length(Ik)
        t_pos = (k - 0.5) * Ts;
        if Ik(k) >= 0
            y_pos = Ik(k) + 0.15;
        else
            y_pos = Ik(k) - 0.2;
        end
        text(t_pos, y_pos, gray_secuencia{k}, 'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 8);
    end
end

% 3. Señal en cuadratura xq(t)
subplot(8,1,3);
stairs(t, x_q, 'LineWidth', 1.5, 'Color', [0.04 0.97 0]);
hold on;
area(t, x_q, 'FaceColor', [0.04 0.97 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('Señal en cuadratura: xq(t)');
ylabel('Amplitud');
grid on;
xlim([0 t(end)]);
% Añadir códigos Gray
for k = 1:min(length(gray_secuencia), num_simbolos)
    if k <= length(Qk)
        t_pos = (k - 0.5) * Ts;
        if Qk(k) >= 0
            y_pos = Qk(k) + 0.15;
        else
            y_pos = Qk(k) - 0.2;
        end
        text(t_pos, y_pos, gray_secuencia{k}, 'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 8);
    end
end

% 4. Portadoras
subplot(8,1,4);
plot(t, portadora_cos, 'b', 'LineWidth', 1.2);
hold on;
plot(t, portadora_sin, 'r', 'LineWidth', 1.2);
title('Portadoras: cos(wc·t) y -sin(wc·t)');
ylabel('Amplitud');
grid on;
legend('cos(wc·t)', '-sin(wc·t)', 'Location', 'northeast');
xlim([0 t(end)]);

% 5. Producto xi(t)cos(wc·t)
subplot(8,1,5);
plot(t, xi_cos, 'b', 'LineWidth', 0.8);
title('Producto: xi(t)·cos(wc·t)');
ylabel('Amplitud');
grid on;
xlim([0 t(end)]);

% 6. Producto -xq(t)sin(wc·t)
subplot(8,1,6);
plot(t, xq_sin, 'r', 'LineWidth', 0.8);
title('Producto: -xq(t)·sin(wc·t)');
ylabel('Amplitud');
grid on;
xlim([0 t(end)]);

% 7. Componentes de la señal modulada
subplot(8,1,7);
plot(t, xi_cos, 'b', 'LineWidth', 0.8);
hold on;
plot(t, xq_sin, 'r', 'LineWidth', 0.8);
title('Componentes de la señal modulada');
ylabel('Amplitud');
legend('xi(t)·cos(wc·t)', '-xq(t)·sin(wc·t)', 'Location', 'northeast', 'FontSize', 8);
xlim([0 t(end)]);

% 8. Señal modulada
subplot(8,1,8);
plot(t, xi_cos + xq_sin, 'k', 'LineWidth', 1.2);
hold on;
title(sprintf('Señal Modulada %d-PSK', M));
ylabel('Amplitud');
legend('Suma', 'Location', 'northeast', 'FontSize', 8);
xlim([0 t(end)]);
% Añadir códigos Gray
y_max = max(xi_cos + xq_sin);
y_min = min(xi_cos + xq_sin);
y_range = y_max - y_min;
for k = 1:min(length(gray_secuencia), num_simbolos)
    t_pos = (k - 0.5) * Ts;
    y_pos = y_min - 0.15 * y_range;
    text(t_pos, y_pos, gray_secuencia{k}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'yellow', 'EdgeColor', 'black', 'FontSize', 8);
end

%% Gráfica de Constelación M-PSK
figure('Position', [150 100 800 800]);

% Calcular todos los puntos de la constelación
I_constellation = zeros(1, M);
Q_constellation = zeros(1, M);
labels_constellation = cell(1, M);

for i = 1:M
    I_constellation(i) = tabla(i).Ik;
    Q_constellation(i) = tabla(i).Qk;
    labels_constellation{i} = tabla(i).gray;
end

% Graficar constelación
plot(I_constellation, Q_constellation, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;

% Añadir etiquetas con los códigos Gray
for i = 1:M
    % Calcular offset para la etiqueta
    offset_factor = 0.15;
    I_offset = I_constellation(i) * offset_factor;
    Q_offset = Q_constellation(i) * offset_factor;

    text(I_constellation(i) + I_offset, Q_constellation(i) + Q_offset, ...
        labels_constellation{i}, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10, 'FontWeight', 'bold');
end

% Dibujar círculo unitario
theta_circle = linspace(0, 2*pi, 100);
plot(cos(theta_circle), sin(theta_circle), 'b--', 'LineWidth', 1);

% Dibujar ejes
plot([-1.5 1.5], [0 0], 'k-', 'LineWidth', 0.5);
plot([0 0], [-1.5 1.5], 'k-', 'LineWidth', 0.5);

% Configuración de la gráfica
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
xlabel('En fase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cuadratura (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Constelación %d-PSK', M), 'FontSize', 14, 'FontWeight', 'bold');

% Añadir información adicional
text(-1.4, 1.4, sprintf('M = %d\nN = %d', M, 1), 'FontSize', 10, ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');

fprintf('\n¡Gráficas generadas exitosamente!\n');
