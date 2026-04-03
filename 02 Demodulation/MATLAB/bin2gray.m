function gray_code = bin2gray(binary_matrix)
    % Convierte una matriz de números binarios a código Gray
    
    % Verifica si la entrada es una matriz binaria
    if ~ismatrix(binary_matrix) || any(binary_matrix(:) ~= 0 & binary_matrix(:) ~= 1)
        error('La entrada debe ser una matriz binaria (conteniendo solo 0s y 1s).');
    end
    
    [m, n] = size(binary_matrix);
    gray_code = false(m, n);
    gray_code(:, 1) = binary_matrix(:, 1);
    
    for j = 2:n
        gray_code(:, j) = xor(binary_matrix(:, j-1), binary_matrix(:, j));
    end
end

