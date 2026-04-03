function subvectores = sep(vector, n)
    % Obtener la longitud del vector original
    len = length(vector);
    
    % Calcular la cantidad de subvectores
    numSubvectores = ceil(len / n);
    
    % Calcular la cantidad de ceros que se deben agregar al final
    cerosFaltantes = n * numSubvectores - len;
    
    % Agregar ceros al final del vector original
    vector = [vector, zeros(1, cerosFaltantes)];
    
    % Inicializar la matriz de subvectores
    subvectores = zeros(numSubvectores, n);
    
    % Rellenar la matriz de subvectores
    for i = 1:numSubvectores
        startIdx = (i - 1) * n + 1;
        endIdx = i * n;
        subvectores(i, :) = vector(startIdx:endIdx);
    end
end
