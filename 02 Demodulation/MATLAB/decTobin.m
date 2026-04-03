function binaryMatrix = decTobin(decimalArray)
    % Convierte todos los números decimales en binario como cadenas de 2 dígitos
    binaryStrArray = dec2bin(decimalArray, 1);
    
    % Convierte las cadenas binarias en una matriz de números binarios
    binaryMatrix = double(binaryStrArray) - '0';
end
