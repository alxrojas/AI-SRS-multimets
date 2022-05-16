

%--------Funcion traslacion en x INICIO
function [des, mascara] = traslacion_x(dimensiones, dimCT, origen, estudio3D_CT, grados_num, mascaracontorno3D, isocentro)

%mascaracontorno3D2 = mascaracontorno3D;  % Por si no entra en el primer modificador y no se generara 
anchopixel = dimCT(1);
dimens = size (mascaracontorno3D);
columnas_empujar = abs ( round ( grados_num / anchopixel));

if grados_num > 0
    % Desplazamos hacia la derecha
    mascaracontorno3D2 = cat (2, zeros(dimens(1), columnas_empujar, dimens(3)), mascaracontorno3D(:, 1:dimens(2)-columnas_empujar, :));
elseif grados_num < 0
    % Desplazamos hacia la izquierda
    mascaracontorno3D2 = cat (2, mascaracontorno3D(:, columnas_empujar+1:end, :), zeros(dimens(1), columnas_empujar, dimens(3)));
end

mascara = mascaracontorno3D2;
%Obtener el desplazamiento total de la primera traslacion en X
estructura_inicial = regionprops3(mascaracontorno3D, "Centroid");
centro = estructura_inicial.Centroid;
estructura_modificada = regionprops3(mascaracontorno3D2, "Centroid");
centro1 = estructura_modificada.Centroid;
V = centro1-centro;
V = V.* dimensiones;
des = sqrt(V*V');

%--------Funcion traslacion en x FIN
