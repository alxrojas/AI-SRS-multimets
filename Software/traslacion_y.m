

%--------Funcion traslacion en y INICIO
function [des, mascara] = traslacion_y(dimensiones, dimCT, origen, estudio3D_CT, grados_num, mascaracontorno3D, isocentro)

%mascaracontorno3D2 = mascaracontorno3D;  % Por si no entra en el primer modificador y no se generara 
altopixel = dimCT(2);
dimens = size (mascaracontorno3D);
filas_empujar = abs ( round ( grados_num / altopixel));

if grados_num > 0
    % Desplazamos hacia el piso:
    mascaracontorno3D2 = cat (1, zeros(filas_empujar, dimens(2), dimens(3)), mascaracontorno3D(1:dimens(1)-filas_empujar, :, :));
elseif grados_num < 0 
    % Desplazamos hacia el techo
    mascaracontorno3D2 = cat (1, mascaracontorno3D(filas_empujar+1:end, :, :), zeros(filas_empujar, dimens(2), dimens(3)));
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

%--------Funcion traslacion en y FIN
