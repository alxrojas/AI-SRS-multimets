
%--------Funcion traslacion en z INICIO
function [des, mascara] = traslacion_z(dimensiones, dimCT, origen, estudio3D_CT, grados_num, mascaracontorno3D, isocentro)

%mascaracontorno3D2 = mascaracontorno3D;  % Por si no entra en el primer modificador y no se generara 
grosor_corte = dimCT(3);
dimens = size (mascaracontorno3D);
planos_empujar = abs ( round ( grados_num / grosor_corte));

if grados_num > 0
    % Desplazamos hacia la cabeza
    mascaracontorno3D2 = cat (3, zeros(dimens(1), dimens(2), planos_empujar), mascaracontorno3D(:, :, 1:dimens(3)-planos_empujar));
elseif grados_num < 0
    % Desplazamos hacia los pies
    mascaracontorno3D2 = cat (3, mascaracontorno3D(:, :, planos_empujar+1:end), zeros(dimens(1), dimens(2), planos_empujar));
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

%--------Funcion traslacion en z FIN
