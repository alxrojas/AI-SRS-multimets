
function [des, mascara] = rotation_pitch(dimensiones, dimCT, origen, estudio3D_CT, grados_num, mascaracontorno3D, isocentro)
x0 = origen(1);
y0 = origen(2);
z0 = origen(3);
[filas, columnas, planos] = size(mascaracontorno3D);
IsocentroX = isocentro(1);
IsocentroY = isocentro(2);
IsocentroZ = isocentro(3);
anchopixelCT = dimCT(1);
altopixelCT = dimCT(2);
grosorcorteCT = dimCT(3);
columnaIso = round ( abs( IsocentroX - x0)/anchopixelCT);
filaIso = round ( abs( IsocentroY - y0)/altopixelCT);
planoIso = round ( abs( IsocentroZ - z0)/grosorcorteCT);
colsDer = columnas - columnaIso;
columnas_completar = abs ( 2 * columnaIso - 1 - columnas);
filas_completar = abs ( 2 * filaIso - 1 - filas);
planos_completar = abs ( 2 * planoIso - 1 - planos);
planos_completar;
if columnaIso > columnas/2   % Completar a la derecha
    mascaracontorno3D20 = cat (2, mascaracontorno3D, zeros (filas, columnas_completar, planos));
else % Completar a la izquierda
    mascaracontorno3D20 = cat (2, zeros(filas, columnas_completar, planos), mascaracontorno3D); 
end
nuevascolumnas = columnas + columnas_completar;

if filaIso > filas/2   % Completar por abajo
    mascaracontorno3D20 = cat ( 1, mascaracontorno3D20, zeros( filas_completar, nuevascolumnas, planos));
else   % Completar arriba
    mascaracontorno3D20 = cat ( 1, zeros(filas_completar, nuevascolumnas, planos), mascaracontorno3D20);
end
nuevasfilas = filas + filas_completar;

if planoIso > planos/2   % Completar desde la direccion cabeza
    mascaracontorno3D20 = cat ( 3, mascaracontorno3D20, zeros( nuevasfilas, nuevascolumnas, planos_completar));
else    % Completar por el lado de los pies
    mascaracontorno3D20 = cat ( 3, zeros(nuevasfilas, nuevascolumnas, planos_completar), mascaracontorno3D20);
end
i = 1;
terminator = 0;
while terminator == 0
    planoRoiExploratorio = mascaracontorno3D20 (:, i, :);
    if any ( planoRoiExploratorio(:))   % Entra en este if, si hay algun valor distinto de cero
        primer_plano_pitcheable = i;
        terminator = 1;
    end
    i = i + 1;
end
terminator = 0;
while terminator == 0
    nuevoPlanoRoiExploratorio = mascaracontorno3D20 (:, i, :);
    if all (nuevoPlanoRoiExploratorio(:) == 0)
        ultimo_plano_pitcheable = i-1;
        terminator = 1;
    end
    i = i + 1;
end
total_planos_pitchear = ultimo_plano_pitcheable - primer_plano_pitcheable + 1;
plano_medio = round (0.5 * (primer_plano_pitcheable + ultimo_plano_pitcheable));
for i = primer_plano_pitcheable : ultimo_plano_pitcheable
    plano_imrotatear = mascaracontorno3D20 (:, i, :);
    plano_imrotatear = permute ( plano_imrotatear, [1 3 2 ]);
    plano_imrotateado = imrotate ( plano_imrotatear, grados_num, 'crop');
    plano_insertar = permute ( plano_imrotateado, [ 1 3 2]);
    mascaracontorno3D20 (:, i, :) = plano_insertar;
    pause (0.05)
end

% Ya rotó en la dimension que fuera, ahora hay que quitar los pedazos que se añadieron
if columnaIso > columnas/2  % Se completo de la derecha y hay que quitar de ah
    mascaracontorno3D20 = mascaracontorno3D20 (:, 1:end-columnas_completar, :);
else % Se completo de la izquierda 
    mascaracontorno3D20 = mascaracontorno3D20 (:, columnas_completar+1:end, :);
end
if filaIso > filas/2 % puesta y quita de abajo
    mascaracontorno3D20 = mascaracontorno3D20 (1:end-filas_completar, :, :);
else % de arriba
    mascaracontorno3D20 = mascaracontorno3D20 (filas_completar+1:end, :, :);
end
if planoIso > planos/2  % de cabeza
    mascaracontorno3D20 = mascaracontorno3D20 (:, :, 1:end-planos_completar);
else  % de pies
    mascaracontorno3D20 = mascaracontorno3D20 (:, :, planos_completar+1:end);
end

mascara = mascaracontorno3D20;
% obtiene el desplazamiento total de la primera rotacion en ROLL
estructura_inicial = regionprops3(mascaracontorno3D, "Centroid");
centro = estructura_inicial.Centroid;
estructura_modificada = regionprops3(mascaracontorno3D20, "Centroid");
centro1 =  estructura_modificada.Centroid;
V = centro1 - centro;
V = V .* dimensiones;
des = sqrt(V*V');
