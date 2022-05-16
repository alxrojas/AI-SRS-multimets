function varargout = VisorIPRO(varargin)
% VISORIPRO MATLAB code for VisorIPRO.fig
%      VISORIPRO, by itself, creates a new VISORIPRO or raises the existing
%      singleton*.
%
%      H = VISORIPRO returns the handle to a new VISORIPRO or the handle to
%      the existing singleton*.
%
%      VISORIPRO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISORIPRO.M with the given input arguments.
%
%      VISORIPRO('Property','Value',...) creates a new VISORIPRO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisorIPRO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisorIPRO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisorIPRO

% Last Modified by GUIDE v2.5 25-Apr-2022 09:35:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VisorIPRO_OpeningFcn, ...
                   'gui_OutputFcn',  @VisorIPRO_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VisorIPRO is made visible.
function VisorIPRO_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for VisorIPRO
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes3)
logo = imread('LogoZunino.png');
image(logo)
axis off
axis image

% --- Outputs from this function are returned to the command line.
function varargout = VisorIPRO_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

function pushbutton15_Callback(hObject, eventdata, handles)
%Search all posible structures, read file paths
cadenaString = uigetdir();
f = waitbar(0, 'Loading files');
pause(0.5)
[namestructure, pathdose] = ruta (cadenaString, 1);
[namedose, pathdose] = ruta (cadenaString, 2);
[namect,pathct] = ruta(cadenaString,3);
[nameplan, pathdose] = ruta (cadenaString, 4);
cd (pathdose);
waitbar(0.25,f, 'Loading files');
pause(0.5)
RS = dicominfo (namestructure)
name_patient = RS.PatientName;
family = name_patient.FamilyName;
aux_isocenter = dicominfo(nameplan);
isocenter = aux_isocenter.BeamSequence.Item_1.ControlPointSequence.Item_1.IsocenterPosition;
waitbar(0.5,f, 'Reading DICOM files');
pause(0.5)
RS.StructureSetROISequence.Item_10;
%for i=1:length(fieldnames(RS.StructureSetROISequence))
for i=1:length(fieldnames(RS.ROIContourSequence))
    %eval (['nombre{i} = RS.StructureSetROISequence.Item_', num2str(i);])
    eval (['nombre{i} = RS.StructureSetROISequence.Item_', num2str(i), '.ROIName;'])
    numero{i} = i;
end
waitbar(0.75,f, 'Reading structures');
pause(0.5)
handles.nombre = nombre;
handles.numero = numero;
handles.cadenaString = cadenaString;
handles.pathdose = pathdose;
handles.namestructure = namestructure;
handles.pathct = pathct;
handles.namect = namect;
handles.namedose = namedose;
handles.nameplan = nameplan;
handles.isocenter = isocenter;
waitbar(1,f, 'Finishing');
pause(0.5)
close(f);
guidata (hObject, handles)
set(handles.popupmenu2,'string',nombre)
set (handles.edit56, 'String', sprintf(family))

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)

set(hObject, 'String', handles.nombre)
lista = get(hObject, 'String');
region = get(hObject, 'Value');
handles.region = region;
guidata (hObject, handles)

function popupmenu2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton1_Callback(hObject, eventdata, handles)
% Place CT on the interface
axes ( handles.axes1)
cadenaString = handles.cadenaString;
pathdose = handles.pathdose;
namestructure = handles.namestructure;
namedose = handles.namedose;
pathct = handles.pathct;
namect = handles.namect;
estudio4D = dicomreadVolume (pathct);
estudio3D = squeeze ( estudio4D);
estudio3D = double ( estudio3D);
dimens = size (estudio3D);
corte_test = estudio3D (:, :, 100);
imagesc (corte_test)
colormap (gray)
set (handles.edit1, 'String', '100');
arreglozeros = zeros (size (dimens));
filas = dimens (1);
columnas = dimens (2);
planos  = dimens (3);
cd (pathct);
dir_carpeta = dir (pathct);
elementos = length (dir_carpeta);
% Ghost directories: there are 2 (WINDOWS) and 3 (MAC)
fichero3 = dir_carpeta(3).name;
primerC = fichero3 (1);
if primerC == 'C'
    fantasmas = 2;
else
    fantasmas = 3;
end
primer_Fich = dir_carpeta( fantasmas+1).name;
info_primerFich = dicominfo (primer_Fich);
ImagePosPat = info_primerFich.ImagePositionPatient;
ubicacion = ImagePosPat (3);
ubicaciones (1) = ubicacion;
for i = fantasmas+2 : fantasmas+planos
    ficheroi = dir_carpeta(i).name;
    info_ficheroi = dicominfo ( ficheroi);
    ImagePositionPatienti = info_ficheroi.ImagePositionPatient;
    ubicaciones (i - fantasmas) = ImagePositionPatienti ( 3);
end
diferencias = 1000 * ones (planos-1);
diferencias = abs ( ubicaciones(2:end) - ubicacion);
grosor_corte = min (diferencias);
slicelocationmax = ImagePosPat (3);
slicelocationmin = ImagePosPat (3);
for i = fantasmas+2 : planos
    fichero = dir_carpeta(i).name;
    infofichero = dicominfo (fichero);
    ImagePosPat = infofichero.ImagePositionPatient;
    slicelocation = ImagePosPat ( 3);
    if slicelocation > slicelocationmax
        slicelocationmax = slicelocation;
    elseif slicelocation < slicelocationmin
        slicelocationmin = slicelocation;
    end
end
% Final slice of the study for locating xo, yo, zo
cd (pathct);
ficherio = dir (pathct); 
terminator = 0;
equis0 = ImagePosPat ( 1);
Ye0 = ImagePosPat ( 2);
anchopixel = info_primerFich.PixelSpacing ( 2);
altopixel = info_primerFich.PixelSpacing ( 1);
cd (pathdose);
infostr = dicominfo (namestructure);
numero_ROI_str = num2str ( get ( handles.popupmenu2, 'Value'));
nombre_estructura = eval (['infostr.StructureSetROISequence.Item_', numero_ROI_str, '.ROIName']);
set ( handles.edit2, 'String', nombre_estructura)
i = 1;
while terminator == 0
    try
        %contorno_cortei{i} = eval (['infostr.StructureSetROISequence.Item_', numero_ROI_str,'.ContourSequence.Item_', ...
        contorno_cortei{i} = eval (['infostr.ROIContourSequence.Item_', numero_ROI_str,'.ContourSequence.Item_', ...
              num2str(i), '.ContourData;']);
    catch
        terminator = 1;
        cantidad_islas = i - 1;
    end
    i = i + 1;
end                 % Contorno_cortei{i} contains isle contour i of ROI selected
planozeros = zeros (filas, columnas);
i = 1;
mascaracontorno3D = zeros ( dimens);
%offset
zplanoroi1 = contorno_cortei{1}(3);
planoroi1 = 1 + round ( ( zplanoroi1 - slicelocationmin)/ grosor_corte);    % Plane/CT space where the isle is
while i < cantidad_islas
    j = 0;
    contorno_en_uso = contorno_cortei{i};
    contorno_siguiente = contorno_cortei{i+1};
    largo = length (contorno_en_uso);
    zeta = contorno_en_uso( 3);
    zeta_Siguiente = contorno_siguiente( 3);
    deltazeta = zeta - zplanoroi1;
    deltaplanos = deltazeta / grosor_corte;
    planoi = planoroi1 + deltaplanos;
    Equises = contorno_en_uso (1 : 3 : largo-2);
    Yes = contorno_en_uso (2:3:largo-1);
    Columnas = round ( (Equises - equis0)/ anchopixel);
    Filas = round ( (Yes - Ye0) / altopixel);
    mascaracontornoi = roipoly (planozeros, Columnas, Filas);
    mascaracontornoi = double (mascaracontornoi);
    if zeta == zeta_Siguiente
        largo = length (contorno_siguiente);
        Equises_siguientes = contorno_siguiente (1 : 3 : largo-2);
        Columnas_siguientes = round ( (Equises_siguientes - equis0)/ anchopixel);
        Yes_siguientes = contorno_siguiente (2:3:largo-1);
        Filas_siguientes = round ( (Yes_siguientes - Ye0) / altopixel);
        mascaracontornoi_siguiente = roipoly (planozeros, Columnas_siguientes, Filas_siguientes);
        mascaracontornoi_siguiente = double (mascaracontornoi_siguiente);
        mascaracontornoi = mascaracontornoi + mascaracontornoi_siguiente;
        j = 1;
    end
    set (handles.edit1, 'String', num2str (planoi))
    pause (0.1)
    planoenmasc = estudio3D(:, :, planoi) .* mascaracontornoi;
    planodobleimagen = planoenmasc + estudio3D(:, :, planoi);
    imagesc ( planodobleimagen)
    axis off
    mascaracontorno3D(:, :, planoi) = mascaracontornoi;
    i = i + 1 + j;
end
imagesc ( planodobleimagen);
axis off
archivo_dosis = dicomread (namedose);
info_dosis = dicominfo (namedose);
dosis3D = squeeze ( double ( archivo_dosis * info_dosis.DoseGridScaling));
dimens = size (dosis3D);
anchopixel_dosis = info_dosis.PixelSpacing(2);
altopixel_dosis = info_dosis.PixelSpacing(1);
origen_dosis = info_dosis.ImagePositionPatient;
estudioCTmasROI = estudio3D + estudio3D .* mascaracontorno3D; 
% GENERATE DVH FROM MATLAB
estudio4D = dicomreadVolume (pathct);
estudio3D = squeeze ( estudio4D);
estudio3D = double ( estudio3D);
dimensCT = size (estudio3D);
corte_test = estudio3D (:, :, 100);
arreglozeros = zeros (size (dimensCT));
filas = dimensCT ( 1);
columnas = dimensCT ( 2);
planos  = dimensCT ( 3);
cd (pathct);
dir_carpeta = dir (pathct);
elementos = length (dir_carpeta);
h = 1;
while h < elementos
    p = dir_carpeta(h).name;
    pC = p(1);
    if pC == 'C'
        init = h;
        h = elementos + 1;
    end
    h = h + 1;
end
fich4 = dir_carpeta(init+1).name;
infofich4 = dicominfo (fich4)
infofich4
slicelocationmax = infofich4.SliceLocation;
%slicelocationmax = 81;
%slicelocationmin = 81;
slicelocationmin = infofich4.SliceLocation;
infofichero
for i = init : planos
    fichero = dir_carpeta(i).name;
    infofichero = dicominfo (fichero);
    if infofichero.SliceLocation > slicelocationmax
        slicelocationmax = infofichero.SliceLocation;
        numfichmax = i;
    elseif infofichero.SliceLocation < slicelocationmin
        slicelocationmin = infofichero.SliceLocation;
       numfichmin = i;
    end
end
fichmax = dir_carpeta ( numfichmax).name;
infomax = dicominfo ( fichmax);
esquinamax = infomax.ImagePositionPatient;
fichmin = dir_carpeta ( numfichmin).name;
infomin = dicominfo ( fichmin);
esquinamin = infomin.ImagePositionPatient;
dimspixelCT = infofich4.PixelSpacing;
anchopixelCT = dimspixelCT ( 2);
altopixelCT = dimspixelCT ( 1);
grosorcorteCT = infofich4.SliceThickness;
cd (pathdose);
RDdosis = dicomread (namedose);
RDdosis = double ( squeeze (RDdosis));
RDinfo = dicominfo (namedose);
plan_name = RDinfo.SeriesDescription;
UnidadesDosis = RDinfo.DoseUnits;
FactorEscala = RDinfo.DoseGridScaling;
DimsPixelsD = RDinfo.PixelSpacing;
anchopixelD = DimsPixelsD ( 2);
altopixelD = DimsPixelsD ( 1);
Filas = RDinfo.Rows;
Columnas = RDinfo.Columns;
OrigenDosis = RDinfo.ImagePositionPatient;
PlanosDosis = RDinfo.NumberOfFrames;
VectorOffsetPlanos = RDinfo.GridFrameOffsetVector;  % Z for each plane
grosorSliceD = abs( VectorOffsetPlanos( 2) - VectorOffsetPlanos( 1));
vectorxdose = double ( anchopixelD * [ 0 : Columnas-1]);% Interpolation. vectors on dose space
vectorydose = double ( altopixelD * [ 0 : Filas-1]);
vectorzdose = double ( grosorSliceD * [ 0 : PlanosDosis-1]);
[Xd, Yd, Zd] = meshgrid ( vectorxdose, vectorydose, vectorzdose);
vectorXD_ect = double ( anchopixelD * [0 : anchopixelCT/anchopixelD : double(Columnas-1)]); % Dose vector for CT
vectorYD_ect = double ( altopixelD * [0 : altopixelCT/altopixelD : double(Filas-1)]);
vectorZD_ect = double ( grosorSliceD * [0 : grosorcorteCT/grosorSliceD : double(PlanosDosis-1)]);
[Xect, Yect, Zect] = meshgrid ( vectorXD_ect, vectorYD_ect, vectorZD_ect);
Dosis3D_ect = interp3 ( Xd, Yd, Zd, RDdosis, Xect, Yect, Zect, 'linear'); % Dose on CT
Xod = OrigenDosis ( 1);
Yod = OrigenDosis ( 2);
Zod = OrigenDosis ( 3);
Xoct = esquinamin ( 1);
Yoct = esquinamin ( 2);
Zoct = esquinamin ( 3) - grosorcorteCT;  % Sup temporal: el origen de la CT es en la cara superior del plano y el origen de la dosis en la cara inferior. 
columnasD_ect = int16 ( length (vectorXD_ect));
filasD_ect = int16 ( length ( vectorYD_ect));
planosD_ect = int16 ( length ( vectorZD_ect));
planoi1 = int16 (round (abs( Zod - Zoct) / grosorcorteCT));
if planoi1 == 0
    planoi1 = int16 (round (1 + abs( Zod - Zoct) / grosorcorteCT));
end
columnai1 = int16 ( round (abs( Xod - Xoct) / anchopixelCT));
if columnai1 == 0 
    columnai1 = int16 ( round (1+ abs( Xod - Xoct) / anchopixelCT));
end
filai1 = int16 ( round ( abs( Yod - Yoct) / altopixelCT));
if filai1 == 0
    filai1 = int16 ( round (1 + abs( Yod - Yoct) / altopixelCT));
end
BaseDosisInterp = zeros ( dimensCT);
%OJO
BaseDosisInterp ( filai1:filai1+filasD_ect-1, columnai1:columnai1+columnasD_ect-1, planoi1:planosD_ect-2) = Dosis3D_ect(:, :, 3:end);
u = size(mascaracontorno3D);
%axes ( handles.axes2)
%imagesc ( BaseDosisInterp (:, :, 100))
%axis off
%colormap (handles.axes2, hsv)
BaseDosisInterp = imresize3(BaseDosisInterp, [u(1) u(2) u(3)], 'Linear');
%%%%
%%%%
set ( handles.edit1, 'String', 100)
Dosis_solo_en_ROI = BaseDosisInterp .* mascaracontorno3D;
Dosis_solo_en_ROI = FactorEscala * Dosis_solo_en_ROI ( Dosis_solo_en_ROI>0);
%Dosis_solo_en_ROI = Dosis_solo_en_ROI ( Dosis_solo_en_ROI>0);
CantidadV = length ( Dosis_solo_en_ROI);
VolumenROI = CantidadV * anchopixelCT * altopixelCT * grosorcorteCT;
Dmean = mean ( Dosis_solo_en_ROI);
Dmin = min ( Dosis_solo_en_ROI);
Dmax = max ( Dosis_solo_en_ROI);
Grosor_bin = (Dmax - Dmin) / 40;
[ N, Edges] = histcounts ( Dosis_solo_en_ROI, 200); % eg 18 bin height and 19 edges. Edges: left bins limit
% histogram (Dosis_solo_en_ROI)     plot diferential histogram
num_barritas = length ( N);
total_elementos = sum ( N);
ancho_bin = Edges ( 2) - Edges ( 1);  % assumption: bins and edges equidistributed
abcisas = [ancho_bin/2 : ancho_bin : Edges( end)];
bines_antes = round ( Edges( 1)/ ancho_bin);
ordenadas = ones ( 1, bines_antes + num_barritas);
ordenadas ( 1 : bines_antes) = 100;
for i = 1 : num_barritas
    total_disminuyente = sum ( N (i:end));
    porciento_disminuyente = 100 * total_disminuyente / total_elementos;
    ordenadas (bines_antes + i) = porciento_disminuyente;
end
histograma_sinmodificar = [abcisas; ordenadas];
%ordenadas(bines_antes+1:bines_antes + num_barritas);
vectordiferencia = abs(ordenadas-98);
minimodiferencia = min(vectordiferencia);
posicion = find(vectordiferencia == minimodiferencia);
dosiscercana = abcisas(posicion);
length(dosiscercana);
if length(dosiscercana) == 1
    D98 = dosiscercana;
elseif length(dosiscercana) == 2
    D98 = mean(dosiscercana);
elseif length(dosiscercana) == 3
    D98 = mean(dosiscercana);
else minimodiferencia == 0
    D98 = dosiscercana;
end
axes ( handles.axes2)
set ( handles.axes2, 'NextPlot', 'Replace')
plot (abcisas, ordenadas)
xlabel ('Dose [Gy]');
ylabel ('Volume [%]');
handles.VolumenROI = VolumenROI;
handles.estudio3D = BaseDosisInterp;
handles.Dmean = Dmean;
handles.D98 = D98;
handles.Dmax = Dmax;
handles.Dmin = Dmin;
handles.FactorEscala = FactorEscala;
handles.anchopixelCT = anchopixelCT;
handles.altopixelCT = altopixelCT;
handles.grosorcorteCT = grosorcorteCT;
handles.OrigenXYZUniversal = [Xoct, Yoct, Zoct];
handles.histograma_sinmodificar = histograma_sinmodificar;
%Este handles.r servira para acumular los depslazamientos de la estructura
%en terminos de fracciones de pixel
handles.r = [0,0,0];
handles.abcisas = abcisas;
handles.ordenadas = ordenadas;
handles.desplaza = [0,0,0];
handles.Zmin = slicelocationmin;
handles.nombre_estructura = nombre_estructura;
handles.numero_ROI = str2num ( numero_ROI_str);
handles.estudio3D_CT = estudio3D;
handles.estudioCTmasROI = estudioCTmasROI;
handles.mascaracontorno3D = mascaracontorno3D;
%Extrae la informacion de las coordenadas iniciales del volumen
coord_i = regionprops3(mascaracontorno3D, "Centroid");
sur = regionprops3(mascaracontorno3D,"SurfaceArea")
centro_i = coord_i.Centroid;
handles.centro_i = centro_i;
handles.dimensiones = [anchopixel, altopixel, grosor_corte];
set (handles.edit57, 'String', sprintf(plan_name))
guidata (hObject, handles)

% --------------------------------------------------------------------
function archivo_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function abrir_estudio_Callback(hObject, eventdata, handles)

function edit1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

% El de ver un corte para atras
estudio3D = handles.estudio3D;
mascaracontorno3D = handles.mascaracontorno3D;
estudioCTmasROI = estudio3D + estudio3D .* mascaracontorno3D;
planoanterior = str2num ( get (handles.edit1, 'String'));
plano = planoanterior - 1;
corte_test = estudioCTmasROI(:, :, plano);
axes ( handles.axes1)
set ( handles.axes1, 'NextPlot', 'Replace')
imagesc (corte_test)
axis off
colormap (gray)
set (handles.edit1, 'String', num2str(plano));

try
    estudio3D2 = handles.estudio3D2_CT;
    corte_test2 = estudio3D2 (:, :, plano);
    axes ( handles.axes2)
    set ( handles.axes2, 'NextPlot', 'Replace')
    imagesc ( corte_test2)
    axis off
    colormap (gray)
catch  % Que todavia no se han efectuado desplazamientos
    lasterr;
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

% el de ver un corte para adelante
estudio3D = handles.estudio3D;
mascaracontorno3D = handles.mascaracontorno3D;
estudioCTmasROI = estudio3D + estudio3D .* mascaracontorno3D;
planoanterior = str2num ( get (handles.edit1, 'String'));
plano = planoanterior + 1;
corte_test = estudioCTmasROI (:, :, plano);
set ( handles.axes1, 'NextPlot', 'Replace')
axes ( handles.axes1)
imagesc (corte_test)
axis off
colormap (gray)
set (handles.edit1, 'String', num2str(plano));

try
    estudio3D2 = handles.estudio3D2_CT;
    corte_test2 = estudio3D2 (:, :, plano);
    axes ( handles.axes2)
    set ( handles.axes2, 'NextPlot', 'Replace')
    imagesc ( corte_test2)
    axis off
    colormap (gray)
catch  % Que todavia no se han efectuado desplazamientos
    lasterr;
end

%Funcion para convertir DVH integral a diferencial
function diffH = int2diff(intH)
    diffH(:,1) = intH(:,1);
    for i=1:size(intH,1) - 1
        diffH(i, 2) = intH(i, 2) - intH(i+1, 2);
    end

%Funcion para obtener la dosis media
function DVHmean = meanDVHdose(diffH)
DVHmean = 0;
Vtotal = sum(diffH(:,2));
for j = 1:size(diffH,1)
    DVHmean = DVHmean + diffH(j,1) * diffH(j,2)/Vtotal;
end

%Funcion para rotar
function [imagen_a_graficar, imagen_rotada_a_graficar, columnaIso, planoIso, filaIso, mascaracontorno3D2, tipogiro, plano_medio] = rota(estudio3D_CT, mascaracontorno3D, anchopixelCT, altopixelCT, grosorcorteCT, origen, grados,i1, i2, i3)
tipogiro = questdlg ('Select the rotation', 'Rotate', 'Roll', 'Pitch', 'Yaw', 'Roll');
switch length (tipogiro)
    case 4
        tipogiro = 'Rol';
    case 5
        tipogiro = 'Pit';
end
x0 = origen ( 1);
y0 = origen ( 2);
z0 = origen ( 3);
[filas, columnas, planos] = size ( mascaracontorno3D);
IsocentroX = i1;
IsocentroY = i2;
IsocentroZ = i3;
columnaIso = round ( abs( IsocentroX - x0)/anchopixelCT);
filaIso = round ( abs( IsocentroY - y0)/altopixelCT);
planoIso = round ( abs( IsocentroZ - z0)/grosorcorteCT);
colsDer = columnas - columnaIso;
columnas_completar = abs ( 2 * columnaIso - 1 - columnas);
filas_completar = abs ( 2 * filaIso - 1 - filas);
planos_completar = abs ( 2 * planoIso - 1 - planos);
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

% Ahora vamos a hacer el giro. En dependencia del tipo de giro, acotar los
% planos necesarios e imrotatearlos, de uno en uno.
grados_num = str2num(grados);
% gira grados_num en sentido antihorario (>0) y en sentido horario (<0)
%Crop  el volumen de salida sea el mismo tamano que el volumen de entrada, recortando el volumen girado para ajustarlo.
% [0 0 1] es el eje de rotacion

% Luego hay que quitar lo que se sumo para que el isocentro quedara en el
% medio
dimensiones2 = size (mascaracontorno3D20); % Recordar que la matriz dosis aumentó de tamaño
switch tipogiro   % Se gira y despues hay que recortar para regresar el estudio 3D al tamaño original
    case 'Rol'  % Se gira como mirando la axial de frente, por el eje Z
        i = 1;
        terminator = 0;
        while terminator == 0
            planoRoiExploratorio = mascaracontorno3D20 (:, :, i);
            if any ( planoRoiExploratorio(:))   % Entra en este if, entrando al volumen a rolear, si hay algun valor distinto de cero
                primer_plano_roleable = i;
                terminator = 1;
            end
            i = i + 1;
        end
        terminator = 0;
        while terminator == 0
            nuevoPlanoRoiExploratorio = mascaracontorno3D20 (:, :, i);  % Y entra en este if, saliendo del volumen a rolear, si ya todos son de nuevo ceros
            if all (nuevoPlanoRoiExploratorio(:) == 0)
                ultimo_plano_roleable = i-1;
                terminator = 1;
            end
            i = i + 1;
        end
        total_planos_rolear = ultimo_plano_roleable - primer_plano_roleable + 1;
        plano_medio = round (0.5 * (primer_plano_roleable + ultimo_plano_roleable));
        for i = primer_plano_roleable : ultimo_plano_roleable
                plano_imrotatear = mascaracontorno3D20 (:, :, i);
                plano_rotado = imrotate (plano_imrotatear, grados_num, 'crop');
                mascaracontorno3D20 (:, :, i) = plano_rotado;
                pause (0.05)
        end
        pause (0.5)
        
            
    case 'Yaw'    % Se gira como mirando el coronal de frente, por el eje Y.
        i = 1;
        terminator = 0;
        while terminator == 0
            planoRoiExploratorio = mascaracontorno3D20 (i, :, :);
            if any ( planoRoiExploratorio(:))   % Entra en este if, si hay algun valor distinto de cero
                primer_plano_yaweable = i;
                terminator = 1;
            end
            i = i + 1;
        end
        terminator = 0;
        while terminator == 0
            nuevoPlanoRoiExploratorio = mascaracontorno3D20 (i, :, :);
            if all (nuevoPlanoRoiExploratorio(:) == 0)
                ultimo_plano_yaweable = i-1;
                terminator = 1;
            end
            i = i + 1;
        end
        total_planos_yawear = ultimo_plano_yaweable - primer_plano_yaweable + 1;
        plano_medio = round (0.5 * (primer_plano_yaweable + ultimo_plano_yaweable));
        for i = primer_plano_yaweable : ultimo_plano_yaweable
            plano_imrotatear = mascaracontorno3D20 (i, :, :);
            plano_imrotatear = permute ( plano_imrotatear, [3 2 1]);
            plano_imrotatear = flip ( plano_imrotatear, 1);
            plano_imrotateado = imrotate ( plano_imrotatear, grados_num, 'crop');
            plano_imrotateado = flip ( plano_imrotateado, 1);
            plano_insertar = permute ( plano_imrotateado, [ 3 2 1]);
            mascaracontorno3D20 (i, :, :) = plano_insertar;
            pause (0.05)
        end
        pause (0.5)
        
    case 'Pit' % Se gira como mirando el sagital de frente,  por el eje X
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
        pause (0.5)

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

j = 1;
t = 0;
while t == 0
    plano_e = mascaracontorno3D20 (:, :, j);
    if any ( plano_e(:))  
        plano_i = j;
        t = 1;
    end
    j = j + 1;
end
t = 0;
while t == 0
    plano_en = mascaracontorno3D20 (:, :, j);  % Y entra en este if, saliendo del volumen a rolear, si ya todos son de nuevo ceros
    if all (plano_en(:) == 0)
        plano_f = j-1;
        t = 1;
    end
    j = j + 1;
end
planos_tot = plano_f - plano_i + 1;
plano_m = round (0.5 * (plano_i + plano_f));

a = size(estudio3D_CT);
mascaracontorno3D20 = imresize3(mascaracontorno3D20, [a(1) a(2) a(3)]);
mascaracontorno3D = imresize3(mascaracontorno3D, [a(1) a(2) a(3)]);

switch tipogiro
    case 'Rol'
        imagen_a_graficar = estudio3D_CT(:, :, plano_m) + estudio3D_CT(:, :, plano_m) .* mascaracontorno3D (:, :, plano_m);
        imagen_rotada_a_graficar = estudio3D_CT(:, :, plano_m) + estudio3D_CT(:, :, plano_m) .* mascaracontorno3D20(:, :, plano_m)+ estudio3D_CT(:, :, plano_m) .* mascaracontorno3D (:, :, plano_m);
    case 'Yaw'
        imagen_a_graficar = estudio3D_CT(plano_m, :, :) + estudio3D_CT(plano_m, :, :) .* mascaracontorno3D (plano_m, :, :);
        imagen_a_graficar = flip ( permute( imagen_a_graficar, [3 2 1]), 1);
        imagen_rotada_a_graficar = estudio3D_CT(plano_m, :, :) + estudio3D_CT(plano_m, :, :) .* mascaracontorno3D20 (plano_m, :, :)+ estudio3D_CT(plano_m, :, :) .* mascaracontorno3D (plano_m, :, :);
        imagen_rotada_a_graficar = flip ( permute( imagen_rotada_a_graficar, [3 2 1]), 1);
    case 'Pit'
        imagen_a_graficar = estudio3D_CT(:, plano_m, :) + estudio3D_CT(:, plano_m, :) .* mascaracontorno3D (:, plano_m, :);
        imagen_a_graficar = permute( imagen_a_graficar, [1 3 2]);
        imagen_rotada_a_graficar = estudio3D_CT(:, plano_m, :) + estudio3D_CT(:, plano_m, :) .* mascaracontorno3D20 (:, plano_m, :) + estudio3D_CT(:, plano_m, :) .* mascaracontorno3D (:, plano_m, :);
        imagen_rotada_a_graficar = permute( imagen_rotada_a_graficar, [1 3 2]);
end

mascaracontorno3D2 = mascaracontorno3D20;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)

% LEE EL HDV DEL PLAN

numero_ROI = handles.numero_ROI; % numero de la ROI en el fichero RS que no corresponde a esa ROI en el fichero RD
nuevo_numeroROI = 1; % Para encontrar el que corresponde
VolumenROI = 0.001 * handles.VolumenROI;  %Porque el otro boton lo calculo en mm^3 y del plan parece que lo lee en cm^3
nombre_estructura = handles.nombre_estructura;
abcisas = handles.abcisas;
ordenadas = handles.ordenadas;
carpeta_madre = handles.cadenaString;
[fichero_dosis, direccion_dosis_y_estructuras] = ruta (carpeta_madre, 2);
[fichero_estructuras, direccion_dosis_y_estructuras] = ruta (carpeta_madre, 1);
cd (direccion_dosis_y_estructuras);
RDdosis = dicomread (fichero_dosis);

RDinfo = dicominfo (fichero_dosis);
RSinfo = dicominfo (fichero_estructuras);
% el numero de ROI que tenemos esta en otra parte
terminator = 0;
i = 1;
prompt={'Cuantas estructuras trajo en total?', 'Y de esas, cuántas no son del paciente?'};
name='Una ayudita aqui por favor';
numlines=1;
defaultanswer = {'41', '4'}
answer=inputdlg(prompt,name,numlines,defaultanswer);
Estructs_tot = str2num ( answer{ 1});
Estructuras_malas = str2num ( answer{2});
Estructuras_buenas = Estructs_tot - Estructuras_malas;
numero_buscado = Estructuras_buenas - numero_ROI;
while terminator == 0
    eval (['numero_mirado = RDinfo.DVHSequence.Item_', num2str( i), '.DVHReferencedROISequence.Item_1.ReferencedROINumber;'])
    if numero_mirado == numero_buscado
        nuevo_numero_ROI = i;    
        terminator = 1;  % Si es la misma ROI, esta es la tipa
    end
  i = i + 1;
end

dosis3D = squeeze ( double ( RDdosis * RDinfo.DoseGridScaling));
dimens = size (dosis3D);
grosor_corte_dosis = RDinfo.GridFrameOffsetVector; 
anchopixel_dosis = RDinfo.PixelSpacing( 2);
altopixel_dosis = RDinfo.PixelSpacing( 1);
origen_dosis = RDinfo.ImagePositionPatient;

NNRs = num2str (nuevo_numero_ROI);
eval (['Tipo_Histograma = RDinfo.DVHSequence.Item_', NNRs,'.DVHType;'])
eval (['DVHDosisunidades = RDinfo.DVHSequence.Item_', NNRs,'.DoseUnits;'])  % Supondremos momentaneamente Gy
eval (['DVHDosisvolumen = RDinfo.DVHSequence.Item_', NNRs,'.DVHVolumeUnits;'])  % Supondremos sean cc
eval (['DVHData = RDinfo.DVHSequence.Item_', NNRs,'.DVHData;']);  % Leemos la info del DVH organizada de la manera rara esa que tiene el fichero, bins, volumen acumulado
eval (['DVHDosisBins = cumsum(DVHData (1:2:end)) * RDinfo.DVHSequence.Item_', NNRs,'.DVHDoseScaling;']) % La abcisa y ya escalada a la unidad de dosis
%extraer volumenes, las ordenadas
DVHVolumenes = DVHData(2:2:end);  % En cc. Total de cc que tienen esa dosis o menos en el HDV integral. Ordenada en el plot, si fuera de absoluto.
VolumenROIe = DVHVolumenes ( 1);
DVHVolumenes_rel = 100 * DVHVolumenes / VolumenROIe; % En por ciento. Ordenada en el plot que es normalizado al volumen de la roi
DVH(:,1) = DVHDosisBins;
DVH(:,2) = DVHVolumenes_rel;

set ( handles.axes2, 'NextPlot', 'Add')
figure
plot (DVH(:,1), DVH(:,2),'.r')
DVH(:,1);
DVH(:,2);
hx = abcisas.';
hy = ordenadas.';
hold on
plot(abcisas, ordenadas,'.b')
xlabel ('Dose [Gy]');
ylabel ('Volume [%]');
Dmeane = meanDVHdose (int2diff(DVH));
VolumenROI = handles.VolumenROI;
Dmean = handles.Dmean;
D95 = handles.D95;
set (handles.edit14, 'String', sprintf('%.2f', Dmean))
set (handles.edit5, 'String', sprintf('%.3f', VolumenROI/1000))
set (handles.edit18, 'String', sprintf('%.2f', Dmeane))
set (handles.edit17, 'String', sprintf('%.2f', VolumenROIe))
set (handles.edit15, 'String', sprintf('%.2f',  D95))
%set (handles.edit16, 'String', sprintf('%.2f',  D2))

function edit2_Callback(hObject, eventdata, handles)

region = handles.region;
set (handles.edit2, 'String', num2str ( region))

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit4_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [nombre_fichero, direccion_fichero] = ruta (carpeta_madre, selector)
% Funcion para leer el nombre y la direccion del fichero
terminator = 0;
i = 1;
carpeta_hija = [ carpeta_madre, '/Dose'];
cd (carpeta_hija)
ficherio = dir ( carpeta_hija);
elementos = length ( ficherio);
switch selector
    case 1 % entregar carpeta y fichero de estructuras
        while terminator == 0 
            try   % tratar de cazar al fichero que se llama RS.dcm
                fichero_est = ficherio(i).name;
                if fichero_est (1:2) == 'RS'  % este es el tipo
                    nombre_fichero = fichero_est;
                    terminator = 1;
                elseif i == elementos
                    terminator = 1;
                end
            catch
            end
            i = i + 1;
        end
    case 2 % entregar carpeta y fichero de dosis
        while terminator == 0 
            try   % tratar de cazar al fichero que se llama RD.dcm
                fichero_est = ficherio(i).name;
                if fichero_est (1:2) == 'RD'  % este es el tipo
                    nombre_fichero = fichero_est;
                    terminator = 1;
                elseif i == elementos
                    terminator = 1;
                end
            catch
            end
            i = i + 1;
        end
    case 3 % entregar carpeta de CT y el primer nombre para que no se quede vacio
        carpeta_hija = [ carpeta_madre, '/CT'];
        cd (carpeta_hija)
        ficherio = dir ( carpeta_hija);
        elementos = length ( ficherio);
        nombre_fichero = ficherio(4).name;
    case 4 % entregar carpeta y fichero del plan
        while terminator == 0 
            try   % tratar de cazar al fichero que se llama RP.dcm
                fichero_est = ficherio(i).name;
                if fichero_est (1:2) == 'RP'  % este es el tipo
                    nombre_fichero = fichero_est;
                    terminator = 1;
                elseif i == elementos
                    terminator = 1;
                end
            catch
            end
            i = i + 1;
        end
end
direccion_fichero = carpeta_hija;
        

function edit7_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)

% Boton para las rotaciones

mascaracontorno3D = handles.mascaracontorno3D;
%auxiliar temporario para acumular traslacion y rotacion
mascaracontorno3D2 = handles.mascaracontorno3D2;
dimensiones = handles.dimensiones;
mascaracontorno3D20 = mascaracontorno3D;
mascaracontorno3D = mascaracontorno3D2;
r = handles.r;
desplaza = handles.desplaza;
anchopixelCT = handles.anchopixelCT;
altopixelCT = handles.altopixelCT;
grosorcorteCT = handles.grosorcorteCT;
dimCT = [anchopixelCT, altopixelCT, grosorcorteCT];
origen = handles.OrigenXYZUniversal;
estudio3D_CT = handles.estudio3D_CT;
isocenter = handles.isocenter;
isocentro = [isocenter(1), isocenter(2), isocenter(3)];
grados = get(handles.edit7, 'String');
[imagen_a_graficar, imagen_rotada_a_graficar, columnaIso, planoIso, filaIso, mascaracontorno3D2, tipogiro, plano_medio] = rota(estudio3D_CT, mascaracontorno3D, anchopixelCT, altopixelCT, grosorcorteCT, origen, grados,isocenter(1), isocenter(2), isocenter(3));
imagendoble = imagen_a_graficar + imagen_rotada_a_graficar;

axis image
axes ( handles.axes1)
set (handles.axes1, 'NextPlot', 'Replace')
imagesc (imagen_a_graficar)
%imagesc (matrizaux(:, :, plano_medio))
axis image
colormap gray
set ( handles.axes1, 'NextPlot', 'Add')
axes ( handles.axes2)
set ( handles.axes2, 'NextPlot', 'Replace')
imagesc ( imagen_rotada_a_graficar)
axis image
colormap gray

set ( handles.axes2, 'NextPlot', 'Add')
        
switch tipogiro
    case 'Rol'
        axes ( handles.axes1)
        plot ( columnaIso, filaIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
        axes ( handles.axes2)
        plot ( columnaIso, filaIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
    case 'Yaw'
        axes ( handles.axes1)
        plot ( columnaIso, planoIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
        axes ( handles.axes2)
        plot (columnaIso, planoIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
    case 'Pit'
        axes ( handles.axes1)
        plot ( planoIso, filaIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
        axes ( handles.axes2)
        plot ( planoIso, filaIso, 'r+', 'MarkerSize', 20)
        axis off
end


estructura_inicial = regionprops3(mascaracontorno3D, "Centroid");
centro = estructura_inicial.Centroid;
estructura_modificada = regionprops3(mascaracontorno3D2, "Centroid");
centro1 =  estructura_modificada.Centroid;
V = centro1 - centro;
r = r + V;
V = V .* dimensiones;
desplaza = desplaza + V;
distancia_1 = sqrt(V*V');
restpixelesx = r(1)-floor(r(1));
restpixelesy = r(2)-floor(r(2));
restpixelesz = r(3)-floor(r(3));
r = [restpixelesx, restpixelesy, restpixelesz];

handles.mascaracontorno3D2 = mascaracontorno3D2;
handles.plano_medio = plano_medio;
handles.columnaIso = columnaIso;
handles.r = r;
handles.desplaza = desplaza;
handles.planoIso = planoIso;
handles.filaIso = filaIso;
guidata (hObject, handles)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)

% BOTON de las Traslaciones

estudio3D_CT = handles.estudio3D_CT;
mascaracontorno3D = handles.mascaracontorno3D;
r = handles.r;
desplaza = handles.desplaza;
mascaracontorno3D2 = mascaracontorno3D;  % Por si no entra en el primer modificador y no se generara 
anchopixel = handles.anchopixelCT;
altopixel = handles.altopixelCT;
grosor_corte = handles.grosorcorteCT;
dimens = size (mascaracontorno3D);
milimeter = get(handles.edit30, 'String');
milimeter_num = str2num(milimeter);

tipodesp = questdlg ('Select the translation', 'Translate', 'x [mm]', 'y [mm]', 'z [mm]', 'x')

switch tipodesp
    case 'x [mm]'
        columnas_empujar = abs ( round ( milimeter_num / anchopixel));
        if milimeter_num > 0 
        % Desplazamos hacia la derecha
            mascaracontorno3D2 = cat (2, zeros(dimens(1), columnas_empujar, dimens(3)), mascaracontorno3D(:, 1:dimens(2)-columnas_empujar, :));
        else 
        % Desplazamos hacia la izquierda
            mascaracontorno3D2 = cat (2, mascaracontorno3D(:, columnas_empujar+1:end, :), zeros(dimens(1), columnas_empujar, dimens(3)));
        end
    case 'y [mm]'
        filas_empujar = abs ( round ( milimeter_num / altopixel));
        if milimeter_num > 0
        % Desplazamos hacia el piso:
            mascaracontorno3D2 = cat (1, zeros(filas_empujar, dimens(2), dimens(3)), mascaracontorno3D2(1:dimens(1)-filas_empujar, :, :));
        else
        % Desplazamos hacia el techo
            mascaracontorno3D2 = cat (1, mascaracontorno3D2(filas_empujar+1:end, :, :), zeros(filas_empujar, dimens(2), dimens(3)));
        end
    case 'z [mm]'
        planos_empujar = abs ( round ( milimeter_num / grosor_corte));
        if milimeter_num > 0
        % Desplazamos hacia la cabeza
            mascaracontorno3D2 = cat (3, zeros(dimens(1), dimens(2), planos_empujar), mascaracontorno3D2(:, :, 1:dimens(3)-planos_empujar));
        else
        % Desplazamos hacia los pies
            mascaracontorno3D2 = cat (3, mascaracontorno3D2(:, :, planos_empujar+1:end), zeros(dimens(1), dimens(2), planos_empujar));
        end
        
end

a = size(estudio3D_CT);
mascaracontorno3D2 = imresize3(mascaracontorno3D2, [a(1) a(2) a(3)]);
mascaracontorno3D = imresize3(mascaracontorno3D, [a(1) a(2) a(3)]);
estudio_dobleimagen3D2 = estudio3D_CT + estudio3D_CT .* mascaracontorno3D2;
estudio_dobleimagen = estudio3D_CT + estudio3D_CT .* mascaracontorno3D;

axes ( handles.axes2)
imagesc ( estudio_dobleimagen3D2 (:, :, 100))
axis off
axes ( handles.axes1)
imagesc ( estudio_dobleimagen (:, :, 100))
axis off
handles.mascaracontorno3D2 = mascaracontorno3D2;
guidata (hObject, handles)



% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)

% Presenta el DVH desplazado
dimensiones = handles.dimensiones;
mascaracontorno3D2 = handles.mascaracontorno3D2; %Trae la estructura trasladada
BaseDosisInterp = handles.estudio3D;
FactorEscala = handles.FactorEscala;
anchopixelCT = handles.anchopixelCT;
altopixelCT = handles.altopixelCT;
grosorcorteCT = handles.grosorcorteCT;
isocenter = handles.isocenter;
centro_i = handles.centro_i;
histograma_sinmodificar = handles.histograma_sinmodificar;
r = handles.r;
desplaza = handles.desplaza;
a = size(mascaracontorno3D2);
BaseDosisInterp = imresize3(BaseDosisInterp, [a(1) a(2) a(3)]);
coord_mod = regionprops3(mascaracontorno3D2, "Centroid");
centro_m =  coord_mod.Centroid;
V = centro_m - centro_i;
isoc = isocenter';
distcompleta = centro_i - isoc;
disttoiso = sqrt(distcompleta*distcompleta')/10;
distancia_f = sqrt(V*V');
var = size(BaseDosisInterp);
Xq = 1:var(1) + r(1);
Yq = 1:var(2) + r(2);
Zq = 1:var(3) + r(3);
[X, Y, Z] = meshgrid(Xq, Yq, Zq);
NuevaDosis = interp3(BaseDosisInterp, X, Y, Z);
Dosis_solo_en_ROI = BaseDosisInterp .* mascaracontorno3D2;
%Dosis_solo_en_ROI = NuevaDosis .* mascaracontorno3D2;
Dosis_solo_en_ROI = FactorEscala * Dosis_solo_en_ROI (Dosis_solo_en_ROI>0);
CantidadV = length ( Dosis_solo_en_ROI);
VolumenROIr = CantidadV * anchopixelCT * altopixelCT * grosorcorteCT;
Dmeanr = mean (Dosis_solo_en_ROI);
Dminr = min ( Dosis_solo_en_ROI);
Dmaxr = max ( Dosis_solo_en_ROI);
Grosor_bin = (Dmaxr - Dminr) / 50;
[N, Edges] = histcounts (Dosis_solo_en_ROI, 500); % Por ejemplo 19 edges y 18 alturas de barritas
num_barritas = length (N); % el numero de barritas
total_elementos = sum (N);
ancho_bin = Edges(2) - Edges (1);  % asumiendo que sean equidistribuidos
abcisas = [ancho_bin/2 : ancho_bin : Edges( end)];
bines_antes = round ( Edges( 1)/ ancho_bin);
ordenadas = ones ( 1, bines_antes + num_barritas);
ordenadas ( 1 : bines_antes) = 100;
for i = 1 : num_barritas
    total_disminuyente = sum ( N (i:end));
    porciento_disminuyente = 100 * total_disminuyente / total_elementos;
    ordenadas (bines_antes + i) = porciento_disminuyente;
end

vectordiferencia = abs(ordenadas-98);
minimodiferencia = min(vectordiferencia);
posicion = find(vectordiferencia == minimodiferencia);
dosiscercana = abcisas(posicion);
length(dosiscercana);
if length(dosiscercana) == 2
    D98r = mean(dosiscercana);
elseif length(dosiscercana) == 3
    D98r = mean(dosiscercana);
elseif minimodiferencia == 0
    D98r = dosiscercana;
else
    D98r = mean(dosiscercana);
%elseif ordenadas(posicion) > 95
%    vector2abcisas = abcisas(posicion:posicion+1);
%    vector2ordenadas = ordenadas (posicion:posicion+1);
%    D95r = interp1(vector2ordenadas, vector2abcisas, [95])
%elseif ordenadas(posicion) < 95
%    vector2abcisas = abcisas(posicion-1:posicion);
%    vector2ordenadas = ordenadas (posicion-1:posicion);
%    D95r = interp1(vector2ordenadas, vector2abcisas, [95]);
end


cocientedvh = trapz(abcisas, ordenadas)/trapz(histograma_sinmodificar(1,:), histograma_sinmodificar(2,:));

%D95r = interp1(ordenadas(bines_antes+3:bines_antes + num_barritas), abcisas(bines_antes+3:bines_antes + num_barritas), [95]);
%D2r = interp1(ordenadas(bines_antes+3:bines_antes + num_barritas), abcisas(bines_antes+3:bines_antes + num_barritas), [2]);

VolumenROI = handles.VolumenROI;
Dmean = handles.Dmean;
D98 = handles.D98;
Dmax = handles.Dmax;
Dmin = handles.Dmin;

set (handles.edit33, 'String', sprintf('%.2f', Dmean))
set (handles.edit5, 'String', sprintf('%.3f', VolumenROI/1000))
set (handles.edit18, 'String', sprintf('%.2f', Dmeanr))
set (handles.edit25, 'String', sprintf('%.2f', Dmax))
set (handles.edit14, 'String', sprintf('%.2f', cocientedvh))
set (handles.edit19, 'String', sprintf('%.2f',  Dmaxr))
set (handles.edit20, 'String', sprintf('%.2f',  distancia_f))
set (handles.edit50, 'String', sprintf('%.1f', disttoiso));

set ( handles.axes2, 'NextPlot', 'Replace')
figure
plot (abcisas, ordenadas, 'r')
xlabel ('Dose [Gy]');
ylabel ('Volume [%]');
hold on
plot (histograma_sinmodificar(1,:), histograma_sinmodificar(2,:))
legend ({'Displaced DVH', 'Original DVH'},'Location','southwest')

guidata (hObject, handles)


function edit14_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit15_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit16_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit17_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit18_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit19_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit20_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)

%prompt={'Dosis Total', 'Fracciones'};
%name='Índice de gradiente';
%numlines=1;
%defaultanswer = {'21', '1'};
%answer=inputdlg(prompt,name,numlines,defaultanswer);
%Fx = str2num ( answer{ 2});
%Dfx = str2num ( answer{ 1})/Fx;
Dfx = 21;
Fx = 1;
mascaracontorno3D2 = handles.mascaracontorno3D2; %Trae la estructura trasladada
mascaracontorno3D = handles.mascaracontorno3D;
BaseDosisInterp = handles.estudio3D;
FactorEscala = handles.FactorEscala;
anchopixelCT = handles.anchopixelCT;
altopixelCT = handles.altopixelCT;
grosorcorteCT = handles.grosorcorteCT;
estudio3D_CT = handles.estudio3D_CT;
origen = handles.OrigenXYZUniversal;
grados = get(handles.edit30, 'String');
f = waitbar(0, 'Conformity and Gradient Index');
pause(0.5)
se = strel('disk', 30); 
mascaracontorno_extension = zeros(size(mascaracontorno3D2)); % mask with 30 pixels expansion
mascaracontorno_extension1 = zeros(size(mascaracontorno3D));
ma = zeros(size(mascaracontorno3D2));
ma1 = zeros(size(mascaracontorno3D));
waitbar(0.1, f, 'Extracting data');
pause(1)
mascaracontorno_extension = imdilate(mascaracontorno3D2,se); 
mascaracontorno_extension1 = imdilate(mascaracontorno3D,se); 
waitbar(0.2, f, 'Extracting data');
pause(1)
BaseDosisInterp = handles.estudio3D;
FactorEscala = handles.FactorEscala;
a = size(mascaracontorno3D);
BaseDosisInterp = imresize3(BaseDosisInterp, [a(1) a(2) a(3)]);
mascaracontorno3D2 = imresize3(mascaracontorno3D2, [a(1) a(2) a(3)]);
mascaracontorno_extension = imresize3(mascaracontorno_extension, [a(1) a(2) a(3)]);
mascaracontorno_extension1 = imresize3(mascaracontorno_extension1, [a(1) a(2) a(3)]);
Dosis_solo_en_ROI = BaseDosisInterp .* mascaracontorno3D2; 
Dosis_solo_en_ROI1 = BaseDosisInterp .* mascaracontorno3D; 
Dosis_auxiliar = BaseDosisInterp .* mascaracontorno_extension;
Dosis_auxiliar1 = BaseDosisInterp .* mascaracontorno_extension1;
Dosis_solo_en_ROI = FactorEscala * Dosis_solo_en_ROI (Dosis_solo_en_ROI>0);
Dosis_solo_en_ROI1 = FactorEscala * Dosis_solo_en_ROI1 (Dosis_solo_en_ROI1>0);
CantidadV = length ( Dosis_solo_en_ROI);
CantidadV1 = length ( Dosis_solo_en_ROI1);

TV_PI = CantidadV * anchopixelCT * altopixelCT * grosorcorteCT;
TV_PI1 = CantidadV1 * anchopixelCT * altopixelCT * grosorcorteCT;

auxiliar = FactorEscala * Dosis_auxiliar (Dosis_auxiliar>0);
auxiliar1 = FactorEscala * Dosis_auxiliar1 (Dosis_auxiliar1>0);
% Vamos a buscar los voxeles que esten en la ROI y tengan D>=Dfx
waitbar(0.4, f, 'Evaluating data');
pause(1)
i=1;
terminator = 0;
while terminator==0
    rawExplora = mascaracontorno3D2(i,:,:);
    if any(rawExplora(:))
        primerraw = i;
        terminator = 1;
    end
    i = i+1;
end
terminator = 0;
while terminator == 0
    newrawExplora = mascaracontorno3D2(i,:,:);
    if all(newrawExplora(:) == 0)
        lastraw = i-1;
        terminator = 1;
    end
    i = i + 1;
end
waitbar(0.5, f, 'Calculating Conformity Index');
pause(1)
j=1;
terminator = 0;
while terminator==0
    columnExplora = mascaracontorno3D2(:,j,:);
    if any(columnExplora(:))
        primercolumn = j;
        terminator = 1;
    end
    j = j+1;
end
terminator = 0;
while terminator == 0
    newcolumnExplora = mascaracontorno3D2(:,j,:);
    if all(newcolumnExplora(:) == 0)
        lastcolumn = j-1;
        terminator = 1;
    end
    j = j + 1;
end
waitbar(0.6, f, 'Calculating Conformity Index');
pause(1)
k=1;
terminator = 0;
while terminator==0
    planoExplora = mascaracontorno3D2(:,:,k);
    if any(planoExplora(:))
        primerplano = k;
        terminator = 1;
    end
    k = k+1;
end
terminator = 0;
while terminator == 0
    newplanoExplora = mascaracontorno3D2(:,:,k);
    if all(newplanoExplora(:) == 0)
        lastplano = k-1;
        terminator = 1;
    end
    k = k + 1;
end
waitbar(0.7, f, 'Calculating Conformity Index');
pause(1)
Cantidad = 0;
DosisExplora = BaseDosisInterp .* mascaracontorno3D2; 
for i=primerraw:lastraw
    for j=primercolumn:lastcolumn
        for k=primerplano:lastplano
            if DosisExplora(i,j,k) >= (Dfx)
                Cantidad = Cantidad + 1;
                ma(i,j,k) = mascaracontorno3D2(i,j,k);
            end
        end
    end
end
au = regionprops3(ma,"SurfaceArea");
S_PI = au.(1)*anchopixelCT*altopixelCT;
av = regionprops3(mascaracontorno3D2,"SurfaceArea");
S_TV = av.(1)*anchopixelCT*altopixelCT;
waitbar(0.8, f, 'Calculating Gradient Index');
pause(1)
i=1;
terminator = 0;
while terminator==0
    rawExplora1 = mascaracontorno3D(i,:,:);
    if any(rawExplora1(:))
        primerraw1 = i;
        terminator = 1;
    end
    i = i+1;
end
terminator = 0;
while terminator == 0
    newrawExplora1 = mascaracontorno3D(i,:,:);
    if all(newrawExplora1(:) == 0)
        lastraw1 = i-1;
        terminator = 1;
    end
    i = i + 1;
end
waitbar(0.85, f, 'Calculating Gradient Index');
pause(1)
j=1;
terminator = 0;
while terminator==0
    columnExplora1 = mascaracontorno3D(:,j,:);
    if any(columnExplora1(:))
        primercolumn1 = j;
        terminator = 1;
    end
    j = j+1;
end
terminator = 0;
while terminator == 0
    newcolumnExplora1 = mascaracontorno3D(:,j,:);
    if all(newcolumnExplora1(:) == 0)
        lastcolumn1 = j-1;
        terminator = 1;
    end
    j = j + 1;
end
waitbar(0.9, f, 'Calculating Gradient Index');
pause(1)
k=1;
terminator = 0;
while terminator==0
    planoExplora1 = mascaracontorno3D(:,:,k);
    if any(planoExplora1(:))
        primerplano1 = k;
        terminator = 1;
    end
    k = k+1;
end
terminator = 0;
while terminator == 0
    newplanoExplora1 = mascaracontorno3D(:,:,k);
    if all(newplanoExplora1(:) == 0)
        lastplano1 = k-1;
        terminator = 1;
    end
    k = k + 1;
end

Cantidad1 = 0;
DosisExplora1 = BaseDosisInterp .* mascaracontorno3D; 
for i=primerraw1:lastraw1
    for j=primercolumn1:lastcolumn1
        for k=primerplano1:lastplano1
            if DosisExplora1(i,j,k) >= (Dfx)
                Cantidad1 = Cantidad1 + 1;
                ma1(i,j,k) = mascaracontorno3D(i,j,k);
            end
        end
    end
end

au1 = regionprops3(ma1,"SurfaceArea");
S_PI1 = au1.(1)*anchopixelCT*altopixelCT;
av1 = regionprops3(mascaracontorno3D,"SurfaceArea");
S_TV1 = av1.(1)*anchopixelCT*altopixelCT;

CantidadV50 = 0;
CantidadV100 = 0;
CantidadV501 = 0;
CantidadV1001 = 0;
waitbar(1, f, 'Finishing');
pause(1)
for i=1:length(auxiliar) % Cuenta cuantos voxeles tienen la dosis del 50%
    if auxiliar(i) >= (Dfx/2)
        CantidadV50 = CantidadV50 + 1; 
    end
    if auxiliar(i) >= Dfx
        CantidadV100 = CantidadV100 + 1;
    end
end

for i=1:length(auxiliar1) % Cuenta cuantos voxeles tienen la dosis del 50%
    if auxiliar1(i) >= (Dfx/2)
        CantidadV501 = CantidadV501 + 1; 
    end
    if auxiliar1(i) >= Dfx
        CantidadV1001 = CantidadV1001 + 1;
    end
end

VolumenROI50 = CantidadV50 * anchopixelCT * altopixelCT * grosorcorteCT;
PI = CantidadV100 * anchopixelCT * altopixelCT * grosorcorteCT;
Volumen = Cantidad * anchopixelCT * altopixelCT * grosorcorteCT;
CI = round((TV_PI*TV_PI)/(PI*Volumen),2);
GI = round(VolumenROI50 / PI, 2);

DCI = (abs(PI-TV_PI)+abs(Volumen-TV_PI))/(0.5*(S_PI+S_TV));
VolumenROI501 = CantidadV501 * anchopixelCT * altopixelCT * grosorcorteCT;
PI_1 = CantidadV1001 * anchopixelCT * altopixelCT * grosorcorteCT;
Volumen1 = Cantidad1 * anchopixelCT * altopixelCT * grosorcorteCT;
CI1 = round((TV_PI1*TV_PI1)/(PI_1*Volumen1),2);
GI1 = round(VolumenROI501 / PI_1, 2);
DCI1 = (abs(PI_1-TV_PI1)+abs(Volumen1-TV_PI1))/(0.5*(S_PI1+S_TV1));
close(f)
set (handles.edit26, 'String', sprintf('%.2f', CI))
set (handles.edit32, 'String', sprintf('%.2f', GI))
set (handles.edit15, 'String', sprintf('%.2f', CI1))
set (handles.edit16, 'String', sprintf('%.2f', GI1))
set (handles.edit54, 'String', sprintf('%.2f', DCI))
set (handles.edit55, 'String', sprintf('%.2f', DCI1))


function edit25_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit26_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)

%BOTON DE ROTAR DE NUEVO
mascaracontorno3D2 = handles.mascaracontorno3D2;
estructura_inicial = regionprops3(mascaracontorno3D2, "Centroid");
centro = estructura_inicial.Centroid;
centro_i = handles.centro_i;
anchopixelCT = handles.anchopixelCT;
dimensiones = handles.dimensiones;
altopixelCT = handles.altopixelCT;
r = handles.r;
desplaza = handles.desplaza;
grosorcorteCT = handles.grosorcorteCT;
estudio3D_CT = handles.estudio3D_CT;
isocenter = handles.isocenter;
isocentro = [isocenter(1), isocenter(2), isocenter(3)];
origen = handles.OrigenXYZUniversal;
grados = get(handles.edit7, 'String');
[imagen_a_graficar, imagen_rotada_a_graficar, columnaIso, planoIso, filaIso, mascaracontorno3D2, tipogiro, plano_medio] = rota(estudio3D_CT, mascaracontorno3D2, anchopixelCT, altopixelCT, grosorcorteCT, origen, grados, isocenter(1), isocenter(2), isocenter(3));
axes ( handles.axes1)
set (handles.axes1, 'NextPlot', 'Replace')
imagesc (imagen_a_graficar)
axis off
colormap gray
set ( handles.axes1, 'NextPlot', 'Add')
axes ( handles.axes2)
set ( handles.axes2, 'NextPlot', 'Replace')
imagesc ( imagen_rotada_a_graficar)
axis image
colormap gray
set ( handles.axes2, 'NextPlot', 'Add')
        
switch tipogiro
    case 'Rol'
        axes ( handles.axes1)
        plot ( columnaIso, filaIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
        axes ( handles.axes2)
        plot ( columnaIso, filaIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
    case 'Yaw'
        axes ( handles.axes1)
        plot ( columnaIso, planoIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
        axes ( handles.axes2)
        plot (columnaIso, planoIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
    case 'Pit'
        axes ( handles.axes1)
        plot ( planoIso, filaIso, 'r+', 'MarkerSize', 20)  % Para poner una marquita en el iso
        axis off
        axes ( handles.axes2)
        plot ( planoIso, filaIso, 'r+', 'MarkerSize', 20)
        axis off
end

handles.mascaracontorno3D2 = mascaracontorno3D2;
handles.plano_medio = plano_medio;
handles.columnaIso = columnaIso;
handles.planoIso = planoIso;
handles.filaIso = filaIso;
handles.desplaza = desplaza;
guidata (hObject, handles)


function edit30_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)

mascaracontorno3D = handles.mascaracontorno3D;
mascaracontorno3D2 = mascaracontorno3D;
r = handles.r;
desplaza = handles.desplaza;
r = [0,0,0];
desplaza = [0,0,0];
handles.mascaracontorno3D2 = mascaracontorno3D2;
handles.desplaza = desplaza;
handles.r = r;
guidata (hObject, handles)

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit31_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit32_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)

%RECONSTRUCCION 3D DE LAS ESTRUCTURAS
mascaracontorno3D = handles.mascaracontorno3D;
mascaracontorno3D2 = handles.mascaracontorno3D2;
estudio3D_CT = handles.estudio3D_CT;
%columnaIso = handles.columnaIso;
%planoIso = handles.planoIso;
%filaIso = handles.filaIso;
%matrizaux = estudio3D_CT; %de aqui salen las imagenes para 3D graficas
set ( handles.axes2, 'NextPlot', 'Replace')
set ( handles.axes2, 'NextPlot', 'Add')
axes ( handles.axes2)
%X = 1:length(matrizaux(1,:,:));
%Y = 1:length(matrizaux(:,1,:));
%Z = 1:length(matrizaux(:,:,1));
%[X1,Y1,Z1] = meshgrid(X,Y,Z);
%slice(X1,Y1,Z1,matrizaux, 250, [], [])
%arreglito = primerplano1:lastplano1;
%contourslice(mascaracontorno3D2,[],[],arreglito,plano_medio);
%view(3);
%axis image;  
%rotate3d on
dimensiones = handles.dimensiones;
s = regionprops3(mascaracontorno3D, "Centroid");
centro = s.Centroid;
s1 = regionprops3(mascaracontorno3D2, "Centroid");
centro1 = s1.Centroid;
V = centro1 - centro;
V = V .* dimensiones;
distancia = sqrt(V*V');
rotate3d on
%axis image
%quiver3(columnaIso, filaIso, planoIso, centro(1), centro(2), centro(3));
%quiver3(columnaIso, filaIso, planoIso, centro1(1), centro1(2), centro1(3));

Ds = smooth3(mascaracontorno3D2);
fv = isosurface(Ds,0);
hiso = patch(fv, 'FaceColor', [0.6350 0.0780 0.1840], 'Edgecolor', 'none');
isonormals(Ds, hiso)
Ds1 = smooth3(mascaracontorno3D);
fv1 = isosurface(Ds1,0);
hiso1 = patch(fv1, 'FaceColor', [0 0.4470 0.7410], 'Edgecolor', 'none');
isonormals(Ds1, hiso1)
axis tight
camlight
camlight(-80,-10)
lighting gouraud
view(3)
daspect([1,1,1])

function edit33_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit34_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit35_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit36_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
%algoritmo genetico de optimizacion
prompt={'Degree [º]', 'Length [mm]','Number of parents (even)','Cycles','Cross-over', 'Mutation rate'};
name='Parameters';
numlines=1;
defaultanswer = {'0','0','6','20','10','0.3'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
grados = str2num ( answer{ 1});
mm = str2num ( answer{2});
n_cromosomas = str2num( answer{3});
n_poblaciones = str2num( answer{4})-1;
pto_cruce = str2num( answer{5});
mut = str2num( answer{6});
mascaracontorno3D = handles.mascaracontorno3D;
dimensiones = handles.dimensiones;
anchopixelCT = handles.anchopixelCT;
altopixelCT = handles.altopixelCT;
grosorcorteCT = handles.grosorcorteCT;
dimCT = [anchopixelCT, altopixelCT, grosorcorteCT];
origen = handles.OrigenXYZUniversal;
estudio3D_CT = handles.estudio3D_CT;
isocenter = handles.isocenter;
centro_i = handles.centro_i;
IsocentroX = isocenter(1);
IsocentroY = isocenter(2);
IsocentroZ = isocenter(3);
isocentro = [IsocentroX, IsocentroY, IsocentroZ];
tam_cromosoma = 16;
desp_opt = zeros(1, str2num(answer{4}));
mse = zeros(1, str2num(answer{4}));
cromosomas = zeros(tam_cromosoma, 1, n_cromosomas);
for i=1:n_cromosomas
    cromosomas(:,:,i) = randi([0 1], tam_cromosoma, 1);
    i = i+1;
end
desplazamientos = zeros(1, n_cromosomas);
f = waitbar(0, 'Optimization');
pause(.5)
waitbar(0.05, f, sprintf('Generating 1 genetic population from %d for optimization', n_poblaciones+1));
pause(1)
tic
[mascaracontorno3D20, desplazamientos] = ga(n_cromosomas, cromosomas, tam_cromosoma, grados, mm, dimensiones, dimCT, origen, estudio3D_CT, mascaracontorno3D, isocentro, centro_i);
[valores_desplazamientos, indice] = sort(desplazamientos(:), 'descend')
indice1 = indice(1);
indice2 = indice(2);
valores_desplazamientos
desp_opt(1,1) = desplazamientos(indice1);
combinacion = cromosomas(:,1,indice1);
orden_de_giro = transpose(combinacion(1:10,:));
cas = orden_de_giro(1)*(2^9) + (orden_de_giro(2))*(2^8) + (orden_de_giro(3))*(2^7) + (orden_de_giro(4))*(2^6) + (orden_de_giro(5))*(2^5) + (orden_de_giro(6))*(2^4) + (orden_de_giro(7))*(2^3) + (orden_de_giro(8))*(2^2) + (orden_de_giro(9))*(2^1) + (orden_de_giro(10))*(2^0);
padre1 = cromosomas(:, 1, indice1);
padre2 = cromosomas(:, 1, indice2);
%Refresh en la poblacion de cromosomas
cromosomas(:, 1, 1) = padre1;
cromosomas(:, 1, 2) = padre2;
cromosomas(:, 1, 3) = transpose([transpose(padre1(1:pto_cruce)), transpose(padre2(pto_cruce+1:tam_cromosoma))]);
cromosomas(:, 1, 4) = transpose([transpose(padre2(1:pto_cruce)), transpose(padre1(pto_cruce+1:tam_cromosoma))]);
%Los cromosomas que pasan a la siguiente generacion sufren una mutacion aleatoria
for l=1:n_cromosomas
    pto_define_mutacion = rand;
    if pto_define_mutacion >= mut
        pto_mutacion = randi([1 tam_cromosoma]);
        if cromosomas(pto_mutacion,1,l) == 1
            cromosomas(pto_mutacion,1,l) = 0;
        else
            cromosomas(pto_mutacion,1,l) = 1;
        end
    end
    l = l+1;
end
toc
%Aqui inicias las nuevas generaciones mejoradas
for m = 1:n_poblaciones
    waitbar(((m+1)/(n_poblaciones+1))-0.1, f, sprintf('Generating %d genetic populations from %d for optimization', m+1, n_poblaciones+1));
    pause(1)
    m
    tic
    [mascaracontorno3D20, desplazamientos] = ga(n_cromosomas, cromosomas, tam_cromosoma, grados, mm, dimensiones, dimCT, origen, estudio3D_CT, mascaracontorno3D, isocentro, centro_i);
    %Me quedo con los dos valores maximos
    [valores_desplazamientos, indice] = sort(desplazamientos(:), 'descend');
    valores_desplazamientos
    indice1 = indice(1);
    indice2 = indice(2);
    desp_opt(1,m+1) = desplazamientos(indice1);
    combinacion = cromosomas(:,1,indice1);
    orden_de_giro = transpose(combinacion(1:10,:));
    cas = orden_de_giro(1)*(2^9) + (orden_de_giro(2))*(2^8) + (orden_de_giro(3))*(2^7) + (orden_de_giro(4))*(2^6) + (orden_de_giro(5))*(2^5) + (orden_de_giro(6))*(2^4) + (orden_de_giro(7))*(2^3) + (orden_de_giro(8))*(2^2) + (orden_de_giro(9))*(2^1) + (orden_de_giro(10))*(2^0);
    padre1 = cromosomas(:, 1, indice1);
    padre2 = cromosomas(:, 1, indice2);
    %Refresh en la poblacion de cromosomas
    cromosomas(:, 1, 1) = padre1;
    cromosomas(:, 1, 2) = padre2;
    cromosomas(:, 1, 3) = transpose([transpose(padre1(1:pto_cruce)), transpose(padre2(pto_cruce+1:tam_cromosoma))]);
    cromosomas(:, 1, 4) = transpose([transpose(padre2(1:pto_cruce)), transpose(padre1(pto_cruce+1:tam_cromosoma))]);
    %Los cromosomas que pasan a la siguiente generacion sufren una mutacion aleatoria
    for l=1:n_cromosomas
        pto_define_mutacion = rand;
        if pto_define_mutacion >= mut
            pto_mutacion = randi([1 tam_cromosoma]);
            if cromosomas(pto_mutacion,1,l) == 1
                cromosomas(pto_mutacion,1,l) = 0;
            else
                cromosomas(pto_mutacion,1,l) = 1;
            end
        end
    l = l+1;
    end
    m = m+1;
    toc
end
if cromosomas(tam_cromosoma-5,:,indice1) == 0
    roll = -grados;
else
    roll = grados;
end
if cromosomas(tam_cromosoma-4,:,indice1) == 0
    yaw = -grados;
else
    yaw = grados;
end
if cromosomas(tam_cromosoma-3,:,indice1) == 0
    pitch = -grados;
else 
    pitch = grados;
end
if cromosomas(tam_cromosoma-2,:,indice1) == 0
    x = -mm;
else
    x = mm;
end
if cromosomas(tam_cromosoma-1,:,indice1) == 0
    y = -mm;
else
    y = mm;
end
if cromosomas(tam_cromosoma,:,indice1) == 0
    z = -mm;
else
    z = mm;
end
[comb] = print_comb(cas);
waitbar(1, f, 'Finishing');
pause(1)
close(f)
mse = (desp_opt-desp_opt(1,1)).^2;
xn = 1:n_poblaciones+1;
set (handles.edit40, 'String', sprintf('%.1f',  roll))
set (handles.edit41, 'String', sprintf('%.1f',  yaw))
set (handles.edit42, 'String', sprintf('%.1f',  pitch))
set (handles.edit37, 'String', sprintf('%.1f',  x))
set (handles.edit38, 'String', sprintf('%.1f',  y))
set (handles.edit39, 'String', sprintf('%.1f',  z))
set (handles.edit43, 'String', sprintf('%s', comb))
set (handles.edit44, 'String', sprintf('%.2f',  desp_opt(1,n_poblaciones+1)))
set ( handles.axes2, 'NextPlot', 'Replace')
figure
plot (xn, flip(mse), 'ro-')
xlabel ('Generations');
ylabel ('MSE');
load splat.mat
sound(y)
hold on

function edit37_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit38_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit39_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit40_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit41_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit42_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit43_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit44_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit45_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)

% BOTON de volver a trasladar

estudio3D_CT = handles.estudio3D_CT;
mascaracontorno3D2 = handles.mascaracontorno3D2; 
anchopixel = handles.anchopixelCT;
altopixel = handles.altopixelCT;
grosor_corte = handles.grosorcorteCT;
dimens = size (mascaracontorno3D2);
milimeter = get(handles.edit30, 'String');
milimeter_num = str2num(milimeter);

tipodesp = questdlg ('Select the translation', 'Translate', 'x [mm]', 'y [mm]', 'z [mm]', 'x')

switch tipodesp
    case 'x [mm]'
        columnas_empujar = abs ( round ( milimeter_num / anchopixel));
        if milimeter_num > 0 
        % Desplazamos hacia la derecha
            mascaracontorno3D2 = cat (2, zeros(dimens(1), columnas_empujar, dimens(3)), mascaracontorno3D2(:, 1:dimens(2)-columnas_empujar, :));
        else 
        % Desplazamos hacia la izquierda
            mascaracontorno3D2 = cat (2, mascaracontorno3D2(:, columnas_empujar+1:end, :), zeros(dimens(1), columnas_empujar, dimens(3)));
        end
    case 'y [mm]'
        filas_empujar = abs ( round ( milimeter_num / altopixel));
        if milimeter_num > 0
        % Desplazamos hacia el piso:
            mascaracontorno3D2 = cat (1, zeros(filas_empujar, dimens(2), dimens(3)), mascaracontorno3D2(1:dimens(1)-filas_empujar, :, :));
        else
        % Desplazamos hacia el techo
            mascaracontorno3D2 = cat (1, mascaracontorno3D2(filas_empujar+1:end, :, :), zeros(filas_empujar, dimens(2), dimens(3)));
        end
    case 'z [mm]'
        planos_empujar = abs ( round ( milimeter_num / grosor_corte));
        if milimeter_num > 0
        % Desplazamos hacia la cabeza
            mascaracontorno3D2 = cat (3, zeros(dimens(1), dimens(2), planos_empujar), mascaracontorno3D2(:, :, 1:dimens(3)-planos_empujar));
        else
        % Desplazamos hacia los pies
            mascaracontorno3D2 = cat (3, mascaracontorno3D2(:, :, planos_empujar+1:end), zeros(dimens(1), dimens(2), planos_empujar));
        end
        
end
estudio_dobleimagen3D2 = estudio3D_CT + estudio3D_CT .* mascaracontorno3D2;
estudio_dobleimagen = estudio3D_CT + estudio3D_CT .* mascaracontorno3D2;


axes ( handles.axes2)
imagesc ( estudio_dobleimagen3D2 (:, :, 100))
axis off
axes ( handles.axes1)
imagesc ( estudio_dobleimagen (:, :, 100))
axis off
handles.mascaracontorno3D2 = mascaracontorno3D2;
guidata (hObject, handles)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% PERCEPTRON
x = readmatrix('dataset.xls');
x_aux = [x(1,:)',x(2,:)',x(3,:)'];
x = [normalize(x(1,:),'range');normalize(x(2,:),'range');normalize(x(3,:),'range')];
%x = [x(1,:);x(2,:);x(3,:)];
t = readmatrix('margin.xls');
t_aux = t(1,:)';
for i = 1:length(t(1,:))
    if t(1,i) == 1
        t_aux(i) = 1;
    elseif t(2,i) == 1
        t_aux(i) = 2;
    else
        t_aux(i) = 3;
    end
end

d = {};
for i = 1:length(t_aux)
    d(i) = {num2str(t_aux(i))}
    i = i+1;
end

[idx,C] = kmeans(x_aux,3,'Distance','cityblock');
valid = t_aux;
acum_yes = 0;
acum_no = 0;

for i = 1:length(t(1,:))
    if idx(i) == 1
        idx(i) = 2;
    elseif idx(i) == 2
        idx(i) = 1;
    else
        idx(i) = 0.5;
    end
end

for i = 1:length(t(1,:))
    if t_aux(i) == 1
        t_aux(i) = 2;
    elseif t_aux(i) == 2
        t_aux(i) = 1;
    else
        t_aux(i) = 0.5;
    end
end

for i = 1:length(t(1,:))
        valid(i) = (idx(i) == t_aux(i));
        if valid(i) == 1
            acum_yes = acum_yes + 1;
        else 
            acum_no = acum_no + 1;
        end
end

%figure
%cm = confusionchart(idx,t_aux,...
%    'ColumnSummary','column-normalized', ...
%    'RowSummary','row-normalized');
%cm.ColumnSummary = 'column-normalized';
%cm.RowSummary = 'row-normalized';
%cm.Title = 'k-means Confusion Matrix';

porcentaje_si = (acum_yes/(acum_yes+acum_no))*100;
porcentaje_no = (acum_no/(acum_yes+acum_no))*100;
X = table(x(1,:),x(2,:),x(3,:));
numlayer = 20;
numlayer2 = 1*numlayer;
numlayer3 = 1*numlayer;
%net = network(1, 2, [1; 1], [1;0], [0 0; 1 0], [0 1]);
%net = network(1, 3, [1; 1; 1], [1;0;0], [0 0 0; 1 0 0; 0 1 0], [0 0 1]);
%net = network(1, 4, [1; 1; 1; 1], [1;0;0;0], [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0], [0 0 0 1]);%
net = network(1, 8, [1; 1; 1; 1; 1; 1; 1; 1], [1;0;0;0;0;0;0;0], [0 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0], [0 0 0 0 0 0 0 1]);%
net.adaptFcn = 'adaptwb';
%Adapt network with weight and bias learning rules: adaptwb
%Sequential order incremental training w/learning functions
net.divideFcn = 'divideint'; %Set the divide function to dividerand (divide training data randomly).
net.divideParam.trainRatio = 0.8;
net.divideParam.valRatio = 0.0;
net.divideParam.testRatio = 0.2;
net.divideMode =  'sample'; %para redes estáticas. 'time' para redes dinamicas
%DerivFcn: This property defines the derivative function to be used to calculate error gradients and Jacobians 
%when the network is trained using a supervised algorithm, such as backpropagation. 
%You can set this property to the name of any derivative function. 
%For a list of functions, type help nnderivative
net.initFcn = 'initlay';
net.performFcn = 'mse';
net.trainFcn = 'trainrp';%'trainrp';%ADAM % set training function to trainlm (Levenberg-Marquardt backpropagation) dX = lr * dperf/dX
net.trainParam.lr = 0.1;
%net.trainParam.mu = 0.1;
%net.trainParam.mc = 0.1;
net.trainParam.epochs = 1500;
net.trainParam.goal = 0.01;
net.trainParam.delta0 = 0.05;
net.trainParam.min_grad = 0.0;
net.trainParam.deltamax = 100;
net.trainParam.max_fail = 100000000;
%net.trainParam.showCommandLine = 1;
net.plotFcns = {'plotperform', 'plottrainstate', 'ploterrhist', 'plotconfusion', 'plotroc', 'plotregression'};

%set Layer1
net.layers{1}.name = 'Capa 1';
net.layers{1}.dimensions = numlayer;
net.layers{1}.initFcn = 'initnw';
net.layers{1}.transferFcn = 'satlins'; %minmaxscaler
%dropoutLayer(0.3);
%set Layer2
net.layers{2}.name = 'Capa 2';
net.layers{2}.dimensions = numlayer2;%mismo num de neuronas
net.layers{2}.initFcn = 'initnw';
net.layers{2}.transferFcn = 'satlins';
%dropoutLayer(0.3);

%set Layer2
net.layers{3}.name = 'Capa 2';
net.layers{3}.dimensions = numlayer2;%mismo num de neuronas
net.layers{3}.initFcn = 'initnw';
net.layers{3}.transferFcn = 'satlins';

%set Layer2
net.layers{4}.name = 'Capa 2';
net.layers{4}.dimensions = numlayer2;%mismo num de neuronas
net.layers{4}.initFcn = 'initnw';
net.layers{4}.transferFcn = 'satlins';

%set Layer2
net.layers{5}.name = 'Capa 2';
net.layers{5}.dimensions = numlayer2;%mismo num de neuronas
net.layers{5}.initFcn = 'initnw';
net.layers{5}.transferFcn = 'satlins';

%set Layer2
net.layers{6}.name = 'Capa 2';
net.layers{6}.dimensions = numlayer2;%mismo num de neuronas
net.layers{6}.initFcn = 'initnw';
net.layers{6}.transferFcn = 'satlins';

%set Layer2
net.layers{7}.name = 'Capa 2';
net.layers{7}.dimensions = numlayer2;%mismo num de neuronas
net.layers{7}.initFcn = 'initnw';
net.layers{7}.transferFcn = 'satlins';

%set Layer3
net.layers{8}.name = 'Capa 3';
net.layers{8}.dimensions = 3;
net.layers{8}.initFcn = 'initnw';
net.layers{8}.transferFcn = 'softmax';%softmax
[net, tr] = train(net,x, t); %training
tr;
y = net(x); %prediction
a = tr.perf;
b = tr.epoch;
c = tr.vperf;
d = tr.tperf;
%view(net);
%plot(tr.epoch,tr.perf,'r-.')
data = [tr.epoch',tr.perf',tr.tperf'];
writematrix(data,'rn.xls')

outputs = net(x);
errors = gsubtract(t,outputs);
performance = perform(net,t,outputs);

trainTargets = t .* tr.trainMask{1};
%valTargets = t  .* tr.valMask{1};
testTargets = t  .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,outputs);
%valPerformance = perform(net,valTargets,outputs)
testPerformance = perform(net,testTargets,outputs);
%figure, plotperform(tr);
%figure, plottrainstate(tr)
%figure, plotregression(t,outputs)
%figure, ploterrhist(errors)

%figure, plotconfusion(t, outputs);

vol = get(handles.edit47, 'String');
dist = get(handles.edit48, 'String');
nummtts = get(handles.edit51, 'String');
vol =str2num(vol);
dist = str2num(dist);
nummtts = str2num(nummtts);
a = [vol; dist; nummtts];
b = net(a);
[m, index] = max(b);
if index == 1
    ptv = 2.0;
elseif index == 2
    ptv = 1.0;
else 
    ptv = 0.5;
end


function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% Estrategia evolutiva
prompt={'Degree [º]', 'Length [mm]','Number of parents (even)','Cycles','Cross-over', 'Mutation rate'};
name='Parameters';
numlines=1;
defaultanswer = {'0','0','6','20','10','0.3'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
grados = str2num ( answer{ 1});
mm = str2num ( answer{2});
n_cromosomas = str2num( answer{3});
n_poblaciones = str2num( answer{4})-1;
pto_cruce = str2num( answer{5});
mut = str2num(answer{6});
mascaracontorno3D = handles.mascaracontorno3D;
dimensiones = handles.dimensiones;
anchopixelCT = handles.anchopixelCT;
altopixelCT = handles.altopixelCT;
grosorcorteCT = handles.grosorcorteCT;
dimCT = [anchopixelCT, altopixelCT, grosorcorteCT];
origen = handles.OrigenXYZUniversal;
estudio3D_CT = handles.estudio3D_CT;
isocenter = handles.isocenter;
centro_i = handles.centro_i;
IsocentroX = isocenter(1);
IsocentroY = isocenter(2);
IsocentroZ = isocenter(3);
isocentro = [IsocentroX, IsocentroY, IsocentroZ];
tam_cromosoma = 16;
desp_opt = zeros(1, str2num(answer{4}));
desp_ciclo = zeros(1, str2num(answer{3}));
mse = zeros(1, str2num(answer{4}));
cromosomas = zeros(tam_cromosoma, 1, n_cromosomas);
offspring = zeros(tam_cromosoma, 1, n_cromosomas);
aux = zeros(tam_cromosoma, 1, 3);
for i=1:n_cromosomas
    cromosomas(:,:,i) = randi([0 1], tam_cromosoma, 1);
    i = i+1;
end
desplazamientos = zeros(1, n_cromosomas);
for l=1:n_poblaciones
    l
    tic
    for k=1:n_cromosomas
        v = randi([1 n_cromosomas]);
        padre1 = cromosomas(:, :, v);
        pto_define_mutacion = rand;
        if pto_define_mutacion >= mut
            pto_mutacion = randi([1 tam_cromosoma]);
            if padre1(pto_mutacion,1) == 1
                padre1(pto_mutacion,1) = 0;
            else
                padre1(pto_mutacion,1) = 1;
            end
        end
        padre2 = cromosomas(:, 1, k);
        aux(:, 1, 1) = padre2;
        aux(:, 1, 2) = transpose([transpose(padre1(1:pto_cruce)), transpose(padre2(pto_cruce+1:tam_cromosoma))]);
        aux(:, 1, 3) = transpose([transpose(padre2(1:pto_cruce)), transpose(padre1(pto_cruce+1:tam_cromosoma))]);
        [mascaracontorno3D20, desplazamientos] = ga(3, aux, tam_cromosoma, grados, mm, dimensiones, dimCT, origen, estudio3D_CT, mascaracontorno3D, isocentro, centro_i);
        %Me quedo con los dos valores maximos
        [valores_desplazamientos, indice] = sort(desplazamientos(:), 'descend');
        indice1 = indice(1);
        desp_ciclo(1,k) = desplazamientos(indice1);
        offspring(:,1,k) = aux(:,1,indice1);
    end
    cromosomas = offspring;
    [value, ind] = sort(desp_ciclo(:),'descend');
    indice2 = ind(1);
    value(indice2)
    offspring(:,1,indice2);
    l = l+1;
    toc
end
combinacion = offspring(:,1,indice2);
despla = value(indice2)
orden_de_giro = transpose(combinacion(1:10,:));
cas = orden_de_giro(1)*(2^9) + (orden_de_giro(2))*(2^8) + (orden_de_giro(3))*(2^7) + (orden_de_giro(4))*(2^6) + (orden_de_giro(5))*(2^5) + (orden_de_giro(6))*(2^4) + (orden_de_giro(7))*(2^3) + (orden_de_giro(8))*(2^2) + (orden_de_giro(9))*(2^1) + (orden_de_giro(10))*(2^0);
if cas>720
    cas = cas-720;
else
    cas = cas;
end
if combinacion(tam_cromosoma-5,:,indice2) == 0
    roll = -grados;
else
    roll = grados;
end
if combinacion(tam_cromosoma-4,:,indice2) == 0
    yaw = -grados;
else
    yaw = grados;
end
if combinacion(tam_cromosoma-3,:,indice2) == 0
    pitch = -grados;
else 
    pitch = grados;
end
if combinacion(tam_cromosoma-2,:,indice2) == 0
    x = -mm;
else
    x = mm;
end
if combinacion(tam_cromosoma-1,:,indice2) == 0
    y = -mm;
else
    y = mm;
end
if combinacion(tam_cromosoma,:,indice2) == 0
    z = -mm;
else
    z = mm;
end
[comb] = print_comb(cas);
set (handles.edit40, 'String', sprintf('%.1f',  roll))
set (handles.edit41, 'String', sprintf('%.1f',  yaw))
set (handles.edit42, 'String', sprintf('%.1f',  pitch))
set (handles.edit37, 'String', sprintf('%.1f',  x))
set (handles.edit38, 'String', sprintf('%.1f',  y))
set (handles.edit39, 'String', sprintf('%.1f',  z))
set (handles.edit43, 'String', sprintf('%s', comb))
set (handles.edit44, 'String', sprintf('%.2f',  despla))
set ( handles.axes2, 'NextPlot', 'Replace')



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% funcion para alimentar la red neuronal con los parametros para cada mtts
% de cada paciente
dato = [;;];
code_margin = [;;];
prompt={'Number of new datasets'};
name='Data';
numlines=1;
defaultanswer = {'1'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
val = str2num ( answer{ 1});

if val == 1
    prompt={'Total number of mtts', 'Volume [cc]','Distance-to-iso [mm]','Maximum displacement [mm]'};
    name='Dataset 1';
    numlines=1;
    defaultanswer = {' ',' ',' ',' ',' ',' '};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    number = str2num ( answer{ 1});
    volume = str2num ( answer{2});
    dist2iso = str2num( answer{3});
    displacement = str2num( answer{4});
    if displacement > 1.73
        margin = 2;
    elseif displacement >= 1
        margin = 1;
    else
        margin = 0.5;
    end
    radius = ((3*volume)/(4*pi()))^(1/3);
    new_radius = radius + (margin/10);
    new_volume = (4*pi()*new_radius^3)/3;
    first_ratio = (new_volume-volume)/volume;
    if first_ratio > 2 && margin == 2
        margin = 1;
    elseif first_ratio > 2 && margin == 1
        margin = 0.5;
    end
    if margin == 2
        aux_margin = [1;0;0];
        bin_margin = cat(2,code_margin,aux_margin);
    elseif margin == 1
        aux_margin = [0;1;0];
        bin_margin = cat(2,code_margin,aux_margin);
    else 
        aux_margin = [0;0;1];
        bin_margin = cat(2,code_margin,aux_margin);
    end
    value = [number;volume;dist2iso];
    data = cat(2,dato,value);
    handles.data = data;
    handles.bin_margin = bin_margin;

elseif val>1
    for i=1:val
        aa = num2str(i);
        prompt={'Total number of mtts', 'Volume [cc]','Distance-to-iso [mm]','Maximum displacement [mm]'};
        name=aa;
        numlines=1;
        defaultanswer = {' ',' ',' ',' ',' ',' '};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        number = str2num ( answer{ 1});
        volume = str2num ( answer{2});
        dist2iso = str2num( answer{3});
        displacement = str2num( answer{4});
        if displacement > 1.73
            margin = 2;
        elseif displacement >= 1
            margin = 1;
        else
            margin = 0.5;
        end
        radius = ((3*volume)/(4*pi()))^(1/3);
        new_radius = radius + (margin/10);
        new_volume = (4*pi()*new_radius^3)/3;
        first_ratio = (new_volume-volume)/volume;
        if first_ratio > 2 && margin == 2
            margin = 1;
        elseif first_ratio > 2 && margin == 1
            margin = 0.5;
        end
        if margin == 2
            aux_margin = [1;0;0];
            bin_margin = cat(2,code_margin,aux_margin);
        elseif margin == 1
            aux_margin = [0;1;0];
            bin_margin = cat(2,code_margin,aux_margin);
        else 
            aux_margin = [0;0;1];
            bin_margin = cat(2,code_margin,aux_margin);
        end
        value = [number;volume;dist2iso];
        data = cat(2,dato,value);
        dato = data;
        code_margin = bin_margin;
        i = i+1;
        handles.data = data;
        handles.bin_margin = bin_margin;
    end
end
writematrix(data,'dataset.xls')
writematrix(bin_margin,'margin.xls')
guidata (hObject, handles)



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit55_Callback(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit55 as text
%        str2double(get(hObject,'String')) returns contents of edit55 as a double


% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
