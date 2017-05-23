function varargout = TiltedFieldLattice_v2(varargin)
% TILTEDFIELDLATTICE_V2 MATLAB code for TiltedFieldLattice_v2.fig
%      TILTEDFIELDLATTICE_V2, by itself, creates a new TILTEDFIELDLATTICE_V2 or raises the existing
%      singleton*.
%
%      H = TILTEDFIELDLATTICE_V2 returns the handle to a new TILTEDFIELDLATTICE_V2 or the handle to
%      the existing singleton*.
%
%      TILTEDFIELDLATTICE_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TILTEDFIELDLATTICE_V2.M with the given input arguments.
%
%      TILTEDFIELDLATTICE_V2('Property','Value',...) creates a new TILTEDFIELDLATTICE_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TiltedFieldLattice_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TiltedFieldLattice_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TiltedFieldLattice_v2

% Last Modified by GUIDE v2.5 23-May-2017 13:34:01

% Begin initialization code - DO NOT EDIT 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TiltedFieldLattice_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @TiltedFieldLattice_v2_OutputFcn, ...
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

% --- Executes just before TiltedFieldLattice_v2 is made visible.
function TiltedFieldLattice_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TiltedFieldLattice_v2 (see VARARGIN)

% Choose default command line output for TiltedFieldLattice_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TiltedFieldLattice_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = TiltedFieldLattice_v2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%-------------------------------------------------------------------------
%----------------------Magnetic Field-------------------------------------
%-------------------------------------------------------------------------

function MagneticField_Callback(hObject, eventdata, handles)
% hObject    handle to MagneticField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MagneticField as text
%        str2double(get(hObject,'String')) returns contents of MagneticField as a double

% --- Executes during object creation, after setting all properties.
function MagneticField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MagneticField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------------------------------------------------
%---------------Intervortex Distance (nm)----------------------------------
%--------------------------------------------------------------------------
%
function IntervortexDistance_Callback(hObject, eventdata, handles)
% hObject    handle to IntervortexDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IntervortexDistance as text
%        str2double(get(hObject,'String')) returns contents of IntervortexDistance as a double

% --- Executes during object creation, after setting all properties.
function IntervortexDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntervortexDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------------------------------------
%--------------------SIZE WINDOW Vortex Frame (nm)------------------------
%-------------------------------------------------------------------------

function SizeWindowVortexFrame_Callback(hObject, eventdata, handles)
% hObject    handle to SizeWindowVortexFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeWindowVortexFrame as text
%        str2double(get(hObject,'String')) returns contents of SizeWindowVortexFrame as a double

% --- Executes during object creation, after setting all properties.
function SizeWindowVortexFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeWindowVortexFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%---------------------------SIZE WINDOW SURFACE FRAME (nm)-----------------
%--------------------------------------------------------------------------

function SizeWindowSurfaceFrame_Callback(hObject, eventdata, handles)
% hObject    handle to SizeWindowSurfaceFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeWindowSurfaceFrame as text
%        str2double(get(hObject,'String')) returns contents of SizeWindowSurfaceFrame as a double

% --- Executes during object creation, after setting all properties.
function SizeWindowSurfaceFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeWindowSurfaceFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------------------------------------
%--------------------AZIMUTHAL ANGLE (deg)--------------------------------
%-------------------------------------------------------------------------

function Azimuthal_Callback(hObject, eventdata, handles)
% hObject    handle to Azimuthal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Azimuthal as text
%        str2double(get(hObject,'String')) returns contents of Azimuthal as a double

% --- Executes during object creation, after setting all properties.
function Azimuthal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Azimuthal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------------------------------------
%---------------------POLAR ANGLE (deg)-----------------------------------
%-------------------------------------------------------------------------

function PolarAngle_Callback(hObject, eventdata, handles)
% hObject    handle to PolarAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PolarAngle as text
%        str2double(get(hObject,'String')) returns contents of PolarAngle as a double

% --- Executes during object creation, after setting all properties.
function PolarAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PolarAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%-----------SizeWindowSurfaceFrame Lattice Angles: Phi_1, Phi_2 and Phi_3------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------------Phi_1 (deg)--------------------------------------
%--------------------------------------------------------------------------
function Phi1_Callback(hObject, eventdata, handles)
% hObject    handle to Phi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Phi1 as text
%        str2double(get(hObject,'String')) returns contents of Phi1 as a double

% --- Executes during object creation, after setting all properties.
function Phi1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%-------------------------Alpha (deg)--------------------------------------
%--------------------------------------------------------------------------

function Alpha_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Alpha as text
%        str2double(get(hObject,'String')) returns contents of Alpha as a double

% --- Executes during object creation, after setting all properties.
function Alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%-------------CALCULATE ALPHA S -------------------------------------------
%--------------------------------------------------------------------------

function AlphaS_Callback(hObject, eventdata, handles)
% hObject    handle to AlphaS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AlphaS as text
%        str2double(get(hObject,'String')) returns contents of AlphaS as a double

% --- Executes during object creation, after setting all properties.
function AlphaS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%------------Number Neighbors Vortex Frame---------------------------------
%--------------------------------------------------------------------------
function NumberNeighborsVF_Callback(hObject, eventdata, handles)
% hObject    handle to NumberNeighborsVF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumberNeighborsVF as text
%        str2double(get(hObject,'String')) returns contents of NumberNeighborsVF as a double

% --- Executes during object creation, after setting all properties.
function NumberNeighborsVF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberNeighborsVF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%------------Number Neighbors Surface Frame -------------------------------
%--------------------------------------------------------------------------
function NumberNeighborsSF_Callback(hObject, eventdata, handles)
% hObject    handle to NumberNeighborsSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumberNeighborsSF as text
%        str2double(get(hObject,'String')) returns contents of NumberNeighborsSF as a double

% --- Executes during object creation, after setting all properties.
function NumberNeighborsSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberNeighborsSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------------------------
%-------------------------H_c2^c (T)---------------------------------------
%--------------------------------------------------------------------------

function Hc2c_Callback(hObject, eventdata, handles)
% hObject    handle to Hc2c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hc2c as text
%        str2double(get(hObject,'String')) returns contents of Hc2c as a double

% --- Executes during object creation, after setting all properties.
function Hc2c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hc2c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%-------------------------H_c2^ab (T)---------------------------------------
%--------------------------------------------------------------------------

function Hc2ab_Callback(hObject, eventdata, handles)
% hObject    handle to Hc2ab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hc2ab as text
%        str2double(get(hObject,'String')) returns contents of Hc2ab as a double

% --- Executes during object creation, after setting all properties.
function Hc2ab_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hc2ab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------ANISOTROPY VALUE-----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AnisotropyValue_Callback(hObject, eventdata, handles)
% hObject    handle to AnisotropyValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AnisotropyValue as text
%        str2double(get(hObject,'String')) returns contents of AnisotropyValue as a double

% --- Executes during object creation, after setting all properties.
function AnisotropyValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnisotropyValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------------------------------------
%----------------------PUSH BUTTON 1--------------------------------------
%-------------------------------------------------------------------------

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----------------------------------------------------------
%Magnetic Field Features
%----------------------------------------------------------

H=get(handles.MagneticField,'String');     % Magnetic Field (T)

Theta=get(handles.PolarAngle,'String');    % Polar Angle (deg)
ThetaValue=str2num(Theta);

Azimuthal=get(handles.Azimuthal,'String'); % Azimuthal Angle (deg)
AzimuthalValue=str2num(Azimuthal); 

%--------------------------------------------------------
%Windows and Neighbors
%--------------------------------------------------------
% Size Window Vortex Frame (nm)
SWVF=get(handles.SizeWindowVortexFrame,'String'); 
SVF=str2num(SWVF);
% Number of Neighbors Vortex Frame
NeighborsVF=get(handles.NumberNeighborsVF,'String');
NNVF=str2num(NeighborsVF); 

% Size Window Surface Frame (nm)
SWSF=get(handles.SizeWindowSurfaceFrame,'String');  
SSF=str2num(SWSF); 
%Number of Neighbors Surface Frame
NeighborsSF=get(handles.NumberNeighborsSF,'String');
NNSF=str2num(NeighborsSF);
%-------------------------------------------------------

%-------------------------------------------------------
%Material Parameters
%-------------------------------------------------------

Hc2c=get(handles.Hc2c,'String');  % Hc2c(T): Critical Magnetic Field parallel to c-axis
Hperp=str2num(Hc2c);   

Hc2ab=get(handles.Hc2ab,'String'); % Hc2ab(T): Critical Magnetic Field parallel to ab-plane
Hparal=str2num(Hc2ab); 
%-------------------------------------------------------

%-------------------------------------------------------
%Alpha Angle, measured from the Azimuthal angle
%-------------------------------------------------------
Phi1=get(handles.Phi1,'String');
%-------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--Command to obtain the Intervortex Distance--%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d0=48/sqrt(str2num(H));
IntervortexDistanceString=num2str(d0);
set(handles.IntervortexDistance,'String',IntervortexDistanceString);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----Angular positions of all the six vortices of the unit cell ON THE VORTEX FRAME-------% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phi=zeros(1,6);
    Phi(1,1)=str2num(Phi1)+AzimuthalValue;
    Phi(1,2)=str2num(Phi1)+AzimuthalValue+60;
    Phi(1,3)=str2num(Phi1)+AzimuthalValue+120;
    Phi(1,4)=str2num(Phi1)+AzimuthalValue+180;
    Phi(1,5)=str2num(Phi1)+AzimuthalValue+240;
    Phi(1,6)=str2num(Phi1)+AzimuthalValue+300;

for i=1:6
    if(Phi(1,i)<0)
        Phi(1,i)=Phi(1,i)+360;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------Phi corregido---------------------%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phi_corregido=zeros(1,6);% Alpha vector initial correlated to each Phi(1,i)
    for i=1:6
        Phi_corregido(1,i)=Phi(1,i);% Alpha correlated to each Phi(1,i)
    end

Phi_corregido_def=min(Phi_corregido)*(pi/180);%Alpha correlated definitive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----Anisotropy Value--------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Anisotropy=(Hperp/Hparal)^0.5; %Anisotropy ratio
AnisotropyValuestring=num2str(Anisotropy);
set(handles.AnisotropyValue,'String',AnisotropyValuestring);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vrp=[];

%--------------------------------------------------------------------------
%---Anisotropies on each axis From PRB 38, 2439 (1988)---------------------
%--------------------------------------------------------------------------
gammax=sqrt(sqrt(Anisotropy));% Anisotropy on "x-axis"
gammay=1/gammax;% Anisotropy on "y-axis"
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------Basis Vectors on the Vortex Frame---------------------------------
%--------------------------------------------------------------------------
u1y=d0*gammax*cos(Phi_corregido_def);
u1x=d0*gammay*sin(Phi_corregido_def);
u2y=d0*gammax*cos(Phi_corregido_def+pi/3);
u2x=d0*gammay*sin(Phi_corregido_def+pi/3);
%--------------------------------------------------------------------------
nmax= 6*round(SVF/d0);
index=1;
rp1=[];
rp2=[];

for i = -NNVF: NNVF
    for j = -NNVF: NNVF
        if(sqrt((i+j)^2)<=NNVF),
            rp1(index,1)= (i*u1y+j*u2y);
            rp1(index,2)= (i*u1x+j*u2x);
            index=index+1;
        end  
    end
end  

nn=length(rp1);
rp2=[];

for i=1:nn
    rp2(i,1)=rp1(i,1);
    rp2(i,2)=rp1(i,2);
end

rp3=[];
for i=1:nn
    rp3(i,1)=rp2(i,1);
    rp3(i,2)=rp2(i,2);
end

centros2=[];
index=1;

for i=1:nn
    centros2(index,:)=round(rp3(i,:));
    index=index+1;
end

% %[u,v] = EllipseKogan(u1x,u1y,u2x,u2y);%Enclose the unit vortex lattice
%          a=sqrt(((u1x*u2y)^2-(u2x*u1y)^2)/(u2y^2-u1y^2)); % horizontal radius
%         b=sqrt(((u1x*u2y)^2-(u2x*u1y)^2)/(u1x^2-u2x^2)); % vertical radius
%          t=-pi:0.01:pi;
%      	u=a*cos(t);
%          v=b*sin(t);      

%--------------------------------------------------------------------------
%-----VORTEX FRAME PLOT PARAMETERS-----------------------------------------
%--------------------------------------------------------------------------

%plot(handles.graf1,v,u,'g',centros2(:,1),centros2(:,2),'g.','MarkerSize',20)

handles.g1=plot(handles.graf1,centros2(:,1), centros2(:,2),'g.','MarkerSize', 20)

set(handles.graf1,'XLim',[-SVF SVF]);
set(handles.graf1,'YLim',[-SVF SVF]);
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------Alpha Definition---------------------%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Alpha_corr=zeros(1,6);% Alpha vector initial correlated to each Phi(1,i)
    AzimuthalValuePositive=AzimuthalValue;%Converts Azimuthal Angle to a positive value
    if  (AzimuthalValuePositive<0)
        AzimuthalValuePositive=AzimuthalValuePositive+360;
    end    

    for i=1:6
        Alpha_corr(1,i)=abs(AzimuthalValue-Phi(1,i));% Alpha correlated to each Phi(1,i)
        if  (Alpha_corr(1,i)<0)
            Alpha_corr(1,i)=Alpha_corr(1,i)+360;
        end  
    end

Alpha_corr_def=min(Alpha_corr);%Alpha correlated definitive
% Alphastring=num2str(Alpha_corr_def);
% set(handles.Alpha,'String',Alphastring);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----To obtain Anisotropy Value----------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Anisotropy=(Hperp/Hparal)^0.5; %Anisotropy ratio
AnisotropyValuestring=num2str(Anisotropy);
set(handles.AnisotropyValue,'String',AnisotropyValuestring);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------Alpha-S Definition---------------------%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
AlphaSfinal=atand(tand(Alpha_corr_def)*cosd(ThetaValue)/(sqrt(Anisotropy)));%AlphaS definition Paper
AlphaSString=num2str(AlphaSfinal);
set(handles.AlphaS,'String',AlphaSString);
%--------------------------------------------------------------------------

phi_az=Alpha_corr_def*(pi/180);
phi_rot=AzimuthalValue*(pi/180);
vrp=[];

%--------------------------------------------------------------------------
%---Anisotropies on each axis From PRB 38, 2439 (1988)---------------------
%--------------------------------------------------------------------------
gammax=sqrt(sqrt(Anisotropy));% Anisotropy on "x-axis"
gammay=1/gammax;% Anisotropy on "y-axis"
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------Phi corregido---------------------%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1y=d0*gammax*cos(phi_az);
v1x=d0*gammay*sin(phi_az);
v2y=d0*gammax*cos(phi_az+pi/3);
v2x=d0*gammay*sin(phi_az+pi/3);

nmax= 6*round(SSF/d0);
index=1;
rp1=[];
rp2=[];

for i = -NNSF: NNSF
    for j = -NNSF: NNSF
        if(sqrt((i+j)^2)<=NNSF),
            rp1(index,1)= (i*v1y+j*v2y);
            rp1(index,2)= (i*v1x+j*v2x);
            index=index+1;
        end  
    end
end  

nn=length(rp1);

rp2=[];

for i=1:nn
    rp2(i,1)=rp1(i,1)/cosd(ThetaValue);
    rp2(i,2)=rp1(i,2);
end

rp3=[];
for i=1:nn
    rp3(i,1)=rp2(i,1)*cos(phi_rot)-rp2(i,2)*sin(phi_rot);
    rp3(i,2)=rp2(i,1)*sin(phi_rot)+rp2(i,2)*cos(phi_rot);
end

centros2=[];
index=1;

for i=1:nn
    centros2(index,:)=round(rp3(i,:));
    index=index+1;
end
%--------------------------------------------------------------------------
%-----SURFACE FRAME PLOT PARAMETERS----------------------------------------
%--------------------------------------------------------------------------
handles.g2=plot(handles.graf2,centros2(:,1), centros2(:,2),'b.','MarkerSize', 20)
set(handles.graf2,'XLim',[-SSF SSF]);
set(handles.graf2,'YLim',[-SSF SSF]);

