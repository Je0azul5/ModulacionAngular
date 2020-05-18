function varargout = ModulacionAngular(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModulacionAngular_OpeningFcn, ...
                   'gui_OutputFcn',  @ModulacionAngular_OutputFcn, ...
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


% --- Executes just before ModulacionAngular is made visible.
function ModulacionAngular_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;

handles.Amp1 = [];
handles.Amp2 = [];

% Update handles structure
guidata(hObject, handles);


function varargout = ModulacionAngular_OutputFcn(~, eventdata, handles) 
varargout{1} = handles.output;



function Vm_Callback(~, ~, ~)
% hObject    handle to Vm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Vm_CreateFcn(hObject, ~, ~)
% hObject    handle to Vm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vc_Callback(~, eventdata, handles)
% hObject    handle to Vc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function Vc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fm_Callback(hObject, eventdata, handles)
% hObject    handle to Fm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fm as text
%        str2double(get(hObject,'String')) returns contents of Fm as a double


% --- Executes during object creation, after setting all properties.
function Fm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrecPort_Callback(hObject, eventdata, handles)
% hObject    handle to FrecPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function FrecPort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrecPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fc_Callback(hObject, eventdata, handles)
% hObject    handle to Fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fc as text
%        str2double(get(hObject,'String')) returns contents of Fc as a double


% --- Executes during object creation, after setting all properties.
function Fc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function K1_Callback(hObject, eventdata, handles)
% hObject    handle to K1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1 as text
%        str2double(get(hObject,'String')) returns contents of K1 as a double


% --- Executes during object creation, after setting all properties.
function K1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Desv1_Callback(hObject, eventdata, handles)
% hObject    handle to Desv1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Desv1 as text
%        str2double(get(hObject,'String')) returns contents of Desv1 as a double


% --- Executes during object creation, after setting all properties.
function Desv1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Desv1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m1_Callback(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m1P_Callback(hObject, eventdata, handles)
% hObject    handle to m1P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m1P as text
%        str2double(get(hObject,'String')) returns contents of m1P as a double


% --- Executes during object creation, after setting all properties.
function m1P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m1P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bwReal1_Callback(hObject, eventdata, handles)
% hObject    handle to bwReal1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bwReal1 as text
%        str2double(get(hObject,'String')) returns contents of bwReal1 as a double


% --- Executes during object creation, after setting all properties.
function bwReal1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bwReal1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bwMin1_Callback(hObject, eventdata, handles)
% hObject    handle to bwMin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bwMin1 as text
%        str2double(get(hObject,'String')) returns contents of bwMin1 as a double


% --- Executes during object creation, after setting all properties.
function bwMin1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bwMin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DR_Callback(hObject, eventdata, handles)
% hObject    handle to DR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DR as text
%        str2double(get(hObject,'String')) returns contents of DR as a double


% --- Executes during object creation, after setting all properties.
function DR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pprom1_Callback(hObject, eventdata, handles)
% hObject    handle to Pprom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pprom1 as text
%        str2double(get(hObject,'String')) returns contents of Pprom1 as a double


% --- Executes during object creation, after setting all properties.
function Pprom1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pprom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ptotal1_Callback(hObject, eventdata, handles)
% hObject    handle to Ptotal1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ptotal1 as text
%        str2double(get(hObject,'String')) returns contents of Ptotal1 as a double


% --- Executes during object creation, after setting all properties.
function Ptotal1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ptotal1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcularAll.
function calcularAll_Callback(hObject, eventdata, handles)
% hObject    handle to calcularAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Vm = str2double(get(handles.Vm, 'String'));
Vc = str2double(get(handles.Vc, 'String'));
Fm = str2double(get(handles.Fm, 'String'));
K1 = str2double(get(handles.K1,'String'));
K = str2double(get(handles.K,'String'));

if isnan(Vm)  
    f = warndlg('Llenar campos Amplitud Moduladora');
    return;
end
if isnan(Vc)  
    f = warndlg('Llenar campos Amplitud Portadora');
    return;
end
if isnan(Fm)  
    f = warndlg('Llenar campos Frecuencia moduladora');
    return;
end
if isnan(K)  
    f = warndlg('Llenar campos Sensibilidad Fase');
    return;
end
if isnan(get(handles.Fc, 'String'))  
    f = warndlg('Llenar campos Frecuencia portadora');
    return;
end
if isnan(K1)  
    f = warndlg('Llenar campos Sensibilidad Frecuencia');
    return;
end




%Modulación FM
m = (Vm*str2double(get(handles.K1, 'String')))/Fm;
desvFr = m*Fm;
set(handles.Desv1, 'String', desvFr);
set(handles.m1, 'String', m);
Bessel1 = Bessel(m);
set(handles.m1P,'String', (desvFr/str2double(get(handles.FrecMax1, 'String')))*100);
set(handles.bwReal1,'String', 2*((length(Bessel1)-1)*Fm));
set(handles.bwMin1,'String', 2*(desvFr+Fm));
set(handles.DR,'String', desvFr/Fm);
Ampl = [];
for i=1:length(Bessel1)
    Ampl(i) = Bessel1(i)*Vc;
end
handles.Amp1 = Ampl;
Pot = [];
res1 = str2double(get(handles.Res1, 'String'));
Pot(1) = (Ampl(1)^2)/(2*res1);
for i=2:length(Ampl)
    Pot(i) = (Ampl(i)^2)/res1;
end
data = [Bessel1; Ampl; Pot];
set(handles.Table1, 'Data', data);
set(handles.Pprom1,'String', (Vc^2)/(2*res1)); %R=50 Ohms
set(handles.Ptotal1,'String', sum(Pot));

%Modulación PM
m2 = str2double(get(handles.K, 'String'))*Vm;
set(handles.Desv2, 'String', m2);
set(handles.m2, 'String', m2);
Bessel2 = Bessel(m2);
%set(handles.m1P,'String', (desvFr/handles.FrecMax1)*100);
set(handles.bwReal2,'String', 2*((length(Bessel2)-1)*Fm));
set(handles.bwMin2,'String', 2*(m+Fm));
Ampl2 = [];
for i=1:length(Bessel2)
    Ampl2(i) = Bessel2(i)*Vc;
end
handles.Amp2 = Ampl2;
Pot2 = [];
res2 = str2double(get(handles.Res2, 'String'));
Pot2(1) = (Ampl2(1)^2)/(2*res2);
for i=2:length(Ampl2)
    Pot2(i) = (Ampl2(i)^2)/res2;
end
data2 = [Bessel2; Ampl2; Pot2];
set(handles.Table2, 'Data', data2);
set(handles.Pprom2,'String', (Vc^2)/(2*res2));
set(handles.Ptotal2,'String', sum(Pot2));

guidata(hObject, handles);

function B1 = Bessel(x)
J1 = zeros(15,1);
for i = 0:14
    J1(i+1,:) = besselj(i,x);
end
B = J1.';
B1 = [];
for j = 1:15
    if abs(B(j)) > 0.01
        B1(j) = abs(B(j));
    end
end


function K_Callback(hObject, eventdata, handles)
% hObject    handle to K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K as text
%        str2double(get(hObject,'String')) returns contents of K as a double


% --- Executes during object creation, after setting all properties.
function K_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Desv2_Callback(hObject, eventdata, handles)
% hObject    handle to Desv2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Desv2 as text
%        str2double(get(hObject,'String')) returns contents of Desv2 as a double


% --- Executes during object creation, after setting all properties.
function Desv2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Desv2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m2_Callback(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m2 as text
%        str2double(get(hObject,'String')) returns contents of m2 as a double


% --- Executes during object creation, after setting all properties.
function m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m2P_Callback(hObject, eventdata, handles)
% hObject    handle to m2P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m2P as text
%        str2double(get(hObject,'String')) returns contents of m2P as a double


% --- Executes during object creation, after setting all properties.
function m2P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m2P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bwReal2_Callback(hObject, eventdata, handles)
% hObject    handle to bwReal2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bwReal2 as text
%        str2double(get(hObject,'String')) returns contents of bwReal2 as a double


% --- Executes during object creation, after setting all properties.
function bwReal2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bwReal2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bwMin2_Callback(hObject, eventdata, handles)
% hObject    handle to bwMin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bwMin2 as text
%        str2double(get(hObject,'String')) returns contents of bwMin2 as a double


% --- Executes during object creation, after setting all properties.
function bwMin2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bwMin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pprom2_Callback(hObject, eventdata, handles)
% hObject    handle to Pprom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pprom2 as text
%        str2double(get(hObject,'String')) returns contents of Pprom2 as a double


% --- Executes during object creation, after setting all properties.
function Pprom2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pprom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ptotal2_Callback(hObject, eventdata, handles)
% hObject    handle to Ptotal2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ptotal2 as text
%        str2double(get(hObject,'String')) returns contents of Ptotal2 as a double


% --- Executes during object creation, after setting all properties.
function Ptotal2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ptotal2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GraficarFM.
function GraficarFM_Callback(hObject, eventdata, handles)
% hObject    handle to GraficarFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
Vm = str2double(get(handles.Vm, 'String'));
Vc = str2double(get(handles.Vc, 'String'));
Fm = str2double(get(handles.Fm, 'String'));
Fc = str2double(get(handles.Fc, 'String'));
k = str2double(get(handles.K1, 'String'));
m = (k*Vm)/Fm;
t = 0.0001;
t = linspace(0,t,2000);

A = Vm*cos(2*pi*Fm*t); 

B = Vc*cos(2*pi*Fc*t); 

C = Vc*sin((2*pi*Fc*t)+(m.*sin(2*pi*Fm*t))); 

figure(1)
subplot(311);
plot(t,A);
xlabel('T [s]');
ylabel('A [V]');
title('Moduladora')
grid on;
subplot(312)
plot(t,B);
grid on;
xlabel('T [s]');
ylabel('A [V]');
title('Portadora');
subplot(313);
plot(t,C);
grid on;
xlabel('T [s]');
ylabel('A [V]');
title('Señal En Frecuencia');

%Flip array left to right
Amp1 = fliplr(handles.Amp1);

fre = (Fc-Fm*(length(Amp1)-1)):Fm:(Fc+Fm*(length(Amp1)-1));
%Concatenate arrays along specified dimension
y = cat(2,Amp1,handles.Amp1(2:length(Amp1)));
figure(2)
stem(fre,y)
grid on;
xlabel('F [Hz]');
ylabel('A (V)');
title('Espectro Frecuencia');


% --- Executes on button press in GraficarPM.
function GraficarPM_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
Vm = str2double(get(handles.Vm, 'String'));
Vc = str2double(get(handles.Vc, 'String'));
Fm = str2double(get(handles.Fm, 'String'));
Fc = str2double(get(handles.Fc, 'String'));
k  = str2double(get(handles.K, 'String'));
m  = k*Vm;
t  = 0.0001;
t1 = linspace(0,t,2000);

P = Vm*sin(2*pi*Fm*t1); 
M = Vc*cos(2*pi*Fc*t1); 
FM = Vc*cos((2*pi*Fc*t1)+(m.*sin(2*pi*Fm*t1))); % señal modulada en FM
figure(1)
subplot(311);
plot(t1,P);
xlabel('T [s]');
ylabel('A [V]');
title('Moduladora')
grid on;
subplot(312)
plot(t1,M);
xlabel('T [s]');
ylabel('A [V]');
title('Portadora');
grid on;
subplot(313);
plot(t1,FM);
xlabel('T [s]');
ylabel('A [V]');
title('Modulacion En Fase');
grid on;
%Flip array left to right
Amp2 = fliplr(handles.Amp2); 
fre = (Fc-Fm*(length(Amp2)-1)):Fm:(Fc+Fm*(length(Amp2)-1));
%Concatenate arrays along specified dimension
y = cat(2,Amp2,handles.Amp2(2:length(Amp2)));
figure(2)
stem(fre,y)
grid on;
xlabel('F [Hz]');
ylabel('A [V]');
title('Espectro Frecuencias');



function Res1_Callback(hObject, eventdata, handles)
% hObject    handle to Res1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Res1 as text
%        str2double(get(hObject,'String')) returns contents of Res1 as a double


% --- Executes during object creation, after setting all properties.
function Res1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Res1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrecMax1_Callback(hObject, eventdata, handles)
% hObject    handle to FrecMax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrecMax1 as text
%        str2double(get(hObject,'String')) returns contents of FrecMax1 as a double


% --- Executes during object creation, after setting all properties.
function FrecMax1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrecMax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Res2_Callback(hObject, eventdata, handles)
% hObject    handle to Res2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Res2 as text
%        str2double(get(hObject,'String')) returns contents of Res2 as a double


% --- Executes during object creation, after setting all properties.
function Res2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Res2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrecMax2_Callback(hObject, eventdata, handles)
% hObject    handle to FrecMax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrecMax2 as text
%        str2double(get(hObject,'String')) returns contents of FrecMax2 as a double


% --- Executes during object creation, after setting all properties.
function FrecMax2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrecMax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in Table1.
function Table1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Table1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
