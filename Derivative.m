function varargout = Derivative(varargin)
% DERIVATIVE MATLAB code for Derivative.fig
%      DERIVATIVE, by itself, creates a new DERIVATIVE or raises the existing
%      singleton*.
%
%      H = DERIVATIVE returns the handle to a new DERIVATIVE or the handle to
%      the existing singleton*.
%
%      DERIVATIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DERIVATIVE.M with the given input arguments.
%
%      DERIVATIVE('Property','Value',...) creates a new DERIVATIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Derivative_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Derivative_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Derivative

% Last Modified by GUIDE v2.5 25-Aug-2011 22:03:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Derivative_OpeningFcn, ...
                   'gui_OutputFcn',  @Derivative_OutputFcn, ...
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


% --- Executes just before Derivative is made visible.
function Derivative_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Derivative (see VARARGIN)



clc;


% Choose default command line output for Derivative
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Derivative wait for user response (see UIRESUME)
% uiwait(handles.nonlinearint);


% --- Outputs from this function are returned to the command line.
function varargout = Derivative_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectfile.
function selectfile_Callback(hObject, eventdata, handles)
% hObject    handle to selectfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


clc;

ftxt = fopen('defaultfolder.txt');
defaultfolder = fgetl(ftxt);
fclose(ftxt);

[handles.filename, handles.dirname] = uigetfile('*.xls','Open Data File',defaultfolder);

delete('defaultfolder.txt');
ftxt=fopen('defaultfolder.txt','w');
fprintf(ftxt, '%s',handles.dirname);
fclose(ftxt);




if (~isequal(handles.filename,0) && ~isequal(handles.dirname,0))
   
    
Fileinfo = interpretfilename(handles.filename);
handles.fuel = Fileinfo.fuel;
handles.oxygeninert = Fileinfo.oxygeninert;
handles.temperature = Fileinfo.temperature;
handles.ppcm = Fileinfo.ppcm;
handles.framerate = Fileinfo.framerate;
handles.pressure = Fileinfo.pressure;
handles.phi = Fileinfo.phi;

[~, handles.stanwinfilename] = readfuelinfo(Fileinfo);
[handles.rhoratio, ~] = readstanwin(handles);

boxtext = ['Density Ratio:',num2str(handles.rhoratio)];
set(handles.densityratiobox,'String',boxtext);



 set(handles.filenamebox,'String',handles.filename(1:length(handles.filename)-4));

[handles.num,handles.txt,~] = xlsread([handles.dirname , handles.filename]);

handles.x = handles.num(:,1);

handles.datasize = size(handles.num);

handles.numofpoint = handles.datasize(1);
handles.numofcollum = handles.datasize(2)-1;

handles.Y = handles.num(:,2:handles.datasize(2));


set(handles.numberofcollumnsbox,'String',handles.numofcollum);
set(handles.numberofpointsbox,'String',handles.numofpoint);

set(handles.collumnbox,'String',num2str(1));
set(handles.plottingcontentbox,'String',handles.txt(1,2));

end


handles.output = hObject;
guidata(hObject, handles);



function ordervalue_Callback(hObject, eventdata, handles)
% hObject    handle to ordervalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ordervalue as text
%        str2double(get(hObject,'String')) returns contents of ordervalue as a double
calculatestart_Callback(hObject, eventdata, handles);
return

% --- Executes during object creation, after setting all properties.
function ordervalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ordervalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pointsvalue_Callback(hObject, eventdata, handles)
% hObject    handle to pointsvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pointsvalue as text
%        str2double(get(hObject,'String')) returns contents of pointsvalue as a double
calculatestart_Callback(hObject, eventdata, handles);
return

% --- Executes during object creation, after setting all properties.
function pointsvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pointsvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in fittingmethod.
function fittingmethod_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in fittingmethod 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


selectedradio = get(handles.fittingmethod,'SelectedObject');

if selectedradio==handles.centraldifferencefit
    set(handles.orderbox,'Visible','off');
    set(handles.ordervalue,'Visible','off');
    set(handles.pointsbox,'Visible','off');
    set(handles.pointsvalue,'Visible','off');
elseif selectedradio==handles.globalpolynomialfit
    set(handles.orderbox,'Visible','on');
    set(handles.ordervalue,'Visible','on');
    set(handles.pointsbox,'Visible','off');
    set(handles.pointsvalue,'Visible','off');
elseif selectedradio==handles.localpolynomialfit   
    set(handles.orderbox,'Visible','on');
    set(handles.ordervalue,'Visible','on');
    set(handles.pointsbox,'Visible','on');
    set(handles.pointsvalue,'Visible','on');
elseif selectedradio==handles.piecewisesplinefit   
    set(handles.orderbox,'Visible','off');
    set(handles.ordervalue,'Visible','off');
    set(handles.pointsbox,'Visible','off');
    set(handles.pointsvalue,'Visible','off');
end

calculatestart_Callback(hObject, eventdata, handles);


handles.output = hObject;
guidata(hObject, handles);

    
    


% --- Executes on button press in calculatestart.
function calculatestart_Callback(hObject, eventdata, handles)
% hObject    handle to calculatestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


clc;

if str2double(get(handles.numberofcollumnsbox,'String'))~=0 
    
 
    handles.x_diff = [];
    handles.Y_diff = [];
    handles.dYdx_diff = [];
    handles.kappa = [];
    
selectedradio = get(handles.fittingmethod,'SelectedObject');

for j = 1:1:handles.numofcollum
    
if selectedradio==handles.centraldifferencefit
    
    [handles.x_diff(:,1),handles.Y_diff(:,j),handles.dYdx_diff(:,j)] = centraldifference(handles.x,handles.Y(:,j));

elseif selectedradio==handles.globalpolynomialfit
    
    order = str2double(get(handles.ordervalue,'String'));
    [handles.x_diff(:,1),handles.Y_diff(:,j),handles.dYdx_diff(:,j)] = globalpolynomial(handles.x,handles.Y(:,j),order,0);
    
elseif selectedradio==handles.localpolynomialfit   
    
    order = str2double(get(handles.ordervalue,'String'));
    points = str2double(get(handles.pointsvalue,'String'));
    [handles.x_diff(:,1),handles.Y_diff(:,j),handles.dYdx_diff(:,j)] = localpolynomial(handles.x,handles.Y(:,j),order,points);
    
elseif selectedradio==handles.piecewisesplinefit   
   
    [handles.x_diff(:,1),handles.Y_diff(:,j),handles.dYdx_diff(:,j)] = splinefit(handles.x,handles.Y(:,j));

end

handles.kappa(:,j) = handles.dYdx_diff(:,j)*2./handles.Y_diff(:,j);



end
    
handles.numoffitdata = length(handles.x_diff(:,1));

if get(handles.allbutton,'Value') == 0
ylow = str2double(get(handles.lowradius,'String'));
yhigh = str2double(get(handles.highradius,'String')); 
elseif get(handles.allbutton,'Value') == 1
ylow = handles.Y_diff(1,1);
yhigh = handles.Y_diff(handles.numoffitdata,1);
end
   % use 1st collumn to select data range
    [~,lowindex] = min(abs(handles.Y_diff(:,1)-ylow));
    [~,highindex] = min(abs(handles.Y_diff(:,1)-yhigh));
    
    [~,lowindex_ori] = min(abs(handles.Y(:,1)-ylow));
    [~,highindex_ori] = min(abs(handles.Y(:,1)-yhigh));

 handles.lowindex_ori = lowindex_ori;
 handles.highindex_ori = highindex_ori;
  handles.lowindex = lowindex;
 handles.highindex = highindex;
 
 
clear handles.beta
clear handles.alpha
clear handles.B
handles.beta=[];
handles.alpha=[];
handles.B=[];


clear handles.x_fit
clear handles.Y_fit
clear handles.dYdx_fit
clear handles.kappa_fit
handles.x_fit = [];
handles.Y_fit = [];
handles.dYdx_fit = [];
handles.kappa_fit = [];


clear handles.x_diff_select;
clear handles.Y_diff_select;
clear handles.dYdx_diff_select;
clear handles.kappa_select;
handles.x_diff_select = [];
handles.Y_diff_select = [];
handles.dYdx_diff_select = [];
handles.kappa_select = [];

clear handles.x_select;
clear handles.Y_select;
handles.x_select = [];
handles.Y_select = [];


formularadio = get(handles.formulaselect,'SelectedObject');

% handles.KBLforint
% handles.KBLforderi
% handles.nonlinearint
% handles.nonlinearderi
% handles.linearint
% handles.linearderi


j = str2double(get(handles.collumnbox,'String'));

    
    handles.x_select(:,1) = handles.x(lowindex_ori:highindex_ori,1);
    handles.Y_select(:,j) = handles.Y(lowindex_ori:highindex_ori,j);
    handles.x_diff_select(:,1) = handles.x_diff(lowindex:highindex,1);
    handles.Y_diff_select(:,j) = handles.Y_diff(lowindex:highindex,j);
    handles.dYdx_diff_select(:,j) = handles.dYdx_diff(lowindex:highindex,j); 
    
    handles.kappa_select(:,j) = handles.kappa(lowindex:highindex,j); 
    
if formularadio==handles.dRdt_R

    if get(handles.inputcheckbox,'Value') == 0
    [p,s] = polyfit(log(handles.Y_diff_select(:,j)),log(handles.dYdx_diff_select(:,j)),1);
    handles.beta(1,j) = p(1); handles.alpha(1,j) = 1/(1-handles.beta(1,j));
    handles.B(1,j) = exp(p(2));
    handles.R2(1,j)= 1 - s.normr^2 / norm(log(handles.dYdx_diff_select(:,j))-mean(log(handles.dYdx_diff_select(:,j))))^2;
    elseif get(handles.inputcheckbox,'Value') == 1
    handles.alpha(1,j) = str2double(get(handles.alphainput,'String'));
    handles.beta(1,j) = (handles.alpha(1,j)-1)/handles.alpha(1,j);
    handles.B(1,j) = exp(mean(log(handles.dYdx_diff_select(:,j))-handles.beta(1,j)*log(handles.Y_diff_select(:,j))));
    handles.R2(1,j)= 1 - sum((log(handles.dYdx_diff_select(:,j))-handles.beta(1,j)*log(handles.Y_diff_select(:,j))-log(handles.B(1,j))).^2) / norm(log(handles.dYdx_diff_select(:,j))-mean(log(handles.dYdx_diff_select(:,j))))^2;
    end
   %handles.x_fit = handles.x_diff(1):(handles.x_diff(handles.numoffitdata)-handles.x_diff(1))/1000:handles.x_diff(handles.numoffitdata);
   handles.Y_fit(:,j) = handles.Y_diff(1,j):(handles.Y_diff(handles.numoffitdata,j)-handles.Y_diff(1,j))/1000:handles.Y_diff(handles.numoffitdata,j);
   handles.dYdx_fit(:,j) = handles.B(1,j)*handles.Y_fit(:,j).^handles.beta(1,j);
    
elseif formularadio==handles.dRdt_t
    
    if get(handles.inputcheckbox,'Value') == 0
    [p,s] = polyfit(log(handles.x_diff_select(:,1)),log(handles.dYdx_diff_select(:,j)),1);
    handles.beta(1,j) = p(1); handles.alpha(1,j) = 1+handles.beta(1,j);
    handles.B(1,j) = exp(p(2));   
    handles.R2(1,j)= 1 - s.normr^2 / norm(log(handles.dYdx_diff_select(:,j))-mean(log(handles.dYdx_diff_select(:,j))))^2;
    elseif get(handles.inputcheckbox,'Value') == 1
    handles.alpha(1,j) = str2double(get(handles.alphainput,'String'));
    handles.beta(1,j) = handles.alpha(1,j)-1;
    handles.B(1,j) = exp(mean(log(handles.dYdx_diff_select(:,j))-handles.beta(1,j)*log(handles.x_diff_select(:,1))));
    handles.R2(1,j)= 1 - sum((log(handles.dYdx_diff_select(:,j))-handles.beta(1,j)*log(handles.x_diff_select(:,1))-log(handles.B(1,j))).^2) / norm(log(handles.dYdx_diff_select(:,j))-mean(log(handles.dYdx_diff_select(:,j))))^2;
    end
   
   handles.x_fit(:,j) = handles.x_diff(1,1):(handles.x_diff(handles.numoffitdata,1)-handles.x_diff(1,1))/1000:handles.x_diff(handles.numoffitdata,1);
   %handles.Y_fit(:,j) = handles.Y_diff(1,j):(handles.Y_diff(handles.numoffitdata,j)-handles.Y_diff(1,j))/1000:handles.Y_diff(handles.numoffitdata,j);
   handles.dYdx_fit(:,j) = handles.B(1,j)*handles.x_fit(:,1).^handles.beta(1,j);
   
   
elseif formularadio==handles.R_r0plust_a
   
    if get(handles.inputcheckbox,'Value') == 0 && get(handles.givenr0,'Value') == 0
    bet0= [0 1000 1.5];
    [bet,~,~] = nlinfit(handles.x_select(:,1),handles.Y_select(:,j),@rtformula,bet0);
    handles.r0(1,j) = bet(1);
    handles.A(1,j) = bet(2);
    handles.alpha(1,j) = bet(3);
    handles.beta(1,j) = handles.alpha(1,j)-1;
    handles.B(1,j) = handles.A(1,j)*handles.alpha(1,j);
    handles.R2(1,j)= 1 - sum((handles.Y_select(:,j) - handles.r0(1,j)-handles.A(1,j)*handles.x_select(:,1).^handles.alpha(1,j)).^2) / norm(handles.Y_select(:,j)-mean(handles.Y_select(:,j)))^2;
   elseif get(handles.inputcheckbox,'Value') == 1 && get(handles.givenr0,'Value') == 0
    handles.alpha(1,j) = str2double(get(handles.alphainput,'String'));
    handles.beta(1,j) = handles.alpha(1,j)-1;
    [p,s] = polyfit(handles.x_select(:,1).^handles.alpha(1,j),handles.Y_select(:,j),1);
    handles.A(1,j) = p(1);
    handles.r0(1,j) = p(2);
    handles.R2(1,j)= 1 - s.normr^2 / norm(handles.Y_select(:,j)-mean(handles.Y_select(:,j)))^2;
    handles.B(1,j) = handles.A(1,j)*handles.alpha(1,j);
   elseif get(handles.inputcheckbox,'Value') == 0 && get(handles.givenr0,'Value') == 1
    handles.r0(1,j) = str2double(get(handles.r0box,'String'));
    [p,s] = polyfit(log(handles.x_select(:,1)),log(handles.Y_select(:,j)-handles.r0(1,j)),1);    
    handles.A(1,j) = exp(p(2));
    handles.alpha(1,j) = p(1);
    handles.beta(1,j) = handles.alpha(1,j)-1;
    handles.B(1,j) = handles.A(1,j)*handles.alpha(1,j);
    handles.R2(1,j)= 1 - s.normr^2 / norm(log(handles.Y_select(:,j)-handles.r0(1,j))-mean(log(handles.Y_select(:,j)-handles.r0(1,j))))^2;
   elseif get(handles.inputcheckbox,'Value') == 1 && get(handles.givenr0,'Value') == 1
    handles.alpha(1,j) = str2double(get(handles.alphainput,'String'));                
    handles.r0(1,j) = str2double(get(handles.r0box,'String'));                
    handles.A(1,j) = exp(mean(log(handles.Y_select(:,j)-handles.r0(1,j))-handles.alpha(1,j).*log(handles.x_select(:,1))));         
    handles.R2(1,j)= 1 - sum((handles.Y_select(:,j) - handles.r0(1,j)-handles.A(1,j)*handles.x_select(:,1).^handles.alpha(1,j)).^2) / norm(handles.Y_select(:,j)-mean(handles.Y_select(:,j)))^2;          
    handles.beta(1,j) = handles.alpha(1,j)-1;
    handles.B(1,j) = handles.A(1,j)*handles.alpha(1,j);                
   end
   
   
   handles.x_fit(:,j) = handles.x_diff(1,1):(handles.x_diff(handles.numoffitdata,1)-handles.x_diff(1,1))/1000:handles.x_diff(handles.numoffitdata,1);
   handles.Y_fit(:,j) = handles.r0(1,j) + handles.A(1,j)*handles.x_fit(:,1).^handles.alpha(1,j);
   handles.dYdx_fit(:,j) = handles.B(1,j)*handles.x_fit(:,1).^(handles.beta(1,j));
   
   
   
elseif formularadio==handles.linearint
       
   [handles.flamespeed(1,j),handles.markstein(1,j),handles.constant(1,j)]= linearintfit(handles.Y_select(:,j), handles.x_select(:,1));

   handles.Y_fit(:,j) = 0.1*handles.Y_diff(1,j):(handles.Y_diff(handles.numoffitdata,j)-handles.Y_diff(1,j))/1000:handles.Y_diff(handles.numoffitdata,j)*100;
   handles.x_fit(:,j) = (handles.Y_fit(:,j)-handles.constant(1,j) + 2*handles.markstein(1,j)*log(handles.Y_fit(:,j)))/handles.flamespeed(1,j);
   handles.dYdx_fit(:,j) = handles.flamespeed(1,j)./(1+2*handles.markstein(1,j)./handles.Y_fit(:,j));
   handles.kappa_fit(:,j) = handles.dYdx_fit(:,j)*2./handles.Y_fit(:,j);
   
elseif formularadio==handles.linearderi
       
   [handles.flamespeed(1,j),handles.markstein(1,j)]= linearderifit(handles.dYdx_diff_select(:,j), handles.kappa_select(:,j));
   handles.Y_fit(:,j) = 0.1*handles.Y_diff(1,j):(handles.Y_diff(handles.numoffitdata,j)-handles.Y_diff(1,j))/1000:handles.Y_diff(handles.numoffitdata,j)*100;
   %handles.x_fit(:,j) = (handles.Y_fit(:,j)-handles.constant(1,j) + 2*handles.markstein(1,j)*log(handles.Y_fit(:,j)))/handles.flamespeed(1,j);
   handles.dYdx_fit(:,j) = handles.flamespeed(1,j)./(1+2*handles.markstein(1,j)./handles.Y_fit(:,j));
   handles.kappa_fit(:,j) = handles.dYdx_fit(:,j)*2./handles.Y_fit(:,j);
   
elseif formularadio==handles.nonlinearint
   
   [beta0 beta1 beta2] = linearintfit(handles.Y_select(:,j), handles.x_select(:,1));
   [handles.flamespeed(1,j),handles.markstein(1,j),handles.constant(1,j)] = nonlinearintfit(handles.Y_select(:,j),handles.x_select(:,1), [beta0 beta1 beta2]);
    %handles.dYdx_fit(:,j) = handles.dYdx_diff(1,j):(handles.dYdx_diff(handles.numoffitdata,j)-handles.dYdx_diff(1,j))/1000:handles.dYdx_diff(handles.numoffitdata,j);
    handles.dYdx_fit(:,j) = 0:max(handles.dYdx_diff(1,j))/1000:3*max(handles.dYdx_diff(1,j));
    handles.Y_fit(:,j)=-2*handles.markstein(1,j)./((handles.dYdx_fit(:,j)/handles.flamespeed(1,j)).*log(handles.dYdx_fit(:,j)/handles.flamespeed(1,j)));
    tau=findtau(handles.Y_fit((length(handles.Y_fit(:,j))),j)/(-2*handles.markstein(1,j)) );
    ttilde=1/(tau^2*log(tau))-real(expint(2*log(tau)))+handles.constant(1,j);
    time_last = ((-2*handles.markstein(1,j)/handles.flamespeed(1,j))*ttilde);
    handles.x_fit(:,j) = zeros(1,length(handles.Y_fit(:,j)));
    handles.x_fit(length(handles.Y_fit(:,j)),j) = time_last;
    for i=(length(handles.x_fit(:,j))-1):-1:1
        dtime = (handles.Y_fit((i+1),j)-handles.Y_fit((i),j))*(1./handles.dYdx_fit(i+1,j)+1./handles.dYdx_fit(i,j))/2;
        handles.x_fit(i,j) = handles.x_fit(i+1,j) - dtime;
    end
    handles.kappa_fit(:,j) = handles.dYdx_fit(:,j)*2./handles.Y_fit(:,j);
    
elseif formularadio==handles.nonlinearderi
   
   
    
    [beta0 beta1 ~] = linearintfit(handles.Y_select(:,j), handles.x_select(:,1));
    [handles.flamespeed(1,j),handles.markstein(1,j)] = nonlinearderifit(handles.dYdx_diff_select(:,j), handles.kappa_select(:,j), [beta0 beta1]);
    %handles.dYdx_fit(:,j) = handles.dYdx_diff(1,j):(handles.dYdx_diff(handles.numoffitdata,j)-handles.dYdx_diff(1,j))/1000:handles.dYdx_diff(handles.numoffitdata,j);
    handles.dYdx_fit(:,j) = 0:max(handles.dYdx_diff(1,j))/1000:3*max(handles.dYdx_diff(1,j));
    handles.kappa_fit(:,j) = -(handles.flamespeed(1,j)./(2*handles.markstein(1,j))).*((handles.dYdx_fit(:,j)./handles.flamespeed(1,j)).^2).*log((handles.dYdx_fit(:,j)./handles.flamespeed(1,j)).^2);
    handles.Y_fit(:,j) = handles.dYdx_fit(:,j)*2./handles.kappa_fit(:,j);
    
   
elseif formularadio==handles.KBLforint
    
    [beta0 beta1 beta2]= linearintfit(handles.Y_select(:,j), handles.x_select(:,1));
    [handles.flamespeed(1,j),handles.markstein(1,j),handles.constant(1,j)] = KBLforint(handles.Y_select(:,j),handles.x_select(:,1), [beta0 beta1 beta2]);
    handles.Y_fit(:,j) = 0.1*handles.Y_diff(1,j):(handles.Y_diff(handles.numoffitdata,j)-handles.Y_diff(1,j))/1000:handles.Y_diff(handles.numoffitdata,j)*100;
    handles.x_fit(:,j) = (handles.Y_fit(:,j) + 2*handles.markstein(1,j)*log(handles.Y_fit(:,j)) - 4*handles.markstein(1,j)^2./handles.Y_fit(:,j) - 8*handles.markstein(1,j)^3./(3*handles.Y_fit(:,j)) - handles.constant(1,j))/handles.flamespeed(1,j);
    handles.dYdx_fit(:,j) = handles.flamespeed(1,j)./( 1+ 2*handles.markstein(1,j)./handles.Y_fit(:,j) + 4*handles.markstein(1,j)^2./handles.Y_fit(:,j).^2 + (16/3)*handles.markstein(1,j)^3./handles.Y_fit(:,j).^3);
    handles.kappa_fit(:,j) = handles.dYdx_fit(:,j)*2./handles.Y_fit(:,j);

elseif formularadio==handles.KBLforderi
    
    [beta0 beta1 ~] = linearintfit(handles.Y_select(:,j), handles.x_select(:,1));
    [handles.flamespeed(1,j),handles.markstein(1,j)] = KBLforderi(handles.Y_diff_select(:,j),handles.dYdx_diff_select(:,j),[beta0 beta1]);
    handles.Y_fit(:,j) = 0.1*handles.Y_diff(1,j):(handles.Y_diff(handles.numoffitdata,j)-handles.Y_diff(1,j))/1000:handles.Y_diff(handles.numoffitdata,j)*100;
    handles.dYdx_fit(:,j) = handles.flamespeed(1,j)./( 1+ 2*handles.markstein(1,j)./handles.Y_fit(:,j) + 4*handles.markstein(1,j)^2./handles.Y_fit(:,j).^2 + (16/3)*handles.markstein(1,j)^3./handles.Y_fit(:,j).^3);
    handles.kappa_fit(:,j) = handles.dYdx_fit(:,j)*2./handles.Y_fit(:,j);
    
end
    

   makeplot(handles);

end


handles.output = hObject;
guidata(hObject, handles);



function makeplot(handles)

    cla(handles.axes1); 
    cla(handles.axes2); 
    cla(handles.axes3); 
    cla(handles.axes4); 
    
    
    plottingcollumn = str2double(get(handles.collumnbox,'String'));
    shape = get(handles.shapebox,'String');
    color = get(handles.colorbox,'String');
    linecolor = get(handles.linecolorbox,'String');
    jump = str2double(get(handles.jumpstep,'String'));
    formularadio = get(handles.formulaselect,'SelectedObject');
    
    
%   handles.lowindex_ori = lowindex_ori;
% handles.highindex_ori = highindex_ori;
%  handles.lowindex = lowindex;
% handles.highindex = highindex;
 
    
      axes(handles.axes1);
     cmd = ['plot(handles.x(handles.highindex_ori:jump:handles.numofpoint,1),handles.Y(handles.highindex_ori:jump:handles.numofpoint,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
      cmd = ['plot(handles.x(1:jump:handles.lowindex_ori,1),handles.Y(1:jump:handles.lowindex_ori,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['plot(handles.x(handles.lowindex_ori:jump:handles.highindex_ori,1),handles.Y(handles.lowindex_ori:jump:handles.highindex_ori,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
     eval(cmd);
     
   %       axes(handles.axes1);
   %  cmd = ['plot(handles.x(1:jump:handles.numofpoint,1),handles.Y(1:jump:handles.numofpoint,plottingcollumn),''',...
   %      shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
   %  eval(cmd);
   %  hold on;
   %  cmd = ['plot(handles.x_select(1:jump:numofselectpoint,1),handles.Y_select(1:jump:numofselectpoint,plottingcollumn),''',...
    %     shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
    %eval(cmd);
     ylabel('R (cm)');
     xlabel('t (s)');
set(handles.axes1,'Xlim',[handles.x(1,1)-0.1*abs(handles.x(handles.numofpoint,1)-handles.x(1,1)),handles.x(handles.numofpoint,1)+0.1*abs(handles.x(handles.numofpoint,1)-handles.x(1,1))]);
set(handles.axes1,'Ylim',[handles.Y(1,plottingcollumn)-0.1*abs(handles.Y(handles.numofpoint,plottingcollumn)-handles.Y(1,plottingcollumn)),handles.Y(handles.numofpoint,plottingcollumn)+0.1*abs(handles.Y(handles.numofpoint,plottingcollumn)-handles.Y(1,plottingcollumn))]);
    if ((formularadio==handles.R_r0plust_a) || (formularadio==handles.linearint) || (formularadio==handles.nonlinearint) || (formularadio==handles.KBLforint))
      cmd = ['plot(handles.x_fit(:,plottingcollumn),handles.Y_fit(:,plottingcollumn),''',...
         '-',''',''linewidth'',2,''color'',''',linecolor,''');']; 
     eval(cmd);
    end
    hold off;

    
    axes(handles.axes2);
     cmd = ['plot(handles.x_diff(handles.highindex:jump:handles.numoffitdata,1),handles.dYdx_diff(handles.highindex:jump:handles.numoffitdata ,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['plot(handles.x_diff(1:jump:handles.lowindex,1),handles.dYdx_diff(1:jump:handles.lowindex ,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['plot(handles.x_diff(handles.lowindex:jump:handles.highindex,1),handles.dYdx_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
     eval(cmd);  
     if formularadio==handles.dRdt_t || formularadio==handles.R_r0plust_a || formularadio==handles.linearint || formularadio==handles.nonlinearint || formularadio==handles.KBLforint
     cmd = ['plot(handles.x_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
         '-',''',''linewidth'',2,''color'',''',linecolor,''');']; 
     eval(cmd);
     end
 %       axes(handles.axes2);
 %    cmd = ['plot(handles.x_diff(1:jump:handles.numoffitdata,1),handles.dYdx_diff(1:jump:handles.numoffitdata ,plottingcollumn),''',...
 %        shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
 %    eval(cmd);
 %    hold on;
 %    cmd = ['plot(handles.x_diff_select(:,1),handles.dYdx_diff_select(:,plottingcollumn),''',...
 %        shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
 %    eval(cmd);  
 %    if formularadio==handles.dRdt_t || formularadio==handles.R_r0plust_a || formularadio==handles.linearint || formularadio==handles.nonlinearint || formularadio==handles.KBLforint
 %    cmd = ['plot(handles.x_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
 %        '-',''',''linewidth'',2,''color'',''',linecolor,''');']; 
 %    eval(cmd);
 %    end
     ylabel('U (cm/s)');
     xlabel('t (s)');
set(handles.axes2,'Xlim',[handles.x(1,1)-0.1*abs(handles.x(handles.numofpoint,1)-handles.x(1,1)),handles.x(handles.numofpoint,1)+0.1*abs(handles.x(handles.numofpoint,1)-handles.x(1,1))]);
set(handles.axes2,'Ylim',[min(handles.dYdx_diff(:,plottingcollumn))-0.1*(abs(max(handles.dYdx_diff(:,plottingcollumn))-min(handles.dYdx_diff(:,plottingcollumn)))),max(handles.dYdx_diff(:,plottingcollumn))+0.1*(abs(max(handles.dYdx_diff(:,plottingcollumn))-min(handles.dYdx_diff(:,plottingcollumn))))]);
     hold off;
    
     
     handles.lowindex
     axes(handles.axes3);
     cmd = ['plot(handles.Y_diff(handles.lowindex:jump:handles.numoffitdata,plottingcollumn),handles.dYdx_diff(handles.lowindex:jump:handles.numoffitdata ,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
      cmd = ['plot(handles.Y_diff(1:jump:handles.lowindex,plottingcollumn),handles.dYdx_diff(1:jump:handles.lowindex ,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['plot(handles.Y_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),handles.dYdx_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
     eval(cmd);
     ylabel('U (cm/s)');
     xlabel('R (cm)');
     if formularadio==handles.dRdt_R || formularadio==handles.R_r0plust_a ...
             || formularadio==handles.linearint || formularadio==handles.nonlinearint || formularadio==handles.KBLforint...
             || formularadio==handles.linearderi || formularadio==handles.nonlinearderi || formularadio==handles.KBLforderi
     cmd = ['plot(handles.Y_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
         '-',''',''linewidth'',2,''color'',''',linecolor,''');'];

     eval(cmd);
     end
%      axes(handles.axes3);
%      cmd = ['plot(handles.Y_diff(1:jump:handles.numoffitdata,plottingcollumn),handles.dYdx_diff(1:jump:handles.numoffitdata ,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
%      eval(cmd);
%      hold on;
%      cmd = ['plot(handles.Y_diff_select(:,plottingcollumn),handles.dYdx_diff_select(:,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
%      eval(cmd);
%      ylabel('U (cm/s)');
%      xlabel('R (cm)');
%      if formularadio==handles.dRdt_R || formularadio==handles.R_r0plust_a ...
%              || formularadio==handles.linearint || formularadio==handles.nonlinearint || formularadio==handles.KBLforint...
%              || formularadio==handles.linearderi || formularadio==handles.nonlinearderi || formularadio==handles.KBLforderi
%      cmd = ['plot(handles.Y_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
%          '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
% 
%      eval(cmd);
%      end     
      hold off;
set(handles.axes3,'Xlim',[handles.Y(1,plottingcollumn)-0.1*abs(handles.Y(handles.numofpoint,plottingcollumn)-handles.Y(1,plottingcollumn)),handles.Y(handles.numofpoint,plottingcollumn)+0.1*abs(handles.Y(handles.numofpoint,plottingcollumn)-handles.Y(1,plottingcollumn))]);
set(handles.axes3,'Ylim',[min(handles.dYdx_diff(:,plottingcollumn))-0.1*(abs(max(handles.dYdx_diff(:,plottingcollumn))-min(handles.dYdx_diff(:,plottingcollumn)))),max(handles.dYdx_diff(:,plottingcollumn))+0.1*(abs(max(handles.dYdx_diff(:,plottingcollumn))-min(handles.dYdx_diff(:,plottingcollumn))))]);
     
     
     axes(handles.axes4);
       
     if formularadio==handles.dRdt_t
         
     cmd = ['loglog(handles.x_diff(handles.highindex:jump:handles.numoffitdata,1),handles.dYdx_diff(handles.highindex:jump:handles.numoffitdata ,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['loglog(handles.x_diff(1:jump:handles.lowindex,1),handles.dYdx_diff(1:jump:handles.lowindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['loglog(handles.x_diff(handles.lowindex:jump:handles.highindex,1),handles.dYdx_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
     eval(cmd);
     cmd = ['loglog(handles.x_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
         '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
     eval(cmd);
%      cmd = ['loglog(handles.x_diff(1:jump:handles.numoffitdata,1),handles.dYdx_diff(1:jump:handles.numoffitdata ,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
%      eval(cmd);
%      hold on;
%      cmd = ['loglog(handles.x_diff_select(:,1),handles.dYdx_diff_select(:,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
%      eval(cmd);
%      cmd = ['loglog(handles.x_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
%          '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
%      eval(cmd);
     ylabel('U_a_v  (cm/sec)');
     xlabel('t (sec)'); 

     elseif formularadio==handles.dRdt_R
     cmd = ['loglog(handles.Y_diff(handles.highindex:jump:handles.numoffitdata,plottingcollumn),handles.dYdx_diff(handles.highindex:jump:handles.numoffitdata,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['loglog(handles.Y_diff(1:jump:handles.lowindex,plottingcollumn),handles.dYdx_diff(1:jump:handles.lowindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['loglog(handles.Y_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),handles.dYdx_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
     eval(cmd);
     cmd = ['loglog(handles.Y_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
         '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
     eval(cmd);
%           cmd = ['loglog(handles.Y_diff(1:jump:handles.numoffitdata,plottingcollumn),handles.dYdx_diff(1:jump:handles.numoffitdata ,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
%      eval(cmd);
%      hold on;
%      cmd = ['loglog(handles.Y_diff_select(:,plottingcollumn),handles.dYdx_diff_select(:,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
%      eval(cmd);
%      cmd = ['loglog(handles.Y_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
%          '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
%      eval(cmd);
     ylabel('U_a_v (cm/sec)');
     xlabel('R (cm)');  
     
     
     elseif formularadio==handles.linearint || formularadio==handles.nonlinearint || formularadio==handles.KBLforint...
             || formularadio==handles.linearderi || formularadio==handles.nonlinearderi || formularadio==handles.KBLforderi
     cmd = ['plot(handles.kappa(handles.highindex:jump:handles.numoffitdata,plottingcollumn),handles.dYdx_diff(handles.highindex:jump:handles.numoffitdata,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['plot(handles.kappa(1:jump:handles.lowindex,plottingcollumn),handles.dYdx_diff(1:jump:handles.lowindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
     eval(cmd);
     hold on;
     cmd = ['plot(handles.kappa(handles.lowindex:jump:handles.highindex,plottingcollumn),handles.dYdx_diff(handles.lowindex:jump:handles.highindex,plottingcollumn),''',...
         shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
     eval(cmd);
     cmd = ['plot(handles.kappa_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
         '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
     eval(cmd);
%      cmd = ['plot(handles.kappa(1:jump:handles.numoffitdata,plottingcollumn),handles.dYdx_diff(1:jump:handles.numoffitdata ,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''w'');'];
%      eval(cmd);
%      hold on;
%      cmd = ['plot(handles.kappa_select(:,plottingcollumn),handles.dYdx_diff_select(:,plottingcollumn),''',...
%          shape,''',''linewidth'',2,''Markeredgecolor'',''',color,''',''Markerfacecolor'',''',color,''');'];
%      eval(cmd);
%      cmd = ['plot(handles.kappa_fit(:,plottingcollumn),handles.dYdx_fit(:,plottingcollumn),''',...
%          '-',''',''linewidth'',2,''color'',''',linecolor,''');'];
%      eval(cmd);
     ylabel('U_a_v (cm/sec)');
     xlabel('K (cm)');
set(handles.axes4,'Xlim',[0, handles.kappa(4,plottingcollumn)]);  
set(handles.axes4,'Ylim',[min(handles.dYdx_diff(:,plottingcollumn))-0.1*(abs(max(handles.dYdx_diff(:,plottingcollumn))-min(handles.dYdx_diff(:,plottingcollumn)))),max(handles.dYdx_diff(:,plottingcollumn))+0.01*(abs(max(handles.dYdx_diff(:,plottingcollumn))-min(handles.dYdx_diff(:,plottingcollumn))))]);

     end
     hold off;
     
     
     
     
     
     
     
     if formularadio==handles.dRdt_t
     formula = ['U = ',num2str(handles.B(1,plottingcollumn)),'*t^',num2str(handles.beta(1,plottingcollumn))];
     set(handles.formulabox,'String',formula);
     elseif formularadio==handles.dRdt_R
     formula = ['U = ',num2str(handles.B(1,plottingcollumn)),'*R^',num2str(handles.beta(1,plottingcollumn))];
     set(handles.formulabox,'String',formula);
     elseif  formularadio==handles.R_r0plust_a
     formula = ['R = ',num2str(handles.r0(1,plottingcollumn)),'+',num2str(handles.A(1,plottingcollumn)),'*t^',num2str(handles.alpha(1,plottingcollumn))];
     set(handles.formulabox,'String',formula);
     elseif formularadio==handles.linearint
     formula = 'linear extrapolation: R-t';
     set(handles.formulabox,'String',formula);     
     elseif formularadio==handles.linearderi
     formula = 'linear extrapolation: U-R';
     set(handles.formulabox,'String',formula);      
     elseif formularadio==handles.nonlinearint
     formula = 'nonlinear extrapolation: R-t';
     set(handles.formulabox,'String',formula);           
     elseif formularadio==handles.nonlinearderi
     formula = 'nonlinear extrapolation: U-R';
     set(handles.formulabox,'String',formula);           
     elseif formularadio==handles.KBLforint
     formula = 'KBL extrapolation: R-t';
     set(handles.formulabox,'String',formula);               
     elseif formularadio==handles.KBLforderi
     formula = 'KBL extrapolation: U-R';
     set(handles.formulabox,'String',formula);              
     end
     
     if formularadio==handles.dRdt_t || formularadio==handles.dRdt_R || formularadio==handles.R_r0plust_a
     box1text = ['Alpha = ',num2str(handles.alpha(1,plottingcollumn))];
     set(handles.box1,'String',box1text);
     box2text = ['R2 =',num2str(handles.R2(1,plottingcollumn))];
     set(handles.box2,'String',box2text);
     else
     box1text = ['Flamespeed = ',num2str(handles.flamespeed(1,plottingcollumn)/handles.rhoratio)];
     set(handles.box1,'String',box1text);
     box2text = ['Markerstein =',num2str(handles.markstein(1,plottingcollumn))];
     set(handles.box2,'String',box2text);
     end
     
return


% --- Executes on button press in minusbox.
function minusbox_Callback(hObject, eventdata, handles)
% hObject    handle to minusbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if str2double(get(handles.numberofcollumnsbox,'String'))~=0 
    
    i = str2double(get(handles.collumnbox,'String'));

    if i<=1
        i=1;
    elseif i>str2double(get(handles.numberofcollumnsbox,'String'))
        i=str2double(get(handles.numberofcollumnsbox,'String'));
    else
        i=i-1;
    end
    
set(handles.collumnbox,'String',num2str(i));
set(handles.plottingcontentbox,'String',handles.txt(1,1+i)); 

calculatestart_Callback(hObject, eventdata, handles);
 %makeplot(handles);
 
end

return

% handles.output = hObject;
% guidata(hObject, handles);

% --- Executes on button press in plusbutton.
function plusbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if str2double(get(handles.numberofcollumnsbox,'String'))~=0 
    
 

 i = str2double(get(handles.collumnbox,'String'));
    if i<1
        i=1;
    elseif i>=str2double(get(handles.numberofcollumnsbox,'String'))
        i=str2double(get(handles.numberofcollumnsbox,'String'));
    else
        i=i+1;
    end
    
    
set(handles.collumnbox,'String',num2str(i));
set(handles.plottingcontentbox,'String',handles.txt(1,1+i));

     %makeplot(handles);
     calculatestart_Callback(hObject, eventdata, handles);


end



return



% --- Executes during object creation, after setting all properties.
function numberofcollumnsbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numberofcollumnsbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function numberofpointsbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numberofpointsbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in createfig.
function createfig_Callback(hObject, eventdata, handles)
% hObject    handle to createfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.newfig1=figure('InvertHardcopy','off','Color',[1 1 1]);
% Create axes
handles.newaxes1=axes('Parent',handles.newfig1,'Position',[0.1854 0.1476 0.6562 0.7774],...
    'FontSize',14,...
    'FontName','Times New Roman');
box('on');
hold('all');
xlabel('t (s)');
ylabel('R (cm)');
set(handles.newaxes1,'Xlim',get(handles.axes1,'Xlim'));
set(handles.newaxes1,'Ylim',get(handles.axes1,'Ylim'));

handles.newfig2=figure('InvertHardcopy','off','Color',[1 1 1]);
% Create axes
handles.newaxes2=axes('Parent',handles.newfig2,'Position',[0.1854 0.1476 0.6562 0.7774],...
    'FontSize',14,...
    'FontName','Times New Roman');
box('on');
hold('all');
xlabel('t (s)');
ylabel('dR/dt (cm/s)');
set(handles.newaxes2,'Xlim',get(handles.axes2,'Xlim'));
set(handles.newaxes2,'Ylim',get(handles.axes2,'Ylim'));

handles.newfig3=figure('InvertHardcopy','off','Color',[1 1 1]);
% Create axes
handles.newaxes3=axes('Parent',handles.newfig3,'Position',[0.1854 0.1476 0.6562 0.7774],...
    'FontSize',14,...
    'FontName','Times New Roman');
box('on');
hold('all');
xlabel('R (cm)');
ylabel('dR/dt (cm/s)');
set(handles.newaxes3,'Xlim',get(handles.axes3,'Xlim'));
set(handles.newaxes3,'Ylim',get(handles.axes3,'Ylim'));

handles.newfig4=figure('InvertHardcopy','off','Color',[1 1 1]);
% Create axes
handles.newaxes4=axes('Parent',handles.newfig4,'Position',[0.1854 0.1476 0.6562 0.7774],...
    'FontSize',14,...
    'FontName','Times New Roman');
box('on');
hold('all');
set(handles.newaxes4,'Xlim',get(handles.axes4,'Xlim'));
set(handles.newaxes4,'Ylim',get(handles.axes4,'Ylim'));

formularadio = get(handles.formulaselect,'SelectedObject');

     if formularadio==handles.dRdt_t
        xlabel('t (sec)');
        ylabel('U_a_v (cm/sec)');
        set(handles.newaxes4,'yscale','log');
        set(handles.newaxes4,'xscale','log');
     elseif formularadio==handles.dRdt_R
        xlabel('R_a_v (cm)');
        ylabel('U_a_v (cm/sec)');
        set(handles.newaxes4,'yscale','log');
        set(handles.newaxes4,'xscale','log');
     else 
        xlabel('\kappa');
        ylabel('S_b (cm/sec)'); 
     end
     
     
     



handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in addtofig.
function addtofig_Callback(hObject, eventdata, handles)
% hObject    handle to addtofig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


copyobj(get(handles.axes1,'Children'),handles.newaxes1);
copyobj(get(handles.axes2,'Children'),handles.newaxes2);
copyobj(get(handles.axes3,'Children'),handles.newaxes3);
copyobj(get(handles.axes4,'Children'),handles.newaxes4);
findfigs;

handles.output = hObject;
guidata(hObject, handles);



% --- Executes on button press in outputbutton.
function outputbutton_Callback(hObject, eventdata, handles)
% hObject    handle to outputbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if str2double(get(handles.numberofcollumnsbox,'String'))~=0 
    
    
DefaultName = [handles.dirname,'Diff Results_',handles.filename];
[FileName,PathName,FilterIndex] = uiputfile({'*.xls';'*.csv'},'Output File',DefaultName);

if FilterIndex == 0
    
else
    
    sizeoftext = size(handles.txt);
    
   nameunit1 = handles.txt;
   writematrix1 = [handles.x_diff(:,1),handles.Y_diff];
   xlswrite([PathName,FileName],nameunit1,1,'A1');
   xlswrite([PathName,FileName],writematrix1,1,['A', num2str(sizeoftext(1)+1)]);
   
   nameunit2 = handles.txt;
   for j=2:1:sizeoftext(2)
       nameunit2(1,j) = {['d',char(nameunit2(1,j)),'dx']};
   end
   writematrix2 = [handles.x_diff(:,1),handles.dYdx_diff];
   xlswrite([PathName,FileName],nameunit2,2,'A1');
   xlswrite([PathName,FileName],writematrix2,2,['A', num2str(sizeoftext(1)+1)]);
   
end

end

handles.output = hObject;
guidata(hObject, handles);



function shapebox_Callback(hObject, eventdata, handles)
% hObject    handle to shapebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shapebox as text
%        str2double(get(hObject,'String')) returns contents of shapebox as a double


 %makeplot(handles);
 
 calculatestart_Callback(hObject, eventdata, handles);
    
    return
    
    
% --- Executes during object creation, after setting all properties.
function shapebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shapebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowradius_Callback(hObject, eventdata, handles)
% hObject    handle to lowradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowradius as text
%        str2double(get(hObject,'String')) returns contents of lowradius as a double
calculatestart_Callback(hObject, eventdata, handles);
return

% --- Executes during object creation, after setting all properties.
function lowradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function highradius_Callback(hObject, eventdata, handles)
% hObject    handle to highradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highradius as text
%        str2double(get(hObject,'String')) returns contents of highradius as a double
calculatestart_Callback(hObject, eventdata, handles);
return

% --- Executes during object creation, after setting all properties.
function highradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function colorbox_Callback(hObject, eventdata, handles)
% hObject    handle to colorbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of colorbox as text
%        str2double(get(hObject,'String')) returns contents of colorbox as a double

 %makeplot(handles);
 
 calculatestart_Callback(hObject, eventdata, handles);
 
 return;

% --- Executes during object creation, after setting all properties.
function colorbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jumpstep_Callback(hObject, eventdata, handles)
% hObject    handle to jumpstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jumpstep as text
%        str2double(get(hObject,'String')) returns contents of jumpstep as a double


 %makeplot(handles);
 calculatestart_Callback(hObject, eventdata, handles);
 
 return;
% --- Executes during object creation, after setting all properties.
function jumpstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jumpstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in formulaselect.
function formulaselect_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in formulaselect 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

calculatestart_Callback(hObject, eventdata, handles);


handles.output = hObject;
guidata(hObject, handles);



function alphainput_Callback(hObject, eventdata, handles)
% hObject    handle to alphainput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphainput as text
%        str2double(get(hObject,'String')) returns contents of alphainput as a double

calculatestart_Callback(hObject, eventdata, handles)

return

% --- Executes during object creation, after setting all properties.
function alphainput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphainput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inputcheckbox.
function inputcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to inputcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of inputcheckbox

if get(handles.inputcheckbox,'Value') == 1
set(handles.alphainput,'Visible','on');
elseif get(handles.inputcheckbox,'Value') == 0
set(handles.alphainput,'Visible','off');
end

calculatestart_Callback(hObject, eventdata, handles)

return


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

formularadio = get(handles.formulaselect,'SelectedObject');

if formularadio==handles.dRdt_t || formularadio==handles.dRdt_R || formularadio==handles.R_r0plust_a

FileName = ['AccSummary_',date,'.xls'];
PathName = 'D:\Projects\';

    if exist([PathName,FileName]) == 2
       [num,txt,~] = xlsread([PathName,FileName],1);
       numsize = size(num); txtsize = size(txt);
       line = txtsize(1)+1;
    else
       header=[cellstr('File') cellstr('Pres') cellstr('phi')  cellstr('formula') cellstr('Rmin') cellstr('Rmax') ...
        cellstr('alpha')  cellstr('beta')  cellstr('B')   cellstr('R2')];
        xlswrite([PathName,FileName],header,1);
      line = 2;
    end
    
    range = ['A',num2str(line)];

       plottingcollumn = str2double(get(handles.collumnbox,'String'));
       
    formularadio = get(handles.formulaselect,'SelectedObject');
    if formularadio==handles.dRdt_R
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('U-R')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.alpha(1,plottingcollumn) handles.beta(1,plottingcollumn) handles.B(1,plottingcollumn) handles.R2(1,plottingcollumn)];
    elseif formularadio==handles.dRdt_t
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('U-t')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.alpha(1,plottingcollumn) handles.beta(1,plottingcollumn) handles.B(1,plottingcollumn) handles.R2(1,plottingcollumn)];
    end
    
    xlswrite([PathName,FileName],content,1,range);
    
else
    
FileName = ['Flamespeed_',date,'.xls'];
PathName = 'D:\Projects\';

    if exist([PathName,FileName]) == 2
       [num,txt,~] = xlsread([PathName,FileName],1);
       numsize = size(num); txtsize = size(txt);
       line = txtsize(1)+1;
    else
       header=[cellstr('File') cellstr('Pres') cellstr('phi') cellstr('formula') cellstr('Rmin') cellstr('Rmax') ...
        cellstr('Flamespeed')  cellstr('Markerstein Length')];
        xlswrite([PathName,FileName],header,1);
      line = 2;
    end
    
    range = ['A',num2str(line)];

       plottingcollumn = str2double(get(handles.collumnbox,'String'));
       
    formularadio = get(handles.formulaselect,'SelectedObject');
    if formularadio==handles.linearint
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('linear:R-t')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.flamespeed(1,plottingcollumn) handles.markstein(1,plottingcollumn) ];
    elseif formularadio==handles.linearderi
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('linear:U-R')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.flamespeed(1,plottingcollumn) handles.markstein(1,plottingcollumn) ];
    elseif formularadio==handles.nonlinearint
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('nonlinear:R-t')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.flamespeed(1,plottingcollumn) handles.markstein(1,plottingcollumn) ];
    elseif formularadio==handles.nonlinearderi
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('nonlinear:U-R')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.flamespeed(1,plottingcollumn) handles.markstein(1,plottingcollumn) ];
    elseif formularadio==handles.KBLforint
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('KBL:R-t')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.flamespeed(1,plottingcollumn) handles.markstein(1,plottingcollumn) ];
    elseif formularadio==handles.KBLforderi
    content = [cellstr(char(handles.filename)) handles.pressure handles.phi cellstr('KBL:U-R')  str2double(get(handles.lowradius,'String')) str2double(get(handles.highradius,'String')) ...
      handles.flamespeed(1,plottingcollumn) handles.markstein(1,plottingcollumn) ];
    end
    

    xlswrite([PathName,FileName],content,1,range);  
    
end


return


% --- Executes on button press in allbutton.
function allbutton_Callback(hObject, eventdata, handles)
% hObject    handle to allbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allbutton

if get(handles.allbutton,'Value') == 0
    
set(handles.lowerradiusbox,'Visible','on');
set(handles.highradiusbox,'Visible','on');
set(handles.lowradius,'Visible','on');
set(handles.highradius,'Visible','on');
set(handles.text28,'Visible','on');
set(handles.text29,'Visible','on');
    
elseif  get(handles.allbutton,'Value') == 1
 
set(handles.lowerradiusbox,'Visible','off');
set(handles.highradiusbox,'Visible','off');
set(handles.lowradius,'Visible','off');
set(handles.highradius,'Visible','off');
set(handles.text28,'Visible','off');
set(handles.text29,'Visible','off');

end
    


calculatestart_Callback(hObject, eventdata, handles);

return


% --- Executes during object creation, after setting all properties.
function calculatestart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calculatestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in givenr0.
function givenr0_Callback(hObject, eventdata, handles)
% hObject    handle to givenr0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of givenr0

if get(handles.givenr0,'Value') == 0
    

set(handles.r0box,'Visible','off');
    
elseif  get(handles.givenr0,'Value') == 1
set(handles.r0box,'Visible','on');

end
    


calculatestart_Callback(hObject, eventdata, handles);


return

function r0box_Callback(hObject, eventdata, handles)
% hObject    handle to r0box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r0box as text
%        str2double(get(hObject,'String')) returns contents of r0box as a double

calculatestart_Callback(hObject, eventdata, handles);


return



% --- Executes during object creation, after setting all properties.
function r0box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r0box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function linecolorbox_Callback(hObject, eventdata, handles)
% hObject    handle to linecolorbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linecolorbox as text
%        str2double(get(hObject,'String')) returns contents of linecolorbox as a double


% --- Executes during object creation, after setting all properties.
function linecolorbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linecolorbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
