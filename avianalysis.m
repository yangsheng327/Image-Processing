function varargout = avianalysis(varargin)
% avianalysis M-file for avianalysis.fig
%      avianalysis, by itself, creates a new avianalysis or raises the existing
%      singleton*.
%
%      H = avianalysis returns the handle to a new avianalysis or the handle to
%      the existing singleton*.
%
%      avianalysis('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in avianalysis.M with the given input arguments.
%
%      avianalysis('Property','Value',...) creates a new avianalysis or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before avianalysis_OpeningFunction gets called.
%      An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to avianalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help avianalysis

% Last Modified by GUIDE v2.5 22-Aug-2013 11:07:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @avianalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @avianalysis_OutputFcn, ...
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


% --- Executes just before avianalysis is made visible.
function avianalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to avianalysis (see VARARGIN)
global issaved
global deletion
deletion=0;
issaved=1;


% Choose default command line output for avianalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes avianalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = avianalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in selectfilebutton.
function selectfilebutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectfilebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;
global issaved
warning off 'MATLAB:aviread:FunctionToBeRemoved';
warning off 'MATLAB:aviinfo:FunctionToBeRemoved';


if issaved
    
ftxt = fopen('defaultfolder.txt');
defaultfolder = fgetl(ftxt);
fclose(ftxt);

    
    
    [handles.filename, handles.dirname] = uigetfile('*.avi','File to be Analyzed',defaultfolder);

delete('defaultfolder.txt');
ftxt=fopen('defaultfolder.txt','w');
fprintf(ftxt, '%s',handles.dirname);
fclose(ftxt);

    
    
    %If the cancel button has not been pushed,
    if (~isequal(handles.filename,0) && ~isequal(handles.dirname,0))

        %Set default values for saved data location
        handles.savefilename=[char(strread(handles.filename,'%[^.].%*s')),'.csv'];
        handles.savedirname=handles.dirname;

        set(handles.figure1,'Units','pixels');
        figuresize=get(handles.figure1,'Position');
        
        set(handles.currentfile,'String',handles.filename)
        info=aviinfo([handles.dirname,handles.filename]);
        cla(handles.axes1)
         set(handles.axes1,...
             'Position',[350 figuresize(4)/2-info.Height/2 info.Width info.Height],...
             'Units','pixels',...
            'Visible','on');
set(handles.currentfile,'Position',[350 100 110 2]);

        set(handles.currentframe,'String','1');

        %Make Analysis controls invisible
        changevisible(hObject,eventdata,handles,'Off');

        handles.movielength=info.NumFrames;
        analyzeframe(hObject, eventdata, handles);
    else
        changevisible(hObject,eventdata,handles,'Off')
    end
else
   set(handles.yesnopanel,'Visible','On')
end


handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

  



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1




% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.sliderval,'String',get(handles.slider1,'Value'))
analyzeframe(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderval_Callback(hObject, eventdata, handles)
% hObject    handle to sliderval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sliderval as text
%        str2double(get(hObject,'String')) returns contents of sliderval as a double


% --- Executes during object creation, after setting all properties.
set(handles.slider1,'Value',str2double(get(handles.sliderval,'String')))
analyzeframe(hObject, eventdata, handles)


function sliderval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end




function analyzeframe(hObject, eventdata, handles)


%clc;
global forfig2 matrixfig3_1 matrixfig3_2 plotmap
global FIX
global xc_mouse
global yc_mouse
global resultsmatrix
global outputmatrix
global issaved

global x_ignleft
global x_ignright
global x_shiftforign


warning off 'MATLAB:aviread:FunctionToBeRemoved';
warning off 'MATLAB:aviinfo:FunctionToBeRemoved';

issaved=0;

if (~strcmp(get(handles.currentframe,'String'),'None'))% If a file has been selected
    
    %clc;
    i=str2double(get(handles.currentframe,'String'));

    if i==1 %if this is the first frame
        clear img
        clear diffmatrix
        clear  compareframe
        clear map

        [img,map]=frame2im(aviread([handles.dirname,handles.filename],i));
        %[img,map]=frame2im(mmreader([handles.dirname,handles.filename],i));
        colormap(map);
        clear map;
        
        hold off;
        image(img,'Parent',handles.axes1);
        
        %Change status and wait for input of center of ignition
        set(handles.status,'String','Please click on the center of ignition')

        [xc_mouse, yc_mouse]=ginput(1);
        
        set(handles.status,'String','Please select left and right boundaries of ignition wire')
        [x_ignleft, ~]=ginput(1);
        [x_ignright, ~]=ginput(1);
        
        
        set(handles.status,'String','Please select the blank location to cover ignition wire area')
        [x_blankforign, ~]=ginput(1);
        x_shiftforign = x_blankforign - xc_mouse;
        
        [img,map]=frame2im(aviread([handles.dirname,handles.filename],i));
        img = removeign(img,x_ignleft,x_ignright,x_shiftforign);
        image(img,'Parent',handles.axes1);
      
        FIX = {};

        resultsmatrix=zeros(handles.movielength,1); %Creates correctly sized resultsmatrix
        
        outputmatrix=zeros(handles.movielength,8);
                    %    R_area
                    %    R_circ
                    %    R_ellicirc
                    %    R_twocirc
                    %    R_vertical
                    %    R_horizontal
                    %    R_leftcirc
                    %    R_rightcirc

        %Determine window size and set default first frame values
        
     %   left=xcenter-5;
     %   right=xcenter+5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        firstframe = img;
        
       
        
        
        clear img;
             
      
    set(handles.markinframe,'String','2');
    set(handles.markoutframe,'String',num2str(handles.movielength));
    

        %Make analysis controls visible
        changevisible(hObject,eventdata,handles,'On');
        
        set(handles.currentframe,'String','2');
                     
%        savedata(hObject,eventdata,handles,i,left,right)
        analyzeframe(hObject, eventdata, handles);
  
  else %This is not the first frame
         
        set(handles.status,'String','Please Analyze Current Frame');
        
        [img,map]=frame2im(aviread([handles.dirname,handles.filename],i));
        img = removeign(img,x_ignleft,x_ignright,x_shiftforign);
        
        plotmap = map;
        clear map;
        hold off;
        image(img,'Parent',handles.axes1);
        hold on;
        forfig1 = img;
        [height width]= size(img);
        angle=get(handles.slider2,'Value');
        threshold=get(handles.slider1,'Value')/100;

  i0 = str2double(get(handles.backgroundframe,'String'));
  if i<i0
      i0 = i-1;
  end
  if i0<=0
     i0=1; 
  end
  set(handles.backgroundframe,'String',num2str(i0));
  

  compareframe=frame2im(aviread([handles.dirname,handles.filename],i0));
  compareframe = removeign(compareframe,x_ignleft,x_ignright,x_shiftforign);
 
  diffmatrix=abs(double(compareframe)-double(img));
  forfig2 = diffmatrix;

  methodradio = get(handles.trackingmethod,'SelectedObject');
  %methodradio
  %handles.method1
  %handles.method2
  
if methodradio==handles.method1

     BW = edge(diffmatrix,'canny',threshold);
     if  ~isempty(FIX) 
       for j=1:length(FIX)
             fixinfo = FIX{j};
             if size(fixinfo)*[1;0]~=0
                 for p = 1:(size(fixinfo)*[1;0])
                   m=round(fixinfo(p,1));
                   n=round(fixinfo(p,2));
                   oneorzero = fixinfo(p,3);
                   BW(n,m) = oneorzero;
                 end
             end 
       end
     end
     BW = bwmorph(BW,'dilate',1);
     BW = bwmorph(BW,'close',10);
     BW = bwmorph(BW,'thin',1);    
     BW = imfill(BW,'holes'); 
     %BW = bwmorph(BW,'erode');BW = bwmorph(BW,'erode');
     %BW = bwmorph(BW,'skel');
     %BW = bwmorph(BW,'branchpoints');
     BW = edge(BW,'canny',threshold);
     BW = bwmorph(BW,'close',1000);
     B = bwboundaries(BW,8,'noholes'); % [B,L,N,A]
     if size(B)*[1;0]  ~= 0  
         %%%%% find the B with maximum size %%%%%%%%
            QQ = size(B);
            Amax = 0;
            qmax = 1;
            for q = 1:QQ(1)
                qvector = B{q};
                [Qarea,~,~] = polycenter(qvector(:,2),qvector(:,1));
                if (Qarea > Amax) 
                    Amax = Qarea;
                    qmax = q;
                end
            end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       rvector = B{qmax};
       x = rvector(:,2);
       y = rvector(:,1);
     end

elseif methodradio==handles.method2
     
    Threshold = threshold*100;
    size_threshold = str2double(get(handles.clustersize,'String'));

    %%%%%%%%%%%%%% from Saha  %%%%%%%%%%%%%%%%%
    I1=imabsdiff(compareframe,img); % I0 is flame image; B0 is background
    I2=I1;
    i_noise = Threshold; % minimum intensiy value from subtracted image to be considered as flame
    I2(find(I1<i_noise))=0;
    % % % % subplot(2,2,2);imshow(I2);
    pconn = size_threshold; % minimum number of connected pixel to be considered as flame
    %pconn
    I3 = bwareaopen(I2, pconn); 
    % subplot(2,2,3);imshow(I3);
    I4 = imfill(I3,'holes'); 
    % subplot(2,2,4);imshow(I4);
    info1=regionprops(I4,'Area');
    if size(info1)*[1;0]  ~= 0 
      %figure();imagesc(img);
      for j=1:length(info1)
        A(j)=info1(j).Area;
      end
      imp_label=find(A==max(A));

     BWlb = bwlabel(I4);
     % CC=bwconncomp(I4);

    BW=BWlb;
    BW(find(BWlb~=imp_label))=0;
    BW(find(BWlb==imp_label))=1;

    info=regionprops(BW,'Area');
    %R=sqrt(info.Area/pi())/93.2;
    % % % % % hold off;
    % % % % % clf;
    % % % % % figure(1);

    I5=bwboundaries(BW);
    % imshow(I0_pure);

    XY=I5{1,1};
   
    x = XY(:,2);
    y = XY(:,1);
    
    end
         
end
        
         if (exist('x')==1 && exist('y')==1)

                    if get(handles.whitedisplay,'Value') == 1
                    handles.theplot = plot(x,y,'w-','LineWidth',2,'Parent',handles.axes1);
                    end
               
                   % calculate area and perimeter
                   [Area,xc_area,yc_area] = polycenter(x,y);
                   Peri = polyperimeter(x,y);
                   
                   cita0 = [0:(pi/1000):(2*pi)]';
                   R_area = (Area/pi)^0.5;
                   
         
                       
                   
                   
                   if get(handles.bluedisplay,'Value') == 1
                   % plot COM and the circle based on area
                   xfit_area = xc_area + R_area*cos(cita0);
                   yfit_area = yc_area + R_area*sin(cita0);
                   plot(xfit_area,yfit_area,'b.','Markersize',1.5,'Parent',handles.axes1);
                   plot(xc_area,yc_area,'o','Markerfacecolor','blue','Parent',handles.axes1);
                   end
                   
                   
                    % find selected data in left and right for flame speeds
                    % measurements
                    r_mouse = ((x-xc_mouse).^2+(y-yc_mouse).^2).^0.5;
                    cita_mouse = atan2(y - yc_mouse,x - xc_mouse);
                    x_right=[]; y_right=[];
                    x_left=[]; y_left=[];
                    for j=1:length(cita_mouse)
                        if abs(cita_mouse(j))<=pi*angle/360
                            x_right = [x_right;x(j)];
                            y_right = [y_right;y(j)];
                        elseif abs(cita_mouse(j)-pi)<=pi*angle/360 || abs(cita_mouse(j)+pi)<=pi*angle/360
                            x_left = [x_left;x(j)];
                            y_left = [y_left;y(j)];      
                        end
                    end
                    x_both = [x_left;x_right];
                    y_both = [y_left;y_right];
                    
            if ~isempty(y_both)
                    
                
                    if get(handles.yellowdisplay,'Value') == 1
                    plot(x_both,y_both,'y.','Parent',handles.axes1);
                    matrixfig3_1 = [x_both,y_both];
                    end
                    
                    
                    
                    %get left and right index
                    [~,leftindex] = min(abs(abs(cita_mouse)-pi));
                    [~,rightindex] = min(abs(abs(cita_mouse)));
                    
                    

                    %get horizontal and vertical two points radius
                    [ymax, ymaxindex] = max(y);
                    [ymin, yminindex] = min(y);
                    R_vertical = (ymax-ymin)/2;
                    
                    
                    [xmax, xmaxindex] = max(x);
                    [xmin, xminindex] = min(x);  
                    R_horizontal = (xmax-xmin)/2;
                    
                    
                     if get(handles.purpledisplay,'Value') == 1
                     leftliney = [0:1:height]';
                     leftlinex =  xmin + leftliney*0;
                     rightliney = [0:1:height]';
                     rightlinex = xmax + rightliney*0;
             plot(leftlinex,leftliney,'-','Markersize',1.5,'Parent',handles.axes1,...
                 'Color',[0.8,0,0.8]);
             plot(rightlinex,rightliney,'-','Markersize',1.5,'Parent',handles.axes1,...
                 'Color',[0.8,0,0.8]);
                     end
                    
                    if get(handles.browndisplay,'Value') == 1

                     bottomlinex =  [0:1:width]';
                     bottomliney =  bottomlinex*0 + ymin;
                     toplinex =  [0:1:width]';
                     topliney =  toplinex*0 + ymax;                     
             plot(bottomlinex,bottomliney,'-','Markersize',1.5,'Parent',handles.axes1,...
                 'Color',[0.6,0.6,0]);                            
             plot(toplinex,topliney,'-','Markersize',1.5,'Parent',handles.axes1,...
                 'Color',[0.6,0.6,0]);   
                    end
                    
                    % perform circular fitting with selected data
                    [xc_both,yc_both,R_circ_both] = circfit(x_both,y_both);
                    if get(handles.greendisplay,'Value') == 1
                    xfit_both = xc_both + R_circ_both*cos(cita0);
                    yfit_both = yc_both + R_circ_both*sin(cita0);
                    plot(xfit_both,yfit_both,'g.','Markersize',1.5,'Parent',handles.axes1);
                    
                    matrixfig3_2 = [xfit_both, yfit_both];
                    end
                    R_circ = R_circ_both;
                    
                    
                    % perform left and right data locally
                    [xc_left,yc_left,R_leftcirc] = circfit(x_left,y_left);
                    [xc_right,yc_right,R_rightcirc] = circfit(x_right,y_right);
                    if get(handles.magentadisplay,'Value') == 1
                    xfit_left = xc_left + R_leftcirc*cos(cita0);
                    yfit_left = yc_left + R_leftcirc*sin(cita0);
                    xfit_right = xc_right + R_rightcirc*cos(cita0);
                    yfit_right = yc_right + R_rightcirc*sin(cita0);
                    plot(xfit_left,yfit_left,'m.','Markersize',1.5,'Parent',handles.axes1);
                    plot(xfit_right,yfit_right,'m.','Markersize',1.5,'Parent',handles.axes1);
                    end

      

                   % perform elliptical fitting with selected data
                    R_both_mouse = ((x_both-xc_mouse).^2+(y_both-yc_mouse).^2).^0.5;
                    cita_both_mouse = atan2(y_both-yc_mouse, x_both-xc_mouse);
                    beta0 = [R_circ_both,R_circ_both];
                    
                    beta = nlinfit(cita_both_mouse, R_both_mouse, @ellipradius, beta0);
                    Cosradius = beta(1);
                    Sinradius = beta(2);
                    R_ellicirc = Sinradius^2/Cosradius;

                   if get(handles.reddisplay,'Value') == 1
                       xfit_elli = xc_both + Cosradius*cos(cita0);
                       yfit_elli = yc_both + Sinradius*sin(cita0);
                       plot(xfit_elli,yfit_elli,'r.','Markersize',1.5,'Parent',handles.axes1);
                       
                       xfit_elli_leftc = x(leftindex)+ R_ellicirc + R_ellicirc*cos(cita0);
                       yfit_elli_leftc = yc_mouse + R_ellicirc*sin(cita0);
                       xfit_elli_rightc = x(rightindex) - R_ellicirc + R_ellicirc*cos(cita0);
                       yfit_elli_rightc = yc_mouse + R_ellicirc*sin(cita0);
                       %plot(xfit_elli_leftc,yfit_elli_leftc,'r.','Markersize',1.5,'Parent',handles.axes1);
                       %plot(xfit_elli_rightc,yfit_elli_rightc,'r.','Markersize',1.5,'Parent',handles.axes1);
                       
                       matrixfig3_2 = [matrixfig3_2, xfit_elli_leftc, yfit_elli_leftc];
                       matrixfig3_2 = [matrixfig3_2, xfit_elli_rightc, yfit_elli_rightc];
                       
                   end
                  
                   
                   
                   
                  % perform two circular fitting with selected data  
                  Ones = [-ones(length(x_left),1);ones(length(x_right),1)];
                  CI = [cita_both_mouse,Ones];
                  gamma0 = [R_circ_both,0,0];
                  gamma = nlinfit(CI, R_both_mouse, @localradius, gamma0);
                  R_twocirc = gamma(1);
                  x_shift = gamma(2); 
                  y_shift = gamma(3);
                  
                  if get(handles.cyandisplay,'Value') == 1
                  xfit_local_left = xc_mouse + R_twocirc*cos(cita0) - x_shift;
                  yfit_local_left = yc_mouse + R_twocirc*sin(cita0) + y_shift;
                  xfit_local_right = xc_mouse + R_twocirc*cos(cita0) + x_shift;
                  yfit_local_right = yc_mouse + R_twocirc*sin(cita0) + y_shift; 
                  xfit_local = [xfit_local_left, xfit_local_right];
                  yfit_local = [yfit_local_left, yfit_local_right];
                  plot(xfit_local,yfit_local,'c.','Markersize',1.5,'Parent',handles.axes1);
                  end

                  
                  
                  
            else
                
                  R_circ = 0;
                  R_area = 0;
                  R_ellicirc = 0;
                  R_twocirc = 0;
                  R_vertical = 0;
                  R_horizontal = 0;
                  R_leftcirc = 0;
                  R_rightcirc = 0;
                  
                  
            end
            
              

     else

                    
                  R_circ = 0;
                  R_area = 0;
                  R_ellicirc = 0;
                  R_twocirc = 0;
                  R_vertical = 0;
                  R_horizontal = 0;
                  R_leftcirc = 0;
                  R_rightcirc = 0;
                    
     end


                    resultsmatrix(i,:) = [2*R_circ];
                    outputmatrix(i,:) = [R_circ, R_area, R_ellicirc, R_twocirc, R_vertical, R_horizontal, R_leftcirc, R_rightcirc];

                    
                    
                    set(handles.diameterbox,'String',num2str(2*R_circ));
                    
                    Fileinfo = interpretfilename(handles.filename);
                    timevalue= 1000*str2double(get(handles.currentframe,'String'))/Fileinfo.framerate; % unit:ms
                    radiusvalue = R_circ/Fileinfo.ppcm;   % cm
                    set(handles.timebox,'String',num2str(timevalue));
                    set(handles.radiusbox,'String',num2str(radiusvalue));
                    
                    
    end
      

        clear img
end



   
   
handles.output = hObject;
guidata(hObject, handles);
    


function R = ellipradius(beta,ci)
    a = beta(1); b = beta(2);
    R  = ((a.*cos(ci)).^2 + (b.*sin(ci)).^2).^0.5;


    

function R = localradius(beta,CI)

    R0 = beta(1); 
    x0 = beta(2);
    y0 = beta(3);
    ci = CI(:,1);
    leftorright = CI(:,2);
    
    x = R0.*cos(ci) + leftorright*x0;
    y = R0.*sin(ci) + y0;
    R = (x.^2+y.^2).^0.5;
    
    
    
function img = removeign(img,x_ignleft,x_ignright,x_shiftforign)
    
        
%dimen = size(img)
%ylimit = dimen(1);

if x_ignleft <= x_ignright

for i=x_ignleft:x_ignright
    img(:,round(i)) = img(:,round(i+x_shiftforign));
end

else
 
    for i=x_ignright:x_ignleft
     img(:,round(i)) = img(:,round(i+x_shiftforign));
    end

    
end

return;


% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global FIX

FIX = {};

i=str2double(get(handles.currentframe,'String'));
i0=str2double(get(handles.backgroundframe,'String'));


if i< handles.movielength
    i=i+1;
end

set(handles.currentframe,'String',num2str(i));

if get(handles.followingbox,'Value') == 1
    if i0<i
    set(handles.backgroundframe,'String',num2str(i0+1));
    end
end

analyzeframe(hObject, eventdata, handles);




% --- Executes on button press in previousbutton.
function previousbutton_Callback(hObject, eventdata, handles)
% hObject    handle to previousbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global FIX

FIX = {};

i=str2double(get(handles.currentframe,'String'));
i0=str2double(get(handles.backgroundframe,'String'));
if i>2
    i=i-1;
else
    if i>1
        i=i-1;
    end
    changevisible(hObject,eventdata,handles,'Off')
end

if get(handles.followingbox,'Value') == 1
    if i0>1
    set(handles.backgroundframe,'String',num2str(i0-1));
    end
end

set(handles.currentframe,'String',num2str(i));


analyzeframe(hObject, eventdata, handles)





function vertrange_Callback(hObject, eventdata, handles)
% hObject    handle to vertrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vertrange as text
%        str2double(get(hObject,'String')) returns contents of vertrange as a double
analyzeframe(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function vertrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vertrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.vertrangebox,'String',num2str(get(handles.slider2,'Value')))
analyzeframe(hObject, eventdata, handles)



% function slider2_CreateFcn(hObject, eventdata, handles)
%     
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end





function changedefaultbutton_Callback(hObject, eventdata, handles)


[handles.savefilename, handles.savedirname] = uiputfile('*.csv','Data Save Location');

handles.output = hObject;
guidata(hObject, handles);



function changevisible(hObject,eventdata,handles,setting)
    
set(handles.text6,'Visible',setting);
set(handles.slider2,'Visible',setting);
set(handles.slider1,'Visible',setting);
set(handles.vertrangebox,'Visible',setting);
set(handles.text2,'Visible',setting);
set(handles.sliderval,'Visible',setting);
set(handles.previousbutton,'Visible',setting);
set(handles.nextbutton,'Visible',setting);
set(handles.currentframename,'Visible',setting);
set(handles.currentframe,'Visible',setting);
set(handles.diameterbox,'Visible',setting);
set(handles.diametername,'Visible',setting);
set(handles.savebutton,'Visible',setting);
set(handles.backgroundframe,'Visible',setting);
set(handles.text16,'Visible',setting);
set(handles.previousbackground,'Visible',setting);
set(handles.nextbackground,'Visible',setting);
set(handles.backgroundframe,'Visible',setting);
set(handles.followingbox,'Visible',setting);
set(handles.yellowdisplay,'Visible',setting);
set(handles.whitedisplay,'Visible',setting);
set(handles.bluedisplay,'Visible',setting);
set(handles.breakbutton,'Visible',setting);
set(handles.connectbutton,'Visible',setting);
set(handles.clearmanual,'Visible',setting);
set(handles.clearallmanual,'Visible',setting);
set(handles.gotobutton,'Visible',setting);
set(handles.gotoframe,'Visible',setting);
set(handles.markinbutton,'Visible',setting);
set(handles.markoutbutton,'Visible',setting);
set(handles.markinframe,'Visible',setting);
set(handles.markoutframe,'Visible',setting);
set(handles.getfigurebutton,'Visible',setting);
set(handles.text21,'Visible',setting);
set(handles.autoframe,'Visible',setting);
set(handles.autorunbutton,'Visible',setting);
set(handles.clearrambutton,'Visible',setting);
set(handles.greendisplay,'Visible',setting);
set(handles.reddisplay,'Visible',setting);
set(handles.cyandisplay,'Visible',setting);
set(handles.magentadisplay,'Visible',setting);
set(handles.purpledisplay,'Visible',setting);
set(handles.browndisplay,'Visible',setting);
set(handles.getonefigure,'Visible',setting);
set(handles.addfigurebutton,'Visible',setting);
set(handles.outputbox,'Visible',setting);
%set(handles.opendifferentiation,'Visible',setting);


set(handles.timebox,'Visible',setting);
set(handles.radiusbox,'Visible',setting);
set(handles.timeboxtext,'Visible',setting);
set(handles.radiusboxtext,'Visible',setting);

set(handles.clustersize,'Visible',setting);
set(handles.clustersize1,'Visible',setting);

handles.output = hObject;
guidata(hObject, handles);


% function changeframe(hObject,eventdata,handles,setting)
% analyzeframe(hObject, eventdata, handles)

function diameterbox_Callback(hObject, eventdata, handles)
% hObject    handle to diameterbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameterbox as text
%        str2double(get(hObject,'String')) returns contents of diameterbox as a double


% --- Executes during object creation, after setting all properties.
function diameterbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diameterbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global resultsmatrix
global issaved

DefaultName = [handles.savedirname,handles.savefilename];
[FileName,PathName,FilterIndex] = uiputfile('*.csv','Save Data File',DefaultName);

if FilterIndex == 0
    
else
    
issaved=1;

ilow = str2double(get(handles.markinframe,'String'));
ihigh = str2double(get(handles.markoutframe,'String'));

frameswrite = (ilow:ihigh)';
matrixwrite = resultsmatrix(ilow:ihigh,:);


csvwrite([PathName,FileName],[frameswrite matrixwrite],1,0);
%csvwrite([handles.savedirname,handles.savefilename],[(1:1:handles.movielength)' 2*(resultsmatrix(:,2)-xcenter)],1,0);

set(handles.status,'String','Results Saved');

end


return




function keypress(hObject, eventdata, handles)
if (~strcmp(get(handles.currentframe,'String'),'None') && ~strcmp(get(handles.currentframe,'String'),'1'))
    switch double(get(handles.figure1,'CurrentCharacter'))
        case 29
            nextbutton_Callback(hObject, eventdata, handles)
        case 13
            nextbutton_Callback(hObject, eventdata, handles)
        case 28
            previousbutton_Callback(hObject, eventdata, handles)
        otherwise
    end       
end



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global deletion
set(handles.yesnopanel,'Visible','Off')
pause(0.1);
savebutton_Callback(hObject, eventdata, handles)
if ~deletion 
    selectfilebutton_Callback(hObject, eventdata, handles)
    outputbox_Callback(hObject, eventdata, handles)
    
else
    delete(handles.figure1);
end




% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global issaved
global deletion
issaved=1;
set(handles.yesnopanel,'Visible','Off')
%Make sure that yesno is not visible (for some reason there is a delay)
pause(0.1);
if ~deletion 
    selectfilebutton_Callback(hObject, eventdata, handles)
else
    delete(handles.figure1);
end






% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global deletion
global issaved

if issaved
    delete(hObject);
else
    deletion=1;
    set(handles.yesnopanel,'Visible','On')
end

        




% --- Executes during object creation, after setting all properties.
function framecomparison_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framecomparison (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function scalelimitset_Callback(hObject, eventdata, handles)
% hObject    handle to scalelimitset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalelimitset as text
%        str2double(get(hObject,'String')) returns contents of scalelimitset as a double


% --- Executes during object creation, after setting all properties.
function scalelimitset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalelimitset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function backgroundframe_Callback(hObject, eventdata, handles)
% hObject    handle to backgroundframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of backgroundframe as text
%        str2double(get(hObject,'String')) returns contents of backgroundframe as a double


% --- Executes during object creation, after setting all properties.
function backgroundframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to backgroundframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previousbackground.
function previousbackground_Callback(hObject, eventdata, handles)
% hObject    handle to previousbackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i0=str2double(get(handles.backgroundframe,'String'));
if i0>1
    i0=i0-1;
end

set(handles.backgroundframe,'String',num2str(i0));
analyzeframe(hObject, eventdata, handles);


% --- Executes on button press in nextbackground.
function nextbackground_Callback(hObject, eventdata, handles)
% hObject    handle to nextbackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=str2double(get(handles.currentframe,'String'));
i0=str2double(get(handles.backgroundframe,'String'));
if i0<i
    i0=i0+1;
end
set(handles.backgroundframe,'String',num2str(i0));
analyzeframe(hObject, eventdata, handles)



% --- Executes on button press in followingbox.
function followingbox_Callback(hObject, eventdata, handles)
% hObject    handle to followingbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in yellowdisplay.
function yellowdisplay_Callback(hObject, eventdata, handles)
% hObject    handle to yellowdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of yellowdisplay

analyzeframe(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function changedefaultbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to changedefaultbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





% --- Executes on button press in clearmanual.
function clearmanual_Callback(hObject, eventdata, handles)
% hObject    handle to clearmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global FIX

 if ~isempty(FIX)
        newFIX = {};
        for j=1:(length(FIX)-1)
            fixmatrix = FIX{j};
            newFIX = [newFIX,fixmatrix];
        end
        FIX = {};
        FIX = newFIX;
        clear newFIX;
 end

analyzeframe(hObject, eventdata, handles);

% handles.output = hObject;
% guidata(hObject, handles);



% --- Executes on button press in whitedisplay.
function whitedisplay_Callback(hObject, eventdata, handles)
% hObject    handle to whitedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of whitedisplay
analyzeframe(hObject, eventdata, handles);

% --- Executes on button press in bluedisplay.
function bluedisplay_Callback(hObject, eventdata, handles)
% hObject    handle to bluedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bluedisplay
analyzeframe(hObject, eventdata, handles);


% --- Executes on button press in breakbutton.
function breakbutton_Callback(hObject, eventdata, handles)
% hObject    handle to breakbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global FIX

rect = getrect(handles.axes1);

rectxmin = round(rect(1));
rectymin = round(rect(2));
rectwidth = round(rect(3));
rectheight = round(rect(4));

 fixmatrix=[];
for m = rectxmin:(rectxmin+rectwidth)
    for n = rectymin:(rectymin+rectheight)
     fixmatrix = [ fixmatrix; m,n,0];  
    end
end
  FIX = [FIX,fixmatrix];

analyzeframe(hObject, eventdata, handles);

% handles.output = hObject;
% guidata(hObject, handles);





% --- Executes on button press in connectbutton.
function connectbutton_Callback(hObject, eventdata, handles)
% hObject    handle to connectbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global FIX


theline = getline(handles.axes1);


xpts = round(theline(:,1));
ypts = round(theline(:,2));


 fixmatrix=[];
for i=1:(length(xpts)-1)
    
     if xpts(i)==xpts(i+1)
         
           if ypts(i)==ypts(i+1)
               m = xpts(i);
               n = ypts(i);
               fixmatrix = [ fixmatrix; m,n,1]; 
           else
               for n = ypts(i):(sign(ypts(i+1)-ypts(i))):ypts(i+1)
                     m = xpts(i);
                     fixmatrix = [ fixmatrix; m,n,1]; 
               end
           end
     else
           for m = xpts(i):(sign(xpts(i+1)-xpts(i))):xpts(i+1)
               n = round(ypts(i) + (m-xpts(i))*(ypts(i+1)-ypts(i))/(xpts(i+1)-xpts(i)));
               fixmatrix = [ fixmatrix; m,n,1]; 
           end
     end
end


FIX = [FIX,fixmatrix];

%fixmatrix
%plot(fixmatrix(:,1),fixmatrix(:,2),'Parent',handles.axes1);

analyzeframe(hObject, eventdata, handles);


%handles.output = hObject;
%guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function connectbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to connectbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in clearallmanual.
function clearallmanual_Callback(hObject, eventdata, handles)
% hObject    handle to clearallmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global FIX

clc;
        FIX = {};

analyzeframe(hObject, eventdata, handles);

% handles.output = hObject;
% guidata(hObject, handles);



% --- Executes on button press in gotobutton.
function gotobutton_Callback(hObject, eventdata, handles)
% hObject    handle to gotobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


gotoi=str2double(get(handles.gotoframe,'String'));

if handles.movielength<gotoi
gotoi = handles.movielength;
elseif gotoi<2
gotoi = 2;   
end

set(handles.currentframe,'String',num2str(gotoi));

analyzeframe(hObject, eventdata, handles);

% handles.output = hObject;
% guidata(hObject, handles);



function gotoframe_Callback(hObject, eventdata, handles)
% hObject    handle to gotoframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gotoframe as text
%        str2double(get(hObject,'String')) returns contents of gotoframe as a double


% --- Executes during object creation, after setting all properties.
function gotoframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gotoframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in markinbutton.
function markinbutton_Callback(hObject, eventdata, handles)
% hObject    handle to markinbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = str2double(get(handles.currentframe,'String'));

set(handles.markinframe,'String',num2str(i));


handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in markoutbutton.
function markoutbutton_Callback(hObject, eventdata, handles)
% hObject    handle to markoutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = str2double(get(handles.currentframe,'String'));

set(handles.markoutframe,'String',num2str(i));



handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in getfigurebutton.
function getfigurebutton_Callback(hObject, eventdata, handles)
% hObject    handle to getfigurebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global forfig2 plotmap matrixfig3_1 matrixfig3_2

N = 1;

info=aviinfo([handles.dirname,handles.filename]);
 width = info.Width;
 height = info.Height;

handles.firstfigure=figure('Units','pixels','Position',[150 50 N*width N*height]);
handles.firstaxis=axes('Parent',handles.firstfigure,'Units','pixels','Position',[0 0 N*width N*height]); 

cla(handles.firstaxis);
copyobj(get(handles.axes1,'Children'),handles.firstaxis);
colormap(plotmap); axis equal;



handles.secondfigure=figure('Units','pixels','Position',[150 50 N*width N*height]);
handles.secondaxis=axes('Parent',handles.secondfigure,'Units','pixels','Position',[0 0 N*width N*height]);
cla(handles.secondaxis);
image(forfig2,'Parent',handles.secondaxis);
colormap(plotmap); axis equal;


 handles.thirdfigure=figure('Units','pixels','Position',[100 100 N*info.Width N*info.Height]);
 handles.thirdaxis=axes('Parent',handles.thirdfigure,'Units','pixels','Position',[0 0 N*width N*height]);
 cla(handles.thirdaxis);
 x_both=matrixfig3_1(:,1);y_both=matrixfig3_1(:,2);
 xfit_circ=matrixfig3_2(:,1); yfit_circ=matrixfig3_2(:,2);
 xfit_elli_circ_left=matrixfig3_2(:,3);yfit_elli_circ_left=matrixfig3_2(:,4);
 xfit_elli_circ_right=matrixfig3_2(:,5);yfit_elli_circ_right=matrixfig3_2(:,6);
 plot(x_both,y_both,'k+','Parent',handles.thirdaxis); hold on;
 plot(xfit_circ,yfit_circ,'b-','Parent',handles.thirdaxis); 
 plot(xfit_elli_circ_left,yfit_elli_circ_left,'r-','Parent',handles.thirdaxis); 
 plot(xfit_elli_circ_right,yfit_elli_circ_right,'r-','Parent',handles.thirdaxis); 
 axis(handles.thirdaxis,[0 width 0 height]); axis equal;


handles.output = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gotobutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gotobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function clearallmanual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clearallmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function clearmanual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clearmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function selectfilebutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectfilebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function followingbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to followingbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function previousbackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to previousbackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function nextbackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nextbackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function vertrangebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vertrangebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function currentframename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentframename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function currentframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function diametername_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diametername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function yellowdisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yellowdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function bluedisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bluedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function whitedisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whitedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function getfigurebutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to getfigurebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function markinbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markinbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function markinframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markinframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function savebutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function markoutbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markoutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function markoutframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markoutframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function yesnopanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yesnopanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function autominbox_Callback(hObject, eventdata, handles)
% hObject    handle to autominbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of autominbox as text
%        str2double(get(hObject,'String')) returns contents of autominbox as a double


% --- Executes during object creation, after setting all properties.
function autominbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to autominbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function autoframe_Callback(hObject, eventdata, handles)
% hObject    handle to autoframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of autoframe as text
%        str2double(get(hObject,'String')) returns contents of autoframe as a double


% --- Executes during object creation, after setting all properties.
function autoframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to autoframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autorunbutton.
function autorunbutton_Callback(hObject, eventdata, handles)
% hObject    handle to autorunbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

autoframe = abs(round(str2double(get(handles.autoframe,'String'))));

if autoframe < 2
    autoframe = 2;
elseif autoframe > handles.movielength
    autoframe = handles.movielength;
end

i = abs(round(str2double(get(handles.currentframe,'String'))));

if autoframe == i
        analyzeframe(hObject, eventdata, handles);
elseif i < autoframe
        set(handles.currentframe,'String',num2str(i));   
        analyzeframe(hObject, eventdata, handles);
        while abs(round(str2double(get(handles.currentframe,'String')))) < autoframe
            nextbutton_Callback(hObject, eventdata, handles);
        end
elseif i > autoframe
        set(handles.currentframe,'String',num2str(i));   
        analyzeframe(hObject, eventdata, handles);
        while abs(round(str2double(get(handles.currentframe,'String')))) > autoframe
            previousbutton_Callback(hObject, eventdata, handles);
        end
end

handles.output = hObject;
guidata(hObject, handles);



% --- Executes on button press in clearrambutton.
function clearrambutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearrambutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global resultsmatrix
%global issaved


dimen = size(resultsmatrix);
leftzeromatrix = zeros(dimen(1),dimen(1));
rightzeromatrix = zeros(dimen(2),dimen(2));
resultsmatrix = leftzeromatrix*resultsmatrix*rightzeromatrix;




dimen = size(outputmatrix);
leftzeromatrix = zeros(dimen(1),dimen(1));
rightzeromatrix = zeros(dimen(2),dimen(2));
outputmatrix = leftzeromatrix*outputmatrix*rightzeromatrix;



 handles.output = hObject;
 guidata(hObject, handles);


% --- Executes on button press in greendisplay.
function greendisplay_Callback(hObject, eventdata, handles)
% hObject    handle to greendisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of greendisplay
analyzeframe(hObject, eventdata, handles);



% --- Executes on button press in reddisplay.
function reddisplay_Callback(hObject, eventdata, handles)
% hObject    handle to reddisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reddisplay


analyzeframe(hObject, eventdata, handles);




% --- Executes during object creation, after setting all properties.
function reddisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reddisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function nextbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in cyandisplay.
function cyandisplay_Callback(hObject, eventdata, handles)
% hObject    handle to cyandisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cyandisplay


analyzeframe(hObject, eventdata, handles);


% --- Executes on button press in magentadisplay.
function magentadisplay_Callback(hObject, eventdata, handles)
% hObject    handle to magentadisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of magentadisplay


analyzeframe(hObject, eventdata, handles);





% --- Executes on button press in purpledisplay.
function purpledisplay_Callback(hObject, eventdata, handles)
% hObject    handle to purpledisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of purpledisplay

analyzeframe(hObject, eventdata, handles);


% --- Executes on button press in getonefigure.
function getonefigure_Callback(hObject, eventdata, handles)
% hObject    handle to getonefigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N = 1;
info=aviinfo([handles.dirname,handles.filename]);
width = info.Width;
height = info.Height;

handles.thefigure=figure('Units','pixels','Position',[100 100 N*info.Width N*info.Height]);
handles.theaxis=axes('Parent',handles.thefigure,'Units','pixels','Position',[0 0 N*width N*height]);
cla(handles.theaxis);
axis(handles.theaxis,[0 width 0 height]); axis equal;


handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in addfigurebutton.
function addfigurebutton_Callback(hObject, eventdata, handles)
% hObject    handle to addfigurebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


copiedplot = copyobj(handles.theplot,handles.theaxis);
set(copiedplot,'Color','k');


handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in browndisplay.
function browndisplay_Callback(hObject, eventdata, handles)
% hObject    handle to browndisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of browndisplay


analyzeframe(hObject, eventdata, handles);




% --- Executes on button press in outputbox.
function outputbox_Callback(hObject, eventdata, handles)
% hObject    handle to outputbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global outputmatrix


DefaultName = [handles.savedirname,handles.filename(1:(length(handles.filename)-4)),'.xls'];

[FileName,PathName,FilterIndex] = uiputfile({'*.xls';'*.csv'},'Output File',DefaultName);


if FilterIndex == 0
    
else

Fileinfo = interpretfilename(handles.filename);

ilow = str2double(get(handles.markinframe,'String'));
ihigh = str2double(get(handles.markoutframe,'String'));

timewrite = ((ilow:ihigh)')/Fileinfo.framerate;
radiuswrite = outputmatrix(ilow:ihigh,:)/Fileinfo.ppcm;

   nameunit = {'Time','R_circ','R_area','R_ellicirc','R_twocirc','R_vertical','R_horizontal','R_leftcirc','R_rightcirc';
               'sec',  'cm',    'cm',      'cm',           'cm',      'cm',         'cm',        'cm',       'cm'};
   finalwritematrix = [timewrite,radiuswrite];
   xlswrite([PathName,FileName],nameunit,1,'A1');
   xlswrite([PathName,FileName],finalwritematrix,1,'A3');

end


return









% --- Executes on button press in opendifferentiation.
function opendifferentiation_Callback(hObject, eventdata, handles)
% hObject    handle to opendifferentiation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Derivative;


return


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.yesnopanel,'Visible','on');
% 
% 
% return


% --- Executes during object deletion, before destroying properties.
function yesnopanel_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to yesnopanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function clustersize_Callback(hObject, eventdata, handles)
% hObject    handle to clustersize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clustersize1 as text
%        str2double(get(hObject,'String')) returns contents of clustersize1 as a double
analyzeframe(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function clustersize1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustersize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in trackingmethod.
function trackingmethod_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in trackingmethod 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

analyzeframe(hObject, eventdata, handles);




% --- Executes during object creation, after setting all properties.
function clustersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
