function varargout = flatnessPlotter(varargin)
% FlATNESSPLOTTER MATLAB code for flatnessPlotter.fig
%      FLATNESSPLOTTER, by itself, creates a new FLATNESSPLOTTER or raises the existing
%      singleton*.
%
%      H = FLATNESSPLOTTER returns the handle to a new FLATNESSPLOTTER or the handle to
%      the existing singleton*.
%
%      FLATNESSPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLATNESSPLOTTER.M with the given input arguments.
%
%      FLATNESSPLOTTER('Property','Value',...) creates a new FLATNESSPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flatnessPlotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flatnessPlotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
%     Copyright (C) 2015  Devin C Prescott
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     
%     Author:
%     Devin C Prescott
%     devin.c.prescott@gmail.com

% Edit the above text to modify the response to help flatnessPlotter

% Last Modified by GUIDE v2.5 12-Jun-2014 12:19:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flatnessPlotter_OpeningFcn, ...
                   'gui_OutputFcn',  @flatnessPlotter_OutputFcn, ...
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


% --- Executes just before flatnessPlotter is made visible.
function flatnessPlotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flatnessPlotter (see VARARGIN)

% Choose default command line output for flatnessPlotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
handles.drawables = {};
handles.KeepFigs = findall(0,'type','figure');
% UIWAIT makes flatnessPlotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = flatnessPlotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% UI Open The Data File
[data,FileName,PathName]=CmmFileParser();
handles.FileName = FileName;
handles.PathName = PathName;
handles.data = data;
guidata(hObject, handles);

%% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get Variables
try
FileName = handles.FileName;
catch
    msgbox('You need to open a file first');
    return
end
KeepFigs = handles.KeepFigs;
DeleteFigs = findall(0,'type','figure');
PathName = handles.PathName;
drawables = handles.drawables;
data = handles.data;
FlatTol = str2double(get(handles.edit_flattol,'String'));
SD = str2double(get(handles.edit_sd,'String'));

% Close Any New Figures 
close(setdiff(DeleteFigs,KeepFigs))

x = data(:,1);
y = data(:,2);
z = data(:,3);

%% Outlier Z data trimming
% radialmag = 1;
index = true(1,length(x));
z_trim = sortrows(z);
z_trim(1:10)=[];
z_trim(end-11:end)=[];
lower = mean(z_trim)-SD*std(z_trim);
upper = mean(z_trim)+SD*std(z_trim);
for i = 1:length(x)
%     mag = sqrt((x-x(i)).^2+(y-y(i)).^2);
%     
%     while length(z(mag<radialmag))<120
%         radialmag = radialmag+1;
%     end
%     
%     upper = mean(z(mag < radialmag))+SD*std(z(mag < radialmag));
%     lower = mean(z(mag < radialmag))-SD*std(z(mag < radialmag));
%     disp({lower,'<',z(i),'<',upper})
    if z(i) > upper || z(i) < lower
        index(i) = false;
    end
%     radialmag=1;
end
% fprintf('Lower:%1.3f,\t Upper:%1.3f',lower,upper);
x=x(index);
y=y(index);
z=z(index);
z_org = z; % Save original Z data

%% Best Fit Plane - Needed for Flatness Measurement
% Drop NaN elements
x = x(isfinite(x));
y = y(isfinite(x));
z = z(isfinite(x));
xfit = x;
yfit = y;
data = [x,y,z];
% Mean Center the data around zero
data_0 = bsxfun(@minus,data,mean(data));
% Augment the matrix
data_0(:,4) = ones(length(data_0),1);
% Singular Value Decomposition
[~,S,V] = svd(data_0,0);
% Error checking svd
if 100*S(4,4)>S(1,1) % Check last diag elem is >> than first
    msgbox('Data is very non-planar. Data will not be reoriented against a best fit plane.');
else % Best fit plane worked
    % Equation of the best fit plane
    P = V(:,4);
    % Residuals
    zfit = ((P(1).*x + P(2).*y + P(4))/-P(3) + mean(data(:,3),1));
    zRes = z-zfit;

    % Z data is now distance from plane
    z = zRes;
end

% Detrend
x = x-mean(x);
y = y-mean(y);
z = z-mean(z);

%% Plot Probe Points
if get(handles.checkbox_points,'Value')
    figure;
    plot(x,y,'.');
    axis equal
    title(strcat('Probe Points-',FileName),'interpreter','none');
end
% Add drawn shapes
for i = 1:length(drawables)
    hold on
    plot(drawables{i,1},drawables{i,2});
    hold off
end
%% Plot Histogram
if get(handles.checkbox_hist,'Value')
    figure;
    hist(z) 
    hold on
    lsl = mean(z)-FlatTol/2;
    usl = mean(z)+FlatTol/2;
    plot([mean(z) mean(z)],[get(gca,'Ylim')],'--r')
    plot([usl usl],[get(gca,'Ylim')],'--b')
    plot([lsl lsl],[get(gca,'Ylim')],'--g')
    legend('Z Histogram','Z Mean','USL Flatness','LSL Flatness');
    title(strcat('Histogram-',FileName),'interpreter','none');
    xlabel('Residual Height');
    ylabel('Occurances');
    hold off
end

%% Delaunay triangulation - Surface creation
primary = [x,y,z];
if ~isempty(drawables)
secondary{length(drawables)}=[];
    for i = 1:length(drawables)
        secondary{i} = [drawables{i,1}',drawables{i,2}'];
    end
    [DT,xyz,~] = delaunayConstrained(primary,secondary);
else
    DT = delaunayTriangulation([x,y]);
    xyz = [x,y,z];
end
figure
SurfPlot = trisurf(DT.ConnectivityList,xyz(:,1),xyz(:,2),xyz(:,3));
% SurfPlot = get(SurfPlot, 'Parent');

if get(handles.checkbox_points,'Value')
    hold on
    plot3(x,y,z,'.k','MarkerSize',12);
    hold off
end

title(strcat('Surface Plot-',FileName),'interpreter','none');
xlabel('X Axis')
ylabel('Y Axis')
zlim(5*[min(z) max(z)]);
axis square
axis vis3d
colormap(jet)
% l = light('Position',[25 -10 1],'Style', 'local');
lighting phong
shading interp
colorbar East

if get(handles.checkbox_minmax,'Value')
    max_ind = find(z==max(z));
    min_ind = find(z==min(z));
    xmx = x(max_ind);
    ymx = y(max_ind);
    zmx = z(max_ind);
    xmn = x(min_ind);
    ymn = y(min_ind);
    zmn = z(min_ind);
    hold on
    offset = 1.25;
    plot3(xmx,ymx,zmx,'ok','LineWidth',2)
    plot3(xmn,ymn,zmn,'ok','LineWidth',2)
    plot3([xmx,offset*xmx],[ymx,offset*ymx],[zmx,offset*zmx],...
        '-k','LineWidth',2)
    plot3([xmn,offset*xmn],[ymn,offset*ymn],[zmn,offset*zmn],...
        '-k','LineWidth',2)
    text(offset*xmx,offset*ymx,offset*zmx,num2str(max(z),3),...
        'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0],...
        'FontWeight','bold','FontSize',8);
    text(offset*xmn,offset*ymn,offset*zmn,num2str(min(z),3),...
        'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0],...
        'FontWeight','bold','FontSize',8);
    hold off
end
%% Best Fit Plot
if get(handles.checkbox_best,'Value')
    figure;
    n = length(x);
    skip = 1;
    if n > 500
        skip = 2;
    elseif n > 1000
        skip = 5;
    elseif n > 10000
        skip = 10;
    end
    for i = 1:skip:n
        hold on
        plot3([x(i),x(i)],[y(i),y(i)],[zfit(i),z_org(i)],'-k');
    end
    trisurf(tri,xfit,yfit,zfit,ones(length(xfit),1),'EdgeColor','none');
    title(strcat('Best Fit Residuals-',FileName),'interpreter','none');
    xlabel('X Axis')
    ylabel('Y Axis')
    zlim(5*[min(z) max(z)]);
    axis square
    view(3);
    grid on
    hold off
end
%% Dump Variables to Workspace
assignin('base','x',x)
assignin('base','y',y)
assignin('base','z',z)

%% Print to 3D PDF
% PDFname = strcat(PathName,datestr(clock,'ddmmmyy-HHMMSS')); 
% if get(handles.checkbox_pdf,'Value')
%     if exist('fig2u3d','dir')==7
%         try
%             fig2pdf3d(SurfPlot, PDFname, 'media9', 'xelatex');
%         catch
%             disp('You need to run Initiaize_FP.m first');
%         end
%     else
%         disp('You may need to update your path.');
%         disp('Use addpath(genpath(cd));savepath()');
%     end
% end

% --- Executes on button press in checkbox_points.
function checkbox_points_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in checkbox_hist.
function checkbox_hist_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in checkbox_pdf.
function checkbox_pdf_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_pdf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit_rad_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rad as text
%        str2double(get(hObject,'String')) returns contents of edit_rad as a double


% --- Executes during object creation, after setting all properties.
function edit_rad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_center_Callback(hObject, eventdata, handles)
% hObject    handle to edit_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_center as text
%        str2double(get(hObject,'String')) returns contents of edit_center as a double


% --- Executes during object creation, after setting all properties.
function edit_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in pushbutton_circle.
function pushbutton_circle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
drawables = handles.drawables;
center = str2num(get(handles.edit_center,'String'));
radius = str2num(get(handles.edit_rad,'String'));
sgn = sign(radius);
radius = abs(radius);
[x,y]=pol2cart(0:sgn*pi/100:sgn*2*pi,radius);
x(end) = [];
y(end) = [];
hold on
plot(x+center(1),y+center(2));
hold off
axis equal;
drawables{end+1,1} = x+center(1);
drawables{end,2} = y+center(2);
handles.drawables = drawables;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_poly.
function pushbutton_poly_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
drawables = handles.drawables;
PolyPoints = get(handles.uitable_poly,'data');
x = [PolyPoints{:,1},PolyPoints{end,1},PolyPoints{1,1}];
y = [PolyPoints{:,2},PolyPoints{end,2},PolyPoints{1,2}];
hold on
plot(x,y);
hold off
drawables{end+1,1} = x;
drawables{end,2} = y;
handles.drawables = drawables;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uitable_poly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in uitable_poly.
function uitable_poly_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_poly (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

function edit_flattol_Callback(hObject, eventdata, handles)
% hObject    handle to edit_flattol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_flattol as text
%        str2double(get(hObject,'String')) returns contents of edit_flattol as a double

% --- Executes during object creation, after setting all properties.
function edit_flattol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_flattol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_minmax.
function checkbox_minmax_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_minmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_minmax


% --- Executes on button press in checkbox_best.
function checkbox_best_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_best (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_best


% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcbf)
clear
Flatness_Plotter()
return



function edit_sd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sd as text
%        str2double(get(hObject,'String')) returns contents of edit_sd as a double


% --- Executes during object creation, after setting all properties.
function edit_sd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
