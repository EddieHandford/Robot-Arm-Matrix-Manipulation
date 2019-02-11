function varargout = RobotGUI(varargin)
% ROBOTGUI MATLAB code for RobotGUI.fig
%      ROBOTGUI, by itself, creates a new ROBOTGUI or raises the existing
%      singleton*.
%
%      H = ROBOTGUI returns the handle to a new ROBOTGUI or the handle to
%      the existing singleton*.
%
%      ROBOTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROBOTGUI.M with the given input arguments.
%
%      ROBOTGUI('Property','Value',...) creates a  new ROBOTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RobotGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RobotGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RobotGUI

% Last Modified by GUIDE v2.5 09-Oct-2013 15:34:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RobotGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RobotGUI_OutputFcn, ...
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

% --- Executes just before RobotGUI is made visible.
function RobotGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RobotGUI (see VARARGIN)

% Choose default command line output for RobotGUI
handles.output = hObject;
handles.thetas = zeros(1,6);
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
guidata(hObject, handles);

% UIWAIT makes RobotGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = RobotGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes on slider movement.
function angle1_Callback(hObject, eventdata, handles)
% hObject    handle to angle1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.thetas(1) = get(hObject,'Value');
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
IK_solution = IK(endEffPos, R(:,:,7));
[~, XYZ2, R2] = FK(IK_solution(1:6,handles.IK_select));
drawRobot(XYZ2,R2, handles.axes2);
upDateTables(handles, IK_solution, handles.thetas, endEffPos, XYZ(6,:), R(:,:,7));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function angle1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on slider movement.
function angle2_Callback(hObject, eventdata, handles)
% hObject    handle to angle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.thetas(2) = get(hObject,'Value');
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
IK_solution = IK(endEffPos, R(:,:,7));
[~, XYZ2, R2] = FK(IK_solution(1:6,handles.IK_select));
drawRobot(XYZ2,R2, handles.axes2);
upDateTables(handles, IK_solution, handles.thetas, endEffPos, XYZ(6,:), R(:,:,7));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function angle2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on slider movement.
function angle3_Callback(hObject, eventdata, handles)
% hObject    handle to angle3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.thetas(3) = get(hObject,'Value');
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
IK_solution = IK(endEffPos, R(:,:,7));
[~, XYZ2, R2] = FK(IK_solution(1:6,handles.IK_select));
drawRobot(XYZ2,R2, handles.axes2);
upDateTables(handles, IK_solution, handles.thetas, endEffPos, XYZ(6,:), R(:,:,7));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function angle3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on slider movement.
function angle4_Callback(hObject, eventdata, handles)
% hObject    handle to angle4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.thetas(4) = get(hObject,'Value');
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
IK_solution = IK(endEffPos, R(:,:,7));
[~, XYZ2, R2] = FK(IK_solution(1:6,handles.IK_select));
drawRobot(XYZ2,R2, handles.axes2);
upDateTables(handles, IK_solution, handles.thetas, endEffPos, XYZ(6,:), R(:,:,7));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function angle4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on slider movement.
function angle5_Callback(hObject, eventdata, handles)
% hObject    handle to angle5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.thetas(5) = get(hObject,'Value');
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
IK_solution = IK(endEffPos, R(:,:,7));
[~, XYZ2, R2] = FK(IK_solution(1:6,handles.IK_select));
drawRobot(XYZ2,R2, handles.axes2);
upDateTables(handles, IK_solution, handles.thetas, endEffPos, XYZ(6,:), R(:,:,7));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function angle5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on slider movement.
function angle6_Callback(hObject, eventdata, handles)
% hObject    handle to angle6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.thetas(6) = get(hObject,'Value');
[endEffPos, XYZ, R] = FK(handles.thetas);
drawRobot(XYZ,R, handles.axes1);
IK_solution = IK(endEffPos, R(:,:,7));
[~, XYZ2, R2] = FK(IK_solution(1:6,handles.IK_select));
drawRobot(XYZ2,R2, handles.axes2);
upDateTables(handles, IK_solution, handles.thetas, endEffPos, XYZ(6,:), R(:,:,7));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function angle6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on selection change in IK_selector.
function IK_selector_Callback(hObject, eventdata, handles)
% hObject    handle to IK_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IK_selector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IK_selector

handles.IK_select = get(hObject,'Value');
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function IK_selector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IK_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.IK_select = 1;
guidata(hObject, handles);
function drawRobot(xyzCoor,R, axesHandles)

%vector length
vl = 0.15;
xplot = 3;
yplot = 3;
zplot = 3;

%Base frame axis direction
a(1,1,:) = R(:,:,1)*[vl 0 0]';
a(1,2,:) = R(:,:,1)*[0 vl 0]';
a(1,3,:) = R(:,:,1)*[0 0 vl]';

%Axis direction frame 1
a(2,1,:) = R(:,:,2)*[vl 0 0]';
a(2,2,:) = R(:,:,2)*[0 vl 0]';
a(2,3,:) = R(:,:,2)*[0 0 vl]';

%Axis direction frame 2
a(3,1,:) = R(:,:,3)*[vl 0 0]';
a(3,2,:) = R(:,:,3)*[0 vl 0]';
a(3,3,:) = R(:,:,3)*[0 0 vl]';

%Axis direction frame 3
a(4,1,:) = R(:,:,4)*[vl 0 0]';
a(4,2,:) = R(:,:,4)*[0 vl 0]';
a(4,3,:) = R(:,:,4)*[0 0 vl]';

%Axis direction frame 4
a(5,1,:) = R(:,:,5)*[vl 0 0]';
a(5,2,:) = R(:,:,5)*[0 vl 0]';
a(5,3,:) = R(:,:,5)*[0 0 vl]';

%Axis direction frame 5
a(6,1,:) = R(:,:,6)*[vl 0 0]';
a(6,2,:) = R(:,:,6)*[0 vl 0]';
a(6,3,:) = R(:,:,6)*[0 0 vl]';

%Axis direction frame 6
a(7,1,:) = R(:,:,7)*[vl 0 0]';
a(7,2,:) = R(:,:,7)*[0 vl 0]';
a(7,3,:) = R(:,:,7)*[0 0 vl]';


axes(axesHandles);
hold on;
cla;

plot3(xyzCoor(:,1),xyzCoor(:,2),xyzCoor(:,3), 'r*');

%plot lines between the axis origins to represent the robot's body
plot3(xyzCoor(:,1),xyzCoor(:,2),xyzCoor(:,3), 'k');

for i=1:7
    plot3([xyzCoor(i,1), xyzCoor(i,1)+ a(i,1,1)], [xyzCoor(i,2), xyzCoor(i,2)+a(i,1,2)], [xyzCoor(i,3),  xyzCoor(i,3)+a(i,1,3)], 'b');
    plot3([xyzCoor(i,1), xyzCoor(i,1)+ a(i,2,1)], [xyzCoor(i,2), xyzCoor(i,2)+a(i,2,2)], [xyzCoor(i,3),  xyzCoor(i,3)+a(i,2,3)], 'r');
    plot3([xyzCoor(i,1), xyzCoor(i,1)+ a(i,3,1)], [xyzCoor(i,2), xyzCoor(i,2)+a(i,3,2)], [xyzCoor(i,3),  xyzCoor(i,3)+a(i,3,3)], 'g');
end

axis([-1*xplot 1*xplot -1*yplot 1*yplot -1*zplot 1*zplot]) 
function upDateTables(handles, IK, thetas, o6, o5, R06)

set(handles.angleTable, 'data', thetas');
set(handles.endEfftable, 'data', o6);
set(handles.wristCtable, 'data', o5');
set(handles.endEffRottable, 'data', R06);

set(handles.IKtable1, 'data', IK(:,1));
set(handles.IKtable2, 'data', IK(:,2));
set(handles.IKtable3, 'data', IK(:,3));
set(handles.IKtable4, 'data', IK(:,4));


%EDIT THE FOLLOWING TWO FUNCTIONS ONLY
function [effpos, xyzCoor, R] = FK(angles)
 
%Robot parameters (lengths in m).

a = [.6 1.4 0.065 0 0 0];
d = [1.100 0 0 1.600 0 0.372];
alphas = [pi/2, 0 , pi/2 , -pi/2, pi/2 , 0];
%angles(2) = angles(2) + pi/2; 
effpos = [0, 0, 0.123]';  % ( is this meant to be all zeros, this is the source of the difference errror between mine and theos ) 


%Generic T matrix
for i = 1:6
    T(:,:,i) = [cos(angles(i)), -sin(angles(i))*cos(alphas(i)), sin(angles(i))*sin(alphas(i)), a(i)*cos(angles(i));
        sin(angles(i)), cos(angles(i))*cos(alphas(i)),  -(cos(angles(i))*sin(alphas(i))) , a(i)*sin(angles(i));
        0 ,  sin(alphas(i)) , cos(alphas(i)) , d(i);
        0 0 0 1]; 
end
%transfering T01 , T12, T23 .... T56  to T01, T02, T03 ... T06
T0(:, :, 1) = T(:, : , 1);
T0(:, :, 2) = T0(:, :, 1)*T(:, : , 2);
T0(:, :, 3) = T0(:, :, 2)*T(:, : , 3);
T0(:, :, 4) = T0(:, :, 3)*T(:, : , 4);
T0(:, :, 5) = T0(:, :, 4)*T(:, : , 5);
T0(:, :, 6) = T0(:, :, 5)*T(:, : , 6);


%spitting out just the rotaion matrices now, defining R01 as an identity
R(:,:,1) = eye(3);
for i = 2:7
    R(:, :, i) = T0(1:3, 1:3, i-1);
end

xyzCoor = zeros(7,3);

for i = 2:7
xyzCoor(i, :) = T0( 1:3 , 4, i-1 );
end

   effpos = xyzCoor(7,:) ;
   

   
   %so far the forward kinematics seem to work
 
    
    


function solution_IK = IK(effPos, R06)
%A's, D's and alphas
a = [.6 1.4 0.065 0 0 0];
d = [1.100 0 0 1.600 0 0.372];
alphas = [pi/2, 0 , pi/2 , -pi/2, pi/2 , 0];
angles = zeros(6);
WC = zeros(3);
solution_IK = zeros(9,4);
Theta = zeros(4,9);
Basepos = zeros(3,1);



WC = effPos' - d(6)*R06(:,3);


JointOneForward2D = [ a(1) , 0 , d(1) ]'  ;
JointOneBackward2D = [ -a(1), 0 , d(1)]'  ;

WristCenterForward2D = [(sqrt(WC(1)^2 + WC(2)^2)), 0 , WC(3)]'   ;
WristCenterBackward2D = [(sqrt(WC(1)^2 + WC(2)^2)), 0 , WC(3)]'  ;

Theta(1:2,1) = atan2( WC(2),WC(1) );
Theta(3:4,1) = pi + atan2( WC(2),WC(1) );
Theta1 = Theta(1,1);
Theta2 = Theta(3,1);
solution_IK(1 , 1:2) = Theta1 ;
solution_IK(1 , 3:4) = Theta2 ;

%Up to here spits out Theta 1 , both versions for forward and back. Dont
%touch ! Okay maybe you can change the + pi, to be - pi. I dunno ?!


 %was fukin doing this wrong !!!! fucking fixed it now
JointOnetoWCForward = sqrt((WristCenterForward2D(1) - JointOneForward2D(1))^2 + (WristCenterForward2D(3) - JointOneForward2D(3))^2) 
JointOnetoWCBackward = sqrt((WristCenterBackward2D(1) - JointOneBackward2D(1))^2 + (WristCenterBackward2D(3) - JointOneBackward2D(3))^2) 

Link3Length = 1.601319768; %from sqrt of 1.6^2 + 0.065^2
Link2Length = 1.4; %a(2)

if JointOnetoWCForward >= (1.601319768 + 1.4);
  
   disp('ERROR, YOU CANT REACH THAT, LIKE, EVER, WHAT WERE YOU THINKING!'); 
end

if JointOnetoWCBackward >= (1.601319768 + 1.4);
   
   disp('ERROR, YOU CANT REACH THAT FLOPPING YOUR HAND BACKWARDS, YOU AINT GOT ARMS LONG ENOUGH!'); 
end

%Up to here spits out theta 1 , and the two WC - joint 1 lenghts, and some
%aggressive error messages ! 



joint1posforward = [-sin(Theta1)*0.6 , (cos(Theta1)*0.6) , 1.1]
joint1posbackward = [-sin(Theta2)*0.6 , (cos(Theta2)*0.6) , 1.1]

%Up to here, ive now decided to express the x,y,z coordinates of joint 1 in
%true 3d.

s = (WristCenterForward2D(3) - JointOneForward2D(3));
s2 = (WristCenterBackward2D(3) - JointOneBackward2D(3));
r = (WristCenterForward2D(1) - JointOneForward2D(1));
r2 = (WristCenterBackward2D(1) - JointOneBackward2D(1));
%ZZ = WC - joint1pos ; 
%zz = (ZZ(1)) ;





%if sign(Theta2) == 1;
 %   Theta2 = Theta2 - 3.141 ;
%end




%if WC(1) <= 0  ;
    
 %   if sign(WC(1) - joint1pos(1))* sign(sin(Theta1)) == 1 ;
  %      solution_IK(1 , 1:2) = Theta2 ;
  %      solution_IK(1 , 3:4) = Theta1 ;
  %  else
  %      solution_IK(1 , 1:2) = Theta1
  %      solution_IK(1 , 3:4) = Theta2  
  %  end

 %   else
  %  if sign(WC(2) - joint1pos(2))* sign(cos(Theta1)) == 1 ;
  %      solution_IK(1 , 1:2) = Theta2 ;
  %      solution_IK(1 , 3:4) = Theta1 ;
  %  else
   %     solution_IK(1 , 1:2) = Theta1 ;
   %     solution_IK(1 , 3:4) = Theta2 ;
   % end 
%end
    
%Theta1 = forward


%floating angle 1 = inverse tan of (p4z - p1z / p4x - p1x)



%floating angle 1 backwards  = inverse tan of p4z - p1z
%backwards / p4x - plxbackwards
%should spit out two angles

%now we have angle of
%                      




%We can now apply the cosine rule, joy !
%the cosine rule , well , needs a triangle. unfurtunelayt i dont have one,
%so this is the point in which i need to work out the lenghts of the sides
%of the triangle, imma do this on paper.


%Cosine rule

%Forward elbow up
%joint1posforward


a = JointOnetoWCForward 
aa = JointOnetoWCBackward 
b = Link2Length ;
c = Link3Length ;
 % a2 = link3Length
   % a1 = Link2Length
D = ((a^2 - b^2 - c^2)/(2*b*c)) 
DD = ((aa^2 - b^2 - c^2)/(2*b*c)) 


if JointOnetoWCForward < (1.601319768 + 1.4);
    elbowup = -atan2 ( sqrt(1-D^2) , D)  ;
    elbowdown = -atan2 ( -sqrt(1-D^2) , D)  ;
    Theta3EU = elbowup - 1.6112485644548 -pi;
    Theta3ED = elbowdown  - 1.6112485644548 -pi;
else
    elbowup = 0 ;
    elbowdown = 0 ;
    Theta3EU = 0 ;
    Theta3ED = 0 ;
end

if JointOnetoWCBackward < (1.601319768 + 1.4);
    backelbowup = -atan2 ( sqrt(1-DD^2) , DD)  ;
    backelbowdown = -atan2 ( -sqrt(1-DD^2) , DD)  ;
    Theta3BEU = backelbowup -  1.6112485644548 -pi ;
    Theta3BED =  backelbowdown  -  1.6112485644548 -pi ;
else
    backelbowup = 0 ;
    backelbowdown = 0 ;
    Theta3BEU = 0 ;
    Theta3BED = 0 ;
end


solution_IK(3 , 1) = Theta3EU ;
solution_IK(3 , 2) = Theta3ED + (2*3.141);
solution_IK(3 , 3) = Theta3BEU + (2*3.141) ;
solution_IK(3 , 4) =  Theta3BED  + (2*3.141);




sintheta23EU  =  atan2(s,r) - atan2( Link3Length * sin(elbowup) , Link2Length + Link3Length*cos(elbowup))
sintheta23ED  = atan2(s,r) - atan2( Link3Length * sin(elbowdown) , Link2Length + Link3Length*cos(elbowdown))
sintheta23BEU = atan2(s,r) - atan2( Link3Length * sin(backelbowup) , Link2Length + Link3Length*cos(backelbowup))
sintheta23BED = atan2(s,r) - atan2( Link3Length * sin(backelbowdown) , Link2Length + Link3Length*cos(backelbowdown))


ANGLE2 = sintheta23EU 
ANGLE2B = sintheta23ED 
ANGLE2C = sintheta23BEU ;
ANGLE2D = sintheta23BED  ;

solution_IK(2 , 1) = ANGLE2 
solution_IK(2 , 2) = ANGLE2B
solution_IK(2 , 3)  = ANGLE2C ;
solution_IK(2 , 4) = ANGLE2D ;




%at this point in the code, angle 3 works,  angle 2 is just random shit


%Joint3 = WC - d(4)*R02




%calculating where the 1st joint would be with the base rotated
% for i = 1
% T1(:,:,1) = [cos(Theta1), -sin(Theta1)*cos(alphas(i)), sin(Theta1)*sin(alphas(i)), a(i)*cos(Theta1);
   %     sin(Theta1), cos(Theta1)*cos(alphas(i)),  -(cos(Theta1)*sin(alphas(i))) , a(i)*sin(Theta1);
  %      0 ,  sin(alphas(i)) , cos(alphas(i)) , d(i);
 %       0 0 0 1];
%end

%calculating where the 1st joint would be with the base rotated (using the
%second possible soltion, named Theta2 just to further confuse you !
%for i = 0
 %T2(:,:,2) = [cos(Theta2), -sin(Theta2)*cos(alphas(i)), sin(Theta2)*sin(alphas(i)), a(i)*cos(Theta2);
  %      sin(Theta2), cos(Theta2)*cos(alphas(i)),  -(cos(Theta2)*sin(alphas(i))) , a(i)*sin(Theta2);
   %     0 ,  sin(alphas(i)) , cos(alphas(i)) , d(i);
    %    0 0 0 1];
%end


%Dit1 = WC - T1(1:3,4);
% Dit2 = WC - T1(1:3,4);



%%%%%% ONLY DONE FORWARDS UP AND BACKWWARDS UP , REDO FOR ELBOW DOWNS BY
%%%%%% REPRODUCING THE R1 MATRICES ETC !!!!! SHIT

R01   = rMatrix ( alphas(1) , solution_IK(1 , 1)) 
R01ED = rMatrix ( alphas(1) , solution_IK(1 , 2))
R01B  = rMatrix ( alphas(1) , solution_IK(1 , 3)) 
R01BED= rMatrix (alphas(1) , solution_IK(1 , 4))

%LIKE THIS

R12   = rMatrix ( alphas(2) , solution_IK(2 , 1))
R12ED = rMatrix ( alphas(2) , solution_IK(2 , 2))
R12B  = rMatrix ( alphas(2) , solution_IK(2 , 3))
R12BED= rMatrix ( alphas(2) , solution_IK(2 , 4))


R23   = rMatrix ( alphas(3) , solution_IK(3 , 1))
R23ED = rMatrix ( alphas(3) , solution_IK(3 , 2))
R23B  = rMatrix ( alphas(3) , solution_IK(3 , 3))
R23BED= rMatrix ( alphas(3) , solution_IK(3 , 4))


R03   = R01*R12*R23
R03ED = R01ED*R12ED*R23ED
R03B  = R01B*R12B*R23B
R03BED= R01BED*R12BED*R23BED

R36   = (inv(R03))*R06
R36ED = (inv(R03ED))*R06
R36B  = (inv(R03B))*R06
R36BED= (inv(R03BED))*R06

       


%1
Beta1    = atan2(  sqrt((1-(R36(3,3) - eps )^2)) , R36(3,3) - eps)
Beta12   = atan2( -sqrt((1-(R36(3,3) - eps )^2)) , R36(3,3) - eps)
%2
BetaED  = atan2(  sqrt(1-(R36ED(3,3) - eps )^2) , R36ED(3,3) - eps)
BetaED2 = atan2( -sqrt(1-(R36ED(3,3) - eps )^2) , R36ED(3,3) - eps)
%3
BetaZ   = atan2(  sqrt(1-(R36B(3,3) - eps )^2) , R36B(3,3) - eps)
BetaZ2  = atan2( -sqrt(1-(R36B(3,3) - eps )^2) , R36B(3,3) - eps)
%4
BetaBED = atan2(  sqrt((1-(R36BED(3,3) - eps )^2)) , R36BED(3,3) - eps)
BetaBED2= atan2( -sqrt((1-(R36BED(3,3) - eps )^2)) , R36BED(3,3) - eps)


%1
Alpha4 = atan2 ( (R36(2,3) - eps)/sin(Beta1) , (R36(1,3)-eps)/sin(Beta1))
Alpha4B = atan2 ( (-R36(2,3) - eps)/sin(Beta12) , (-R36(1,3)-eps)/sin(Beta12))
%2
Alpha4ED = atan2 ( (R36ED(2,3) - eps)/sin(BetaED) , (R36ED(1,3)-eps)/sin(BetaED))
Alpha4BED = atan2 ( (-R36ED(2,3) - eps)/sin(BetaED2) , (-R36ED(1,3)-eps)/sin(BetaED2))
%3
Alpha4Z = atan2 ( (R36B(2,3) - eps)/sin(BetaZ) , (R36B(1,3)-eps)/sin(BetaZ))
Alpha4BZ = atan2 ( (-R36B(2,3) - eps)/sin(BetaZ2) , (-R36B(1,3)-eps)/sin(BetaZ2))
%4
Alpha4Z2 = atan2 ( (R36BED(2,3) - eps)/sin(BetaBED) , (R36BED(1,3)-eps)/sin(BetaBED))
Alpha4BZ2 = atan2 ( (-R36BED(2,3) - eps)/sin(BetaBED2) , (-R36BED(1,3)-eps)/sin(BetaBED2))


%1
Gamma = atan2 ( (R36(3,2) - eps)/sin(Beta1) , (-R36(3,1)-eps)/sin(Beta1))
GammaB = atan2 ( (-R36(3,2) - eps)/sin(Beta12) , (R36(3,1)-eps)/sin(Beta12))
%2
GammaED = atan2 ( (R36ED(3,2) - eps)/sin(BetaED) , (-R36ED(3,1)-eps)/sin(BetaED))
GammaED2 = atan2 ( (-R36ED(3,2) - eps)/sin(BetaED2) , (R36ED(3,1)-eps)/sin(BetaED2))
%3
GammaZ = atan2 ( (R36B(3,2) - eps)/sin(BetaZ) , (-R36B(3,1)-eps)/sin(BetaZ))
GammaBZ = atan2 ( (-R36B(3,2) - eps)/sin(BetaZ2) , (R36B(3,1)-eps)/sin(BetaZ2))
%4
GammaZ2 = atan2 ( (R36BED(3,2) - eps)/sin(BetaBED) , (-R36BED(3,1)-eps)/sin(BetaBED))
GammaBZ2 = atan2 ( (-R36BED(3,2) - eps)/sin(BetaBED2) , (R36BED(3,1)-eps)/sin(BetaBED2))





%1
solution_IK(4 , 1) = Alpha4
solution_IK(5 , 1) = Beta1
solution_IK(6 , 1) = Gamma
solution_IK(7 , 1) = Alpha4B  -pi
solution_IK(8 , 1) = Beta12
solution_IK(9 , 1) = GammaB   -pi
%2
solution_IK(4 , 2) = Alpha4ED
solution_IK(5 , 2) = BetaED
solution_IK(6 , 2) = GammaED
solution_IK(7 , 2) = Alpha4BED -pi
solution_IK(8 , 2) = BetaED2
solution_IK(9 , 2) = GammaED2   -pi
%3
solution_IK(4 , 3) = Alpha4Z
solution_IK(5 , 3) = BetaZ
solution_IK(6 , 3) = GammaZ
solution_IK(7 , 3) = Alpha4BZ -pi
solution_IK(8 , 3) = BetaZ2
solution_IK(9 , 3) = GammaBZ -pi 
%4
solution_IK(4 , 4) = Alpha4Z2
solution_IK(5 , 4) = BetaBED
solution_IK(6 , 4) = GammaZ2
solution_IK(7 , 4) = Alpha4BZ2 -pi
solution_IK(8 , 4) = BetaBED2
solution_IK(9 , 4) = GammaBZ2 






