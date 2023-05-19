%% Variable Initialization
format longG;
% Main 
%Rotating Base
motor1HomeingTick = 2048; % DH joint angle unit rad
%Motor that is connected between Link 1 and Link 2
motor2HomeingTick = 1180;
%Motor between Link2 and link 3
motor3HomeingTick = 3700;
%Motor between Link 3 and gripper
motor4HomeingTick = 2048;
%motor that opens and closes claw (Need experimentation
motor5HomeingTick = 2048;

gripperOpenTick = 2048;
gripperClosedTick = 1400;

convertFromRad2Tick = 4096 / (2*pi);
convertFromTick2Rad = (2*pi) / 4096;
% All measurements are in Millimeters

% Defining the Possible Positions for Bottle
bottleHeight = 250;
bottleHorizontalDisplacement = 254;
bottleVerticalDisplacement = 254;

baseHeight = 114.9;


phiBottlePosition = 0;

%Distance between frame 3 and 2 along x axis
L1 = 250;
%distance between frame 4 and frame 3 along x axis
L2 = 200;
%X distance between frame 4 and gripper position
L3 = 96.22;

alpha1 = pi/2;

%Z link displacement from bottom of base to Link 2
d1 = baseHeight;
%Z link displacement from frame 2 to frame 3
d3 = -62.75;
% Z link displacement from frame 3 to frame 4
d4 = 0;
%Z link displacement from frame 4 to frame 5
d5 = 3.78;
%DH parameter for x distance between frame 3 and 2
a2 = L1;
%DH parameter for x distance between frame 4 and frame 3
a3 = L2;
%Dh parameter for x distance between gripper and frame 4
a4 = L3;

%Getting current motor positioin
currentMotorPosition1 = motor1HomeingTick;
currentMotorPosition2 = motor2HomeingTick;
currentMotorPosition3 = motor3HomeingTick;
currentMotorPosition4 = motor4HomeingTick;

%Getting Position where it should be with respect to frames
currentTheta1 = deg2rad(0);
currentTheta2 = deg2rad(90);
currentTheta3 = deg2rad(-145.1953);
currentTheta4 = deg2rad(0);

% Matrix of motor positions created to aid in code testing
CurrentMotorPositions = [currentMotorPosition1 currentMotorPosition2 currentMotorPosition3 currentMotorPosition4];
CurrentThetas = [currentTheta1 currentTheta2 currentTheta3 currentTheta4];

% Motor offsets implemented due to shape of limbs
motor3TickOffset = -50;
motor4TickOffset = -130;

%Calculates XYZ Position of Homing
HomePosition = forwardKinematics(currentTheta1,currentTheta2,currentTheta3,currentTheta4,a2,a3,a4,d1,d3,d5);
%% Position 1 Calculations
%Inverse Kinematics Portion
% XYZ coordinates for First Position
x = bottleHorizontalDisplacement;
y = 0;
z = bottleHeight;

% Inverse Kinematics for Point 1
[goaltheta1,goaltheta2,goaltheta3,goaltheta4] = inverseKinematics(x,y,z,phiBottlePosition,a2,a3,a4,d1, d3,d5);

% Matrix of motor positions and thetas created to aid in code testing
goalThetas1 = [goaltheta1 goaltheta2 goaltheta3 goaltheta4];

XYZTrajectoryPoint1 = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5);

currentTheta1 = goaltheta1;
currentTheta2 = goaltheta2;
currentTheta3 = goaltheta3;
currentTheta4 = goaltheta4;

currentThetasPoint1 = [currentTheta1 currentTheta2 currentTheta3 currentTheta4];

x = 0;
y = 400;
z = 400;

% Inverse Kinematics for Point 2
[goaltheta1,goaltheta2,goaltheta3,goaltheta4] = inverseKinematics(x,y,z,phiBottlePosition,a2,a3,a4,d1, d3,d5);

% Matrix of motor positions and thetas created to aid in code testing
goalThetas2 = [goaltheta1 goaltheta2 goaltheta3 goaltheta4];

XYZTrajectoryPoint2 = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5);

RobotLocationWithRelationToBase2 = forwardKinematics(goaltheta1,goaltheta2,goaltheta3,goaltheta4,a2,a3,a4,d1,d3,d5);


currentTheta1 = goaltheta1;
currentTheta2 = goaltheta2;
currentTheta3 = goaltheta3;
currentTheta4 = goaltheta4;

currentThetasPoint2 = [currentTheta1 currentTheta2 currentTheta3 currentTheta4];
%% Position 3 Calculations
pause(3);
x = -400;
y = 0;
z = 400;

% Inverse Kinematics for Point 3
[goaltheta1,goaltheta2,goaltheta3,goaltheta4] = inverseKinematics(x,y,z,phiBottlePosition,a2,a3,a4,d1, d3,d5);


% Matrix of motor positions and thetas created to aid in code testing
goalThetas3 = [goaltheta1 goaltheta2 goaltheta3 goaltheta4];

XYZTrajectoryPoint3 = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5);

RobotLocationWithRelationToBase2 = forwardKinematics(goaltheta1,goaltheta2,goaltheta3,goaltheta4,a2,a3,a4,d1,d3,d5);



currentTheta1 = goaltheta1;
currentTheta2 = goaltheta2;
currentTheta3 = goaltheta3;
currentTheta4 = goaltheta4;

currentThetasPoint3 = [currentTheta1 currentTheta2 currentTheta3 currentTheta4];
%% Position 4 Caluclations
z = bottleHeight;
x = -1.5*bottleHorizontalDisplacement;
y = 0;

% Inverse Kinematics for Point 4
[goaltheta1,goaltheta2,goaltheta3,goaltheta4] = inverseKinematics(x,y,z,phiBottlePosition,a2,a3,a4,d1, d3,d5);


%
% Matrix of motor positions and thetas created to aid in code testing
goalThetas4 = [goaltheta1 goaltheta2 goaltheta3 goaltheta4];

XYZTrajectoryPoint4 = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5);

RobotLocationWithRelationToBase3 = forwardKinematics(goaltheta1,goaltheta2,goaltheta3,goaltheta4,a2,a3,a4,d1,d3,d5);


currentTheta1 = goaltheta1;
currentTheta2 = goaltheta2;
currentTheta3 = goaltheta3;
currentTheta4 = goaltheta4;

currentThetasPoint4 = [currentTheta1 currentTheta2 currentTheta3 currentTheta4];
%% Point 5 Calculations

x = -1.25*bottleHorizontalDisplacement;
y = 0;
z = 400;

% Inverse Kinematics for Point 5
[goaltheta1,goaltheta2,goaltheta3,goaltheta4] = inverseKinematics(x,y,z,phiBottlePosition,a2,a3,a4,d1, d3,d5);

% Matrix of motor positions and thetas created to aid in code testing
goalThetas5 = [goaltheta1 goaltheta2 goaltheta3 goaltheta4];

XYZTrajectoryPoint5 = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5);

RobotLocationWithRelationToBase5 = forwardKinematics(goaltheta1,goaltheta2,goaltheta3,goaltheta4,a2,a3,a4,d1,d3,d5);

%% Back to home
currentTheta1 = goaltheta1;
currentTheta2 = goaltheta2;
currentTheta3 = goaltheta3;
currentTheta4 = goaltheta4;

goaltheta1 = deg2rad(0);
goaltheta2 = deg2rad(90);
goaltheta3 = deg2rad(-145.1953);
goaltheta4 = deg2rad(0);

XYZTrajectoryPointHome = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5);




function PointSet = ObtainPoints(goaltheta1, goaltheta2, goaltheta3, goaltheta4, currentTheta1, currentTheta2,currentTheta3, currentTheta4,a2,a3,a4,d1,d3,d5)

s1 = linspace(currentTheta1,goaltheta1,15);
s2 = linspace(currentTheta2,goaltheta2,15);
s3 = linspace(currentTheta3,goaltheta3,15);
s4 = linspace(currentTheta4,goaltheta4,15);

x = zeros(1,15);
y = zeros(1,15);
z = zeros(1,15);

for i = 1:15
    ForwardKinematicsMatrix = forwardKinematics(s1(i),s2(i),s3(i),s4(i),a2,a3,a4,d1,d3,d5);
    x(i) = ForwardKinematicsMatrix(1,4);
    y(i) = ForwardKinematicsMatrix(2,4);
    z(i) = ForwardKinematicsMatrix(3,4);
end

PointSet = [x;y;z];
end




