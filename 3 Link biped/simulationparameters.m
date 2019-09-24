%% Call this function to get simulation parameters
%%[m,MH,MT,r,L,g,theta3d,theta1d,alphaa,epsilon,dtheta1_]=simulationparameters()
function [m,MH,MT,r,L,g,theta3d,theta1d,alphaa,epsilon,dtheta1_]=simulationparameters()
m=5;   
MH=15; 
MT=10; 
r=1;
L=0.5; 
g=9.8; 
theta3d=pi/6;
theta1d=pi/8;
alphaa=0.9;
epsilon=0.1;
dtheta1_=1.6;
end