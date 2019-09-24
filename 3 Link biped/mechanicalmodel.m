function mod=mechanicalmodel(~,x)
%% initial
%  syms m r MH MT L g theta1 theta2 theta3 dtheta1 dtheta2 dtheta3 fx Lfh theta3d
%x=x';
theta1=x(1);
theta2=x(2);
theta3=x(3);
dtheta1=x(4);
dtheta2=x(5);
dtheta3=x(6);
m=5;   
MH=15; 
MT=10;
r=1;
L=0.5; 
g=9.8; 
theta3d=pi/6;
theta1d=pi/8;

%% D matrix
D=[(5/4*m+MH+MT)*r^2 -1/2*r^2*cos(theta1-theta2)*m MT*r*L*cos(theta1-theta3);
    -1/2*r^2*cos(theta1-theta2)*m 1/4*m*r^2 0; 
    MT*r*L*cos(theta1-theta3) 0 MT*L^2];
%% C matrix
C=[0 -1/2*sin(theta1-theta2)*m*dtheta2*r^2 MT*r*L*sin(theta1-theta3)*dtheta3;
    1/2*m*sin(theta1-theta2)*dtheta1*r^2 0 0;
    -MT*r*L*sin(theta1-theta3)*dtheta1 0 0];
%% G matrix
G=[-1/2*g*(2*MH+3*m+2*MT)*r*sin(theta1);
    1/2*g*m*r*sin(theta2);
    -g*MT*L*sin(theta3)];
%% B matrix
B=[-1 0;
    0 -1;
    1 1];
%% h matrix
h=[theta3-theta3d;
    theta2+theta1;];
%% final function
mod=zeros(6,1);
mod(1:3)=x(4:6);
fx=(D)\(-C*x(4:6)-G);
% Lfh=jacobian(fx,[dtheta2 dtheta3])
% simplify(Lfh)
gx=(D)\(B);
u=[-(theta3d-x(3))*1;-(theta3d+x(3))*5];
mod(4:6)=fx+gx*u;
% mod=mod';

end