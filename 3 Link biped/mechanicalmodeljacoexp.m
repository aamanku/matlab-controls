%% This function makes mechanical model's expression and saves to
%%mechanicalfun file which is later called in ode45 for faster procesesing.
%%it calculates only one time
%%=mechanicalmodeljacoexp()
function funmod=mechanicalmodeljacoexp()
%% initial
syms m r MH MT L g theta1 theta2 theta3 dtheta1 dtheta2 dtheta3 fx Lfh th3d v v1 v2 real 
[m,MH,MT,r,L,g,theta3d,theta1d,alphaa,epsilon,dtheta1_]=simulationparameters();

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
h=[theta3-theta3d ;theta2+theta1];
%% final function
% mod=zeros(6,1);
% mod(1:3)=x(4:6);
fx=[dtheta1;dtheta2;dtheta3;(D)\(-C*[dtheta1 dtheta2 dtheta3]'-G)];
% size(fx)
gx=[zeros(3,2);(D)\(B)];
% size(gx)
Lfh=jacobian(h,[theta1 theta2 theta3 dtheta1 dtheta2 dtheta3])*fx;
% size(Lfh)
Lf2h=jacobian(Lfh,[theta1 theta2 theta3 dtheta1 dtheta2 dtheta3])*fx;
% size(Lf2h)
LgLfh=jacobian(Lfh,[theta1 theta2 theta3 dtheta1 dtheta2 dtheta3])*gx;
% size(LgLfh)
%  simplify(Lfh)
v=PSIi(h,Lfh,alphaa,epsilon);
% simplify(subs(Lf2h))
% v=[v1;v2]
u=LgLfh\(v-Lf2h);
% mod=(fx+gx*u);
% mod=double(subs(fx+gx*u));
% mod=mod';
rr=subs(fx+gx*u);

% funmod=matlabFunction(rr,'File','mechanicalfun','Vars',[t,q]);
funmod=matlabFunction(rr,'File','mechanicalfun','Vars',[theta1 theta2 theta3 dtheta1 dtheta2 dtheta3]);
end