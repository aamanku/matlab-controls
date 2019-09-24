%% This function gives the state variables after the impact of swing foot with ground
%%qn=impacttransition(q)
function qn=impacttransition(q)
% syms m r MH MT L g theta1 theta2 theta3 dtheta1 dtheta2 dtheta3 fx Lfh theta3d real   
[m,MH,MT,r,L,g,theta3d,theta1d,alphaa,epsilon,dtheta1_]=simulationparameters();
theta1=q(1);
theta2=q(2); 
theta3=q(3);

%% De matrix
De=[(5/4*m+MH+MT)*r^2 -1/2*m*cos(theta1-theta2)*r^2 MT*r*L*cos(theta1-theta3) (3/2*m+MH+MT)*r*cos(theta1) -(3/2*m+MH+MT)*r*sin(theta1);
    -1/2*m*cos(theta1-theta2)*r^2 1/4*m*r^2 0 -1/2*m*r*cos(theta2) 1/2*m*r*sin(theta2);
    MT*r*L*cos(theta1-theta3) 0 MT*L^2 MT*L*cos(theta3) -MT*L*sin(theta3);
    (3/2*m+MH+MT)*r*cos(theta1) -1/2*m*r*cos(theta2) MT*L*cos(theta3) 2*m+MH+MT 0;
    -(3/2*m+MH+MT)*r*sin(theta1) 1/2*m*r*sin(theta2) -MT*L*sin(theta3) 0 2*m+MH+MT];
% De=(De+De')
% De=0.5*eye(5).*De
%% E matrix
E=[r*cos(theta1) -r*cos(theta2) 0 1 0;
    -r*sin(theta1) r*sin(theta2) 0 0 1];
%%
qa=([De -E';E zeros(2)])\[De*[q(4:6)';zeros(2,1)];zeros(2,1)];

qn(1)=q(2);
qn(2)=q(1);
qn(3)=q(3);
qn(4)=qa(2);
qn(5)=qa(1);
qn(6)=qa(3);
% rad2deg(qn)


end