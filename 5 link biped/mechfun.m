function [Q,u]=mechfun(t,x)
global g lT lf lt MT Mf Mt IT If It pMT pMf pMt thetap thetam a
x=x(:);
D=D_q_m(x);
% size(D)
G=G_q_m(x);
% size(G)
C=C_q_m(x(1:5),x(6:10));
% size(C)
B = [zeros(1,4);eye(4)];
fx=[x(6:10);D\(-C*x(6:10)-G)];
gx = [zeros(5,4);D\B];
v=PSIi(h(t,x,a,thetam,thetap),Lfh(t,x,a,thetam,thetap),0.9,0.1);
u = LgLfh(t,x,a,thetam,thetap)\(v - Lf2h(t,x,a,thetam,thetap));
Q = fx + gx*u;
end

