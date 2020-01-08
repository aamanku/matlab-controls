function jacsym()
a=sym('a',[4 7]);
x=sym('x',[10 1]);
syms thm thp t
D=D_q_m(x);
%  size(D)
G=G_q_m(x);
% size(G)
C=C_q_m(x(1:5),x(6:10));
% size(C)
B = [zeros(1,4);eye(4)];
% size(B)
invD=(inv(D));
fx=([x(6:10);invD*(-C*x(6:10)-G)]);
% size(fx)
gx = ([zeros(5,4);invD*B]);
c = [1 1 0 0.5 0];
th=c*x(1:5);
s = (th - thp)/(thm - thp);
a1=a(1,:);a2=a(2,:);a3=a(3,:);a4=a(4,:);
y1 = x(2) - bezier(s,6,a1);
y2 = x(3) - bezier(s,6,a2);
y3 = x(4) - bezier(s,6,a3);
y4 = x(5) - bezier(s,6,a4);
h = ([y1; y2; y3; y4]);
Lfh = (jacobian(h,x)*fx);
% size(Lfh)
Lf2h=(jacobian(Lfh,x)*fx);
% size(Lf2h)
LgLfh=(jacobian(Lfh,x)*gx);
matlabFunction(Lfh,'File','Lfh','Vars',{t,x,a,thm,thp},'Optimize',false);
matlabFunction(Lf2h,'File','Lf2h','Vars',{t,x,a,thm,thp},'Optimize',false);
matlabFunction(LgLfh,'File','LgLfh','Vars',{t,x,a,thm,thp},'Optimize',false);
matlabFunction(h,'File','h','Vars',{t,x,a,thm,thp},'Optimize',false);
end

function y = bezier(s,M,a)
y = 0;
for k = 0:M
    coff = factorial(M)/( factorial(k)*factorial(M-k) );
    b = a(k+1)*coff*(s.^k).*((1-s).^(M-k));
    y = y + b;
end
end

