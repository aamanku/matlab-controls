function extended_matrices()
syms q1 q2 q3 q4 q5 dq1 dq2 dq3 dq4 dq5 g lT lf lt MT Mf Mt IT If It pMT pMf pMt z1 z2 z1d z2d real
q = [q1;q2;q3;q4;q5;z1;z2];
dq = [dq1; dq2; dq3;dq4;dq5;z1d;z2d];
% g=9.8;
% lT=0.625;
% lf=0.4;
% lt=0.4;
% MT = 20;
% Mf = 6.8;
% Mt = 3.2;
% IT = 2.22;
% If = 1.08;
% It = 0.93;
% pMT = 0.2;
% pMf = 0.163;
% pMt = 0.128;
global g lT lf lt MT Mf Mt IT If It pMT pMf pMt thetap thetam a
%position of center of mass of each link as a function of states
p_4 = [z1-(lt-pMt)*sin(pi-q1-q2-q4);z2+(lt-pMt)*cos(pi-q1-q2-q4)];
v_4=jacobian(p_4,q)*dq;

p_2=[z1-lt*sin(pi-q1-q2-q4)-(lf-pMf)*sin(pi-q1-q2); z2+lt*cos(pi-q1-q2-q4)+(lf-pMf)*cos(pi-q1-q2)];
v_2=jacobian(p_2,q)*dq;

p_1=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+pMT*sin(q1); z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+pMT*cos(q1)];
v_1=jacobian(p_1,q)*dq;

p_3=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-pMf*sin(q1+q3-pi); z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-pMf*cos(q1+q3-pi)];
v_3=jacobian(p_3,q)*dq;

p_5=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-pMt*sin(q1+q3+q5-pi); z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-pMt*cos(q1+q3+q5-pi)];
v_5=jacobian(p_5,q)*dq;

% kinetic energy of links
KE_1=1/2*MT*v_1.'*v_1+0.5*IT*dq1^2;
KE_2=1/2*Mf*v_2.'*v_2+0.5*If*(dq1+dq2)^2;
KE_3=1/2*Mf*v_3.'*v_3+0.5*If*(dq1+dq3)^2;
KE_4=1/2*Mt*v_4.'*v_4+0.5*It*(dq1+dq2+dq4)^2;
KE_5=1/2*Mt*v_5.'*v_5+0.5*It*(dq1+dq3+dq5)^2;

% total kinetic energy
KE=KE_1+KE_2+KE_3+KE_4+KE_5;
KE=simplify(KE);


%potential energy of links

PE=g*(Mf*(p_2(2)+p_3(2))+Mt*(p_4(2)+p_5(2))+MT*p_1(2));
PE=simplify(PE);

% gravity vector
G=jacobian(PE,q)';
G=simplify(G);

% mass-inertial matrix
D=simplify(jacobian(KE,dq).');
D=jacobian(KE,dq).';
D=simplify(jacobian(D,dq));

% Coriolis and centrigugal matrix

syms C real
% C1=simplify(jacobian(jacobian(KE,dq),q))
% C2=-jacobian(KE,q).'
n=max(size(q));
for k=1:n
    for j=1:n
        C(k,j)=0*g;
        for i=1:n
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+...
                diff(D(k,i),q(j))-...
                diff(D(i,j),q(k)))*dq(i);
        end
    end
end
C = simplify(C);

% input matrix
%B=[eye(4);zeros(1,4)];

p2=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-lt*sin(q1+q3+q5-pi) ; lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-lt*cos(q1+q3+q5-pi)];
v2v=jacobian(p2,q)*dq;
v2v=simplify(v2v);

%Jacobian of force
%P_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)];  % at hip
P_Fext = [-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+lT*sin(q1) lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+lT*cos(q1)];  %at head
J = jacobian(P_Fext,q);
p_st = [z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-lt*sin(q1+q3+q5-pi); z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-lt*cos(q1+q3+q5-pi)];
E = jacobian(p_st,q);
%-------------------------------- Write the results
list_q_e  = {'q1','q(1)'; 'q2','q(2)'; 'q3','q(3)'; 'q4','q(4)'; ...
    'q5','q(5)'; 'dq1','q(6)';'dq2','q(7)';'dq3','q(8)';'dq4','q(9)';'dq5','q(10)'};
write_fcn('D_e.m',{'q'},[list_q_e],{D,'D'});
write_fcn('C_e.m',{'q','dq'},[list_q_e],{C,'C'});
write_fcn('G_e.m',{'q'},[list_q_e],{G,'G'});
write_fcn('E.m',{'q'},[list_q_e],{E,'E'});
end