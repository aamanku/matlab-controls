function F = contact_force(q,u)
dqe = [q(6:10) ; 0 ; 0];
B = [zeros(1,4); eye(4); zeros(2,4)];
D = D_e(q);
C = C_e(q(1:5),q(6:10));
G = G_e(q);
% fx = [dq;D\(-C*dq-G)];
% gx = [zeros(5,4);D\B];
J = [0 0 0 0 0 1 0;
    0 0 0 0 0 0 1];
F = (J*(D\J'))\(J*(D\(C*dqe+G-B*u)));