function q_n = impacttransition(q)
De = D_e(q);
Ce = C_e(q(1:5),q(6:10));
Ge = G_e(q);
Eq = E(q);
dqe_minus = [q(6) q(7) q(8) q(9) q(10) 0 0]';
dqe_plus = [De -Eq';Eq zeros(2,2)]\[De*dqe_minus;0;0];
q_n = [q(1); q(3); q(2); q(5); q(4); dqe_plus(1); dqe_plus(3); dqe_plus(2); dqe_plus(5); dqe_plus(4)];
