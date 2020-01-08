function [G] = G_e(q)

  G(1,1)=(477456*sin(q(1) + q(2) + q(4)))/3125 - (12544*sin(q(1) + q(3) + q(5)))/3125 + (1667421*...
         sin(q(1) + q(2)))/12500 - (292579*sin(q(1) + q(3)))/12500 - (196*sin(q(1)))/5;
  G(2,1)=(477456*sin(q(1) + q(2) + q(4)))/3125 + (1667421*sin(q(1) + q(2)))/12500;
  G(3,1)=- (12544*sin(q(1) + q(3) + q(5)))/3125 - (292579*sin(q(1) + q(3)))/12500;
  G(4,1)=(477456*sin(q(1) + q(2) + q(4)))/3125;
  G(5,1)=-(12544*sin(q(1) + q(3) + q(5)))/3125;
  G(6,1)=0;
  G(7,1)=392;

 