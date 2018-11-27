% clear command window
clc;

% Part I:

L1 = 2000;
a1 = 200;
b1 = 400;
d1 = 100;

L2 = 1200;
a2 = 400;
b2 = 600;
d2 = 150;

% 1.a.

figure();

[x1,y1,z1] = cylinder(1, 100);
c1 = surf(x1*a1, y1*b1, z1*L1-L1/2);
rotate(c1, [1,0,0], 90);
hold on;

[x2,y2,z2] = cylinder(1, 10);
c2 = surf(x2*d1/2+a1/2, y2*d1/2-3*d1, z2*2*a1-a1/2);
rotate(c2, [0,1,0], 90);
hold on;

[x3,y3,z3] = cylinder(1, 10);
c3 = surf(x3*d1/2+a1/2, y3*d1/2-d1, z3*2*a1-a1/2);
rotate(c3, [0,1,0], 90);
hold on;

[x4,y4,z4] = cylinder(1, 10);
c4 = surf(x4*d1/2+a1/2, y4*d1/2+d1, z4*2*a1-a1/2);
rotate(c4, [0,1,0], 90);
hold on;

[x5,y5,z5] = cylinder(1, 10);
c5 = surf(x5*d1/2+a1/2, y5*d1/2+3*d1, z5*2*a1-a1/2);
rotate(c5, [0,1,0], 90);
hold off;

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Plot of Tank and Cylinders for Case 1');

figure();

[x1,y1,z1] = cylinder(1, 100);
c1 = surf(x1*a2, y1*b2, z1*L2-L2/2);
rotate(c1, [1,0,0], 90);
hold on;

[x2,y2,z2] = cylinder(1, 10);
c2 = surf(x2*d2/2, y2*d2/2-2*d2, z2*2*a2-a2);
rotate(c2, [0,1,0], 90);
hold on;

[x3,y3,z3] = cylinder(1, 10);
c3 = surf(x3*d2/2, y3*d2/2-a2/4, z3*2*a2-a2);
rotate(c3, [0,1,0], 90);
hold on;

[x4,y4,z4] = cylinder(1, 10);
c4 = surf(x4*d2/2, y4*d2/2+a2/4, z4*2*a2-a2);
rotate(c4, [0,1,0], 90);
hold on;

[x5,y5,z5] = cylinder(1, 10);
c5 = surf(x5*d2/2, y5*d2/2+2*d2, z5*2*a2-a2);
rotate(c5, [0,1,0], 90);
hold off;

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Plot of Tank and Cylinders for Case 2');

% 1.b.

% Setup triple integral expressions

% 1.c.

syms theta L r x;

Vshell_1 = double(int(int(int(r, r, 0, sqrt((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))), L, 0, L1), theta, 0, 2*pi) / 1000^3);
sprintf('Volume of Shell for Case 1: %f m^3', Vshell_1)

Vtube_1 = double(int(int(int(r, r, 0, d1/2), theta, 0, 2*pi), x, -a1, a1) / 1000^3);
sprintf('Volume of Tube for Case 1: %f m^3', Vtube_1)

Vcap_1 = double(int(int(int(r, x, a1*sqrt((1-(((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))*(cos(theta)^2))/b1^2)), a1), r, 0, d1/2), theta, 0, 2*pi) / 1000^3);
sprintf('Volume of Cap for Case 1: %f m^3', Vcap_1)

Vtotal_1 = Vshell_1 + Vtube_1 + Vcap_1;
sprintf('Total Volume for Case 1: %f m^3', Vtotal_1)

Vshell_2 = double(int(int(int(r, r, 0, sqrt((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))), L, 0, L2), theta, 0, 2*pi) / 1000^3);
sprintf('Volume of Shell for Case 2: %f m^3', Vshell_2)

Vtube_2 = double(int(int(int(r, r, 0, d2/2), theta, 0, 2*pi), x, -a2, a2) / 1000^3);
sprintf('Volume of Tube for Case 2: %f m^3', Vtube_2)

Vcap_2 = double(int(int(int(r, x, a2*sqrt((1-(((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))*(cos(theta)^2))/b2^2)), a2), r, 0, d2/2), theta, 0, 2*pi) / 1000^3);
sprintf('Volume of Cap for Case 2: %f m^3', Vcap_2)

Vtotal_2 = Vshell_2 + Vtube_2 + Vcap_2;
sprintf('Total Volume for Case 2: %f m^3', Vtotal_2)

% 2.a.

Atube_1 = pi * a1 * b1 / 1000^2;
sprintf('Area of Tube Cross Section for Case 1: %f m^2', Atube_1)

Atube_2 = pi * a2 * b2 / 1000^2;
sprintf('Area of Tube Cross Section for Case 1: %f m^2', Atube_2)

% 2.b.

% Setup triple integral expressions

% 2.c.


