% clear command window
clc;

% Part I:

% case 1 parameters
L1 = 2000;
a1 = 200;
b1 = 400;
d1 = 100;

% case 2 parameters
L2 = 1200;
a2 = 400;
b2 = 600;
d2 = 150;

% 1.a.

figure();

% plot tank for case 1 in the x-z plane along y
[x1,y1,z1] = cylinder(1, 100);
c1 = surf(x1*a1, y1*b1, z1*L1-L1/2);
rotate(c1, [1,0,0], 90);
hold on;

% plot cylinder furthest to the left (greatest y value)
[x2,y2,z2] = cylinder(1, 10);
c2 = surf(x2*d1/2+a1/2, y2*d1/2-3*d1, z2*2*a1-a1/2);
rotate(c2, [0,1,0], 90);
hold on;

% plot next cylinder
[x3,y3,z3] = cylinder(1, 10);
c3 = surf(x3*d1/2+a1/2, y3*d1/2-d1, z3*2*a1-a1/2);
rotate(c3, [0,1,0], 90);
hold on;

% plot next cylinder
[x4,y4,z4] = cylinder(1, 10);
c4 = surf(x4*d1/2+a1/2, y4*d1/2+d1, z4*2*a1-a1/2);
rotate(c4, [0,1,0], 90);
hold on;

% plot cylinder furthest to the right (lowest y value)
[x5,y5,z5] = cylinder(1, 10);
c5 = surf(x5*d1/2+a1/2, y5*d1/2+3*d1, z5*2*a1-a1/2);
rotate(c5, [0,1,0], 90);
hold off;

% add labels and title
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');
title('Plot of Tank and Cylinders for Case 1');

figure();

% plot tank for case 2 in the x-z plane along y
[x1,y1,z1] = cylinder(1, 100);
c1 = surf(x1*a2, y1*b2, z1*L2-L2/2);
rotate(c1, [1,0,0], 90);
hold on;

% plot cylinder furthest to the left (greatest y value)
[x2,y2,z2] = cylinder(1, 10);
c2 = surf(x2*d2/2, y2*d2/2-2*d2, z2*2*a2-a2);
rotate(c2, [0,1,0], 90);
hold on;

% plot next cylinder
[x3,y3,z3] = cylinder(1, 10);
c3 = surf(x3*d2/2, y3*d2/2-a2/4, z3*2*a2-a2);
rotate(c3, [0,1,0], 90);
hold on;

% plot next cylinder
[x4,y4,z4] = cylinder(1, 10);
c4 = surf(x4*d2/2, y4*d2/2+a2/4, z4*2*a2-a2);
rotate(c4, [0,1,0], 90);
hold on;

% plot cylinder furthest to the right (lowest y value)
[x5,y5,z5] = cylinder(1, 10);
c5 = surf(x5*d2/2, y5*d2/2+2*d2, z5*2*a2-a2);
rotate(c5, [0,1,0], 90);
hold off;

% add labels and title
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');
title('Plot of Tank and Cylinders for Case 2');

% 1.b.

% triple integral expressions are derived in the report

% 1.c.

syms theta L r x;

% calculate volume of shell, tube, and cap, and total volume for case 1

Vshell_1 = double(int(int(int(r, r, 0, sqrt((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))), L, 0, L1), theta, 0, 2*pi) / 10^9);
sprintf('Volume of Shell for Case 1: %f m^3', Vshell_1)

Vtube_1 = double(int(int(int(r, r, 0, d1/2), theta, 0, 2*pi), x, -a1, a1) / 1000^3);
sprintf('Volume of Tube for Case 1: %f m^3', Vtube_1)

Vcap_1 = double(int(int(int(r, x, a1*sqrt(1-(r^2*(cos(theta)^2)/b1^2)), a1), r, 0, d1/2), theta, 0, 2*pi) / 10^9);
sprintf('Volume of Cap for Case 1: %f m^3', Vcap_1)

Vtotal_1 = Vshell_1 - 4*Vtube_1 + 8*Vcap_1;
sprintf('Total Volume for Case 1: %f m^3', Vtotal_1)

% calculate volume of shell, tube, and cap, and total volume for case 2

Vshell_2 = double(int(int(int(r, r, 0, sqrt((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))), L, 0, L2), theta, 0, 2*pi) / 10^9);
sprintf('Volume of Shell for Case 2: %f m^3', Vshell_2)

Vtube_2 = double(int(int(int(r, r, 0, d2/2), theta, 0, 2*pi), x, -a2, a2) / 1000^3);
sprintf('Volume of Tube for Case 2: %f m^3', Vtube_2)

Vcap_2 = double(int(int(int(r, x, a2*sqrt(1-(r^2*(cos(theta)^2)/b2^2)), a2), r, 0, d2/2), theta, 0, 2*pi) / 10^9);
sprintf('Volume of Cap for Case 2: %f m^3', Vcap_2)

Vtotal_2 = Vshell_2 - 4*Vtube_2 + 8*Vcap_2;
sprintf('Total Volume for Case 2: %f m^3', Vtotal_2)

% 2.a.

% calculate area of tube cross section in x-z plane for case 1
Atube_1 = double(int(int(r, r, 0, sqrt((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))), theta, 0, 2*pi) / 10^6);
sprintf('Area of Tube Cross Section for Case 1: %f m^2', Atube_1)

% calculate area of tube cross section in x-z plane for case 2
Atube_2 = double(int(int(r, r, 0, sqrt((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))), theta, 0, 2*pi) / 10^6);
sprintf('Area of Tube Cross Section for Case 2: %f m^2', Atube_2)

% 2.b.

% triple integral expressions are derived in the report

% 2.c.

% calculate moments about the x-axis and z-axis for case 1
Mx_1 = double(int(int(r^2*sin(theta)*(420+140*(cos(theta)-sin(theta))), r, 0, sqrt((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))), theta, 0, 2*pi));
Mz_1 = double(int(int(r^2*cos(theta)*(420+140*(cos(theta)-sin(theta))), r, 0, sqrt((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))), theta, 0, 2*pi));

% calculate moments about the x-axis and z-axis for case 2
Mx_2 = double(int(int(r^2*sin(theta)*(420+140*(cos(theta)-sin(theta))), r, 0, sqrt((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))), theta, 0, 2*pi));
Mz_2 = double(int(int(r^2*cos(theta)*(420+140*(cos(theta)-sin(theta))), r, 0, sqrt((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))), theta, 0, 2*pi));

% calculate mass for case 1 and 2
M_1 = double(int(int(r*(420+140*(cos(theta)-sin(theta))), r, 0, sqrt((a1^2*b1^2)/(a1^2*(sin(theta))^2+b1^2*(cos(theta))^2))), theta, 0, 2*pi));
M_2 = double(int(int(r*(420+140*(cos(theta)-sin(theta))), r, 0, sqrt((a2^2*b2^2)/(a2^2*(sin(theta))^2+b2^2*(cos(theta))^2))), theta, 0, 2*pi));

% calculate center of mass for case 1
x_bar_1 = Mz_1 / M_1 / 10^3;
z_bar_1 = Mx_1 / M_1 / 10^3;
sprintf('Center of Mass for Case 1: x: %f m^2 z: %f m^2', x_bar_1, z_bar_1)

% calculate center of mass for case 2
x_bar_2 = Mz_2 / M_2 / 10^3;
z_bar_2 = Mx_2 / M_2 / 10^3;
sprintf('Center of Mass for Case 2: x: %f m^2 z: %f m^2', x_bar_2, z_bar_2)
