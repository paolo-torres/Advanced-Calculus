% clear command window
clc;

% Part I:

% 1.a.

% define x and y region domains
x_range = -5:0.1:5;
y_range = -2:0.1:4;

% define terrain function
[x,y] = meshgrid(x_range, y_range);
z = log(x.^4+1).*(4.*x.^4+(2.*y).^2).*exp(-0.5.*x.^2-(y-0.8).^2-3).*(sin(2.*x+0.05.*y.^4)+2.*cos(0.75.*y));

% plot terrain
figure();
surf(x,y,z);

% plot title and labels
title('Plot of Terrain');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

% plot contours
figure();
[c,h] = contour(x,y,z,25);
clabel(c,h);
title('Contours of Terrain');
xlabel('X [km]');
ylabel('Y [km]');
colorbar;

% 1.b.

g_max_x = 0;
g_max_y = 0;
g_max_val = 0;

% loop through entire domain of surface
for x0 = -5:0.1:5
    for y0 = -2:0.1:4
        
        % calculate gradient at every point
        [fx,fy] = gradient(z,0.1);
        t = (x == x0) & (y == y0);
        indt = find(t);
        f_grad = [fx(indt) fy(indt)];
        
        % take the magnitude of the gradient
        norm_grad = norm(f_grad);
        
        % extract maximum gradient and its coordinates
        if (norm_grad > g_max_val)
            g_max_x = x0;
            g_max_y = y0;
            g_max_val = norm_grad;
        end
    end
end
sprintf('Max Slope X: %f, Max Slope Y: %f, Max Slope: %f', g_max_x, g_max_y, g_max_val)

% 1.c.

% define terrain function
[a,b] = meshgrid(x_range, y_range);
syms x y z;
z = log(x.^4+1).*(4.*x.^4+(2.*y).^2).*exp(-0.5.*x.^2-(y-0.8).^2-3).*(sin(2.*x+0.05.*y.^4)+2.*cos(0.75.*y));
z_func = matlabFunction(z);

% take first partial derivatives
dzdx = diff(z, x);
dzdy = diff(z, y);
val_dzdx = matlabFunction(dzdx);
val_dzdy = matlabFunction(dzdy);
c_val_dzdx = val_dzdx(a,b);
c_val_dzdy = val_dzdy(a,b);

% plot the contours of each partial derivative to determine intersections
[c,h] = contour(a, b, c_val_dzdx, [0.0 0.0], 'red');
clabel(c,h);
hold on;
[c,h] = contour(a, b, c_val_dzdy, [0.0 0.0], 'blue');
clabel(c,h);

% use the intersections to determine initial guesses to solve for critical points 
% (-4.2,2.8) (-4.0,1.8) (-3.9,-1.7) (-2.4,0.6) (-1.3,-1.8) (-1.6,2.2)
% (-0.4,2.1) (-0.7,2.5) (-1.8,3.1) (-1.6,3.3) (-1.6,3.8) (-1.3,3.9)
% (2.2,-1.8) (1.6,0.5) (2.3,0.4) (3.0,0.5) (2.0,1.9) (1.8,3.1)
% (2.9,3.0) (2.4,3.7) (2.6,3.9)

% x and y guesses stored as arrays
x_guesses = [-4.2, -4.0, -3.95, -2.4, -1.3, -1.6, -0.4, -0.7, -1.8, -1.6, -1.6, -1.3, 2.2, 1.6, 2.3, 3.0, 2.0, 1.8, 2.9, 2.4, 2.6];
y_guesses = [2.8, 1.8, -1.75, 0.6, -1.8, 2.2, 2.1, 2.5, 3.1, 3.3, 3.8, 3.9, -1.8, 0.5, 0.4, 0.5, 1.9, 3.1, 3.0, 3.7, 3.9];

values = [];
values_x = [];
values_y = [];

max_elevation = 0;
max_elevation_x = 0;
max_elevation_y = 0;
min_elevation = 0;
min_elevation_x = 0;
min_elevation_y = 0;

% iterate through for each guess
for i = 1:length(x_guesses)
    
    % solve partial derivatives for critical points
    S = vpasolve([dzdx, dzdy], [x, y], [x_guesses(i), y_guesses(i)]);
    values_x = [values_x; double(S.x)];
    values_y = [values_y; double(S.y)];
    values = [values; feval(z_func, double(S.x), double(S.y))];
    
    % extract highest and lowest elevations
    if values(i) > max_elevation
        max_elevation = values(i);
        max_elevation_x = double(S.x);
        max_elevation_y = double(S.y);
    end
    if values(i) < min_elevation
        min_elevation = values(i);
        min_elevation_x = double(S.x);
        min_elevation_y = double(S.y);
    end
    
    sprintf('x: %f, y: %f, f(x,y): %f', values_x(i), values_y(i), values(i))
end

sprintf('Highest Elevation (x, y, f(x,y)): %f, %f, %f', max_elevation_x, max_elevation_y, max_elevation)
sprintf('Lowest Elevation (x, y, f(x,y)): %f, %f, %f', min_elevation_x, min_elevation_y, min_elevation)

% take second partial derivatives
A = diff(dzdx, x);
B = diff(dzdx, y);
C = diff(dzdy, y);
A_func = matlabFunction(A);
B_func = matlabFunction(B);
C_func = matlabFunction(C);

% iterate through for each critical point
for i = 1:length(values)
    
    % determine A, B, C, and D for second derivative test
    A_values = feval(A_func, values_x(i), values_y(i));
    B_values = feval(B_func, values_x(i), values_y(i));
    C_values = feval(C_func, values_x(i), values_y(i));
    D_values = B_values.^2 - A_values .* C_values;
    
    % determine if relative maximum or relative minimum
    if D_values < 0
        
        % determine if relative maximum
        if A_values < 0
            sprintf('%s', 'max')
        
        % determine if relative minimum
        else
            sprintf('%s', 'min')
        end
    
    % determine if saddle point
    elseif D_values > 0
        sprintf('%s', 'saddle')
    
    % determine if cannot be determined
    else
        sprintf('%s', 'nothing')
    end
end

% Part II:

% 2.a.

syms X Y Z T

% determine temperature at maximum and minimum elevations
T = -2.*Z.^2-exp(-0.1.*((0.1.*X-2)-(0.05.*Y-3).^2-(Z-1).^2))+10;
T_func = matlabFunction(T);
T_max_elevation = T_func(max_elevation_x, max_elevation_y, max_elevation);
T_min_elevation = T_func(min_elevation_x, min_elevation_y, min_elevation);
sprintf('Temp at Max Elevation: %f, Temp at Min Elevation: %f', T_max_elevation, T_min_elevation)

% 2.b.i.

% determine elevation at (2.2, 0.5) and calculate temperature at this point
hiker_elevation = feval(z_func, 2.2, 0.5);
T_hiker = T_func(2.2, 0.5, hiker_elevation);
sprintf('Temperature at (2.2, 0.5): %f', T_hiker)

% 2.b.ii.

% isotherm plot showing all contours within a +/- 1 tolerance range
[x,y] = meshgrid(x_range, y_range);
z = log(x.^4+1).*(4.*x.^4+(2.*y).^2).*exp(-0.5.*x.^2-(y-0.8).^2-3).*(sin(2.*x+0.05.*y.^4)+2.*cos(0.75.*y));
T = -2.*z.^2-exp(-0.1.*((0.1.*x-2)-(0.05.*y-3).^2-(z-1).^2))+10;
figure();
[c,h] = contour(x, y, T, ([T_hiker-1 T_hiker T_hiker+1]), 'red');
clabel(c,h);
title('Isotherms at Hiker''s Elevation');
xlabel('X [km]');
ylabel('Y [km]');

% 2.c.i.

% calculate the partial derivatives with respect to x and y at the point
hiker_x = feval(val_dzdx, 2.2, 0.5);
hiker_y = feval(val_dzdy, 2.2, 0.5);

% determine rate of change by dotting the hiker's gradient with north direction
hiker_direction = dot([hiker_x hiker_y], [0 1]);

% determine if hiker is ascending or descending based on sign of value
if hiker_direction >= 0
    sprintf('Hiker is ascending at %f', hiker_direction)
else
    sprintf('Hiker is descending at %f', hiker_direction)
end

% 2.c.ii.

% calculate partial derivative of T with respect to x, y, and z
T = -2.*Z.^2-exp(-0.1.*((0.1.*X-2)-(0.05.*Y-3).^2-(Z-1).^2))+10;
dTdx(X, Y, Z) = diff(T, X);
dTdy(X, Y, Z) = diff(T, Y);
dTdz(X, Y, Z) = diff(T, Z);

% evaluate each partial derivative at the point(2.2, 0.5) and at the hiker's elevation
T_gradient = [double(dTdx(2.2, 0.5, hiker_elevation)), double(dTdy(2.2, 0.5, hiker_elevation)), double(dTdz(2.2, 0.5, hiker_elevation))];

% determine the unit vector corresponding to its direction
T_direction = [0 1 hiker_direction];
T_unit_direction = T_direction./norm(T_direction);

% calculate the rate of change of the temperature
T_rate_of_change = dot(T_gradient, T_unit_direction);
sprintf('Rate of Change in Temperature: %f', T_rate_of_change)

% 2.d.i.

% determine rate of change by dotting the hiker's gradient with southwest direction
hiker_direction_2 = dot([hiker_x hiker_y], [-1/sqrt(2) -1/sqrt(2)]);

% determine if hiker is ascending or descending based on sign of value
if hiker_direction_2 >= 0
    sprintf('Hiker is ascending at %f', hiker_direction_2)
else
    sprintf('Hiker is descending at %f', hiker_direction_2)
end

% 2.d.ii.

T_direction_2 = [-1/sqrt(2) -1/sqrt(2) hiker_direction_2];
T_unit_direction_2 = T_direction_2./norm(T_direction_2);

T_rate_of_change_2 = dot(T_gradient, T_unit_direction_2);
sprintf('Rate of Change in Temperature: %f', T_rate_of_change_2)

% 2.e.

% define terrain and temperature functions
[x,y] = meshgrid(x_range, y_range);
z = log(x.^4+1).*(4.*x.^4+(2.*y).^2).*exp(-0.5.*x.^2-(y-0.8).^2-3).*(sin(2.*x+0.05.*y.^4)+2.*cos(0.75.*y));
T = -2.*z.^2-exp(-0.1.*((0.1.*x-2)-(0.05.*y-3).^2-(z-1).^2))+10;

color_map_max = 7.1877;
color_map_find_min = 0;
color_map_min = -22.8003;
color_map_min_x = 0;
color_map_min_y = 0;

% loop through entire domain of surface
for x0 = -5:0.1:5
    for y0 = -2:0.1:4
        
        % calculate elevation at every point
        z = log(x0.^4+1).*(4.*x0.^4+(2.*y0).^2).*exp(-0.5.*x0.^2-(y0-0.8).^2-3).*(sin(2.*x0+0.05.*y0.^4)+2.*cos(0.75.*y0));
        
        % calculate temperature at every point and corresponding elevation
        color_map_value = T_func(x0, y0, z);
        
        % extract lowest temperature (used later to verify Lagrange answer)
        if color_map_value < color_map_find_min
            color_map_find_min = color_map_value;
            color_map_min_x = x0;
            color_map_min_y = y0;
        end
    end
end

sprintf('Lowest Temperature (x, y, T): %f, %f, %f', color_map_min_x, color_map_min_y, color_map_find_min)

z = log(x.^4+1).*(4.*x.^4+(2.*y).^2).*exp(-0.5.*x.^2-(y-0.8).^2-3).*(sin(2.*x+0.05.*y.^4)+2.*cos(0.75.*y));

% plot terrain
figure();
surf(x,y,z,T);

% plot labels
title('Plot of Terrain with Temperature Color Map');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

% 2.f.

% plot terrain
figure();
surf(x,y,T);

% plot title and labels
title('Plot of Temperature');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

% 2.g.

syms lambda

% create Lagrange function of the form: L(x, y, z, lambda) = T(x, y, z) + (lambda * f(x, y) - z)
lagrange(X, Y, Z, lambda) = (-2.*Z.^2-exp(-0.1.*((0.1.*X-2)-(0.05.*Y-3).^2-(Z-1).^2))+10) + (lambda.*(log(X.^4+1).*(4.*X.^4+(2.*Y).^2).*exp(-0.5.*X.^2-(Y-0.8).^2-3).*(sin(2.*X+0.05.*Y.^4)+2.*cos(0.75.*Y)) - Z));

% calculate partial derivatives of the Lagrange function in terms of x, y, z, and lambda and set them all to 0
dLdx = diff(lagrange, X) == 0;
dLdy = diff(lagrange, Y) == 0;
dLdz = diff(lagrange, Z) == 0;
dLdLambda = diff(lagrange, lambda) == 0;

% initial guess
guess = [-2.3, 0.7, 3.5, 1];

% final result for lowest temperature
lowest_temp = vpasolve([dLdx, dLdy, dLdz, dLdLambda], [X, Y, Z, lambda], guess);
sprintf('Lowest Temperature (x, y, z): %f, %f, %f', lowest_temp.X, lowest_temp.Y, lowest_temp.Z)
