clc; clear; close all;

%% Step 1: Create a Uniform 3D Grid Using meshgrid
N = 10;  % Grid resolution (NxNxN points)
L = 10;  % Grid extends from -L to L in all directions
[X, Y, Z] = meshgrid(linspace(-L, L, N), linspace(-L, L, N), linspace(-L, L, N));

%% Step 2: Convert to Spherical Coordinates
Rs = 3; % Schwarzschild radius
R = max(sqrt(X.^2 + Y.^2 + Z.^2), Rs + 1e-6);  % Avoid complex numbers by enforcing R > Rs
Theta = acos(Z ./ R);  % Compute polar angle
Phi = atan2(Y, X);  % Compute azimuthal angle

%% Step 3: Apply Schwarzschild Stretching
R_stretched = (R + Rs * log((sqrt(R) + sqrt(R - Rs)) / sqrt(Rs))); 
R_new = (R.^2) ./ R_stretched;  % Apply radial compression

%% Step 4: Convert Back to Cartesian Coordinates
X_new = R_new .* sin(Theta) .* cos(Phi);
Y_new = R_new .* sin(Theta) .* sin(Phi);
Z_new = R_new .* cos(Theta);

%% Step 5: Visualize the Original vs Stretched Grid
figure;
subplot(1,2,1);
scatter3(X(:), Y(:), Z(:), 10, 'filled'); % Original Grid
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Original 3D Grid');
grid on; axis equal;

subplot(1,2,2);
scatter3(X_new(:), Y_new(:), Z_new(:), 10, 'filled'); % Schwarzschild-Stretched Grid
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Schwarzschild-Stretched 3D Grid');
grid on; axis equal;

%% Step 6: Visualize as a Wireframe to Show Connections
figure;
hold on;
for i = 1:N
    for j = 1:N
        % Z-direction (k varies)
        plot3( squeeze(X_new(i, j, :)), squeeze(Y_new(i, j, :)), squeeze(Z_new(i, j, :)), 'b'); 

        % Y-direction (j varies)
        plot3( squeeze(X_new(i, :, j)), squeeze(Y_new(i, :, j)), squeeze(Z_new(i, :, j)), 'r');

        % X-direction (i varies)
        plot3( squeeze(X_new(:, i, j)), squeeze(Y_new(:, i, j)), squeeze(Z_new(:, i, j)), 'g');
    end
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Wireframe of Schwarzschild-Stretched Grid');
grid on; axis equal;
view(3);            % Ensures a 3D view
rotate3d on;        % Enables 3D rotation with the mouse
hold off;
