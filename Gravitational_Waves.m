clc; clear; close all;

%% =============== Step 0: Grid and Simulation Setup ===============
% Grid size (keep it small to avoid slow plotting)
N = 30;                   
L = 20; 
[X, Y, Z] = meshgrid(linspace(-L, L, N), linspace(-L, L, N), linspace(-L, L, N));

% Spherical coordinates of the grid (for later TT transformation)
R     = sqrt(X.^2 + Y.^2 + Z.^2);
Theta = acos(Z ./ (R + eps));
Phi   = atan2(Y, X);

% Simulation time parameters
dt      = 0.05;         % time step
T_total = 200;           % total simulation time
time    = 0:dt:T_total;
numSteps = length(time);
c       = 1;            % speed of light in toy units

% Total number of grid points
numPoints = numel(X);

%% =============== Step 1: N-Body Setup (Binary Stars + Planet) ===============
% Binary star parameters:
m1 = 1;          % mass of star 1
m2 = 1;          % mass of star 2
M_tot = m1 + m2; % total mass
mu    = (m1*m2) / M_tot; % reduced mass

% Initial separation between the stars (distance between them)
r_sep = 10;
% Place stars along the x-axis (center-of-mass at origin)
x1 =  r_sep/2;  y1 = 0;  z1 = 0;
x2 = -r_sep/2;  y2 = 0;  z2 = 0;

% For a circular orbit, relative speed is: v_rel = sqrt(M_tot/r_sep)
v_rel = sqrt(M_tot / r_sep);
% Set initial velocities so that the binary orbits in the xy-plane
% (Using symplectic Euler, we update velocities then positions)
v1x = 0;           v1y =  (m2/M_tot) * v_rel;    v1z = 0;
v2x = 0;           v2y = -(m1/M_tot) * v_rel;    v2z = 0;

% Planet (test particle) initial conditions:
xp  = 1;  yp  = 0;  zp = 0;
vxp = 0;  vyp = 0;  vzp = 0;
planetMass = 0.01;   % small mass for the planet

% Radiation reaction parameters:
% We approximate energy loss using a drag acceleration.
% Damping coefficient D (from Peters formula scaling)
% D = (32/5)*(mu^2*M_tot^2) / (r12^4), computed later with r12 = separation between stars

%% =============== Step 2: Preallocate Arrays for Saving Data ===============
% Save binary positions
x1save = zeros(1, numSteps);  y1save = zeros(1, numSteps);  z1save = zeros(1, numSteps);
x2save = zeros(1, numSteps);  y2save = zeros(1, numSteps);  z2save = zeros(1, numSteps);

% Save planet positions (background, before TT transformation)
xpSave = zeros(1, numSteps);  ypSave = zeros(1, numSteps);  zpSave = zeros(1, numSteps);

% Save instantaneous orbital frequency (for GW calculations)
omegaSave = zeros(1, numSteps);

%% =============== Step 3: Main Integration Loop (Symplectic Euler) ===============
% In symplectic Euler we update velocities first then use the new velocity to update positions.

for k = 1:numSteps
    t = time(k);
    
    %% (a) --- Compute Gravitational Accelerations for the Binary Stars ---
    % Distance vector between stars:
    r12_vec = [x1 - x2, y1 - y2, z1 - z2];
    r12 = norm(r12_vec) + eps;
    
    % Newtonian acceleration on star 1 due to star 2 (G = 1)
    a1_star2 = - m2 * r12_vec / r12^3;
    % Acceleration on star 1 due to the planet:
    r1p_vec = [x1 - xp, y1 - yp, z1 - zp];
    r1p = norm(r1p_vec) + eps;
    a1_planet = - planetMass * r1p_vec / r1p^3;
    a1_newt = a1_star2 + a1_planet;
    
    % For star 2:
    a2_star1 = - m1 * (-r12_vec) / r12^3;  % note: -r12_vec = r21_vec
    r2p_vec = [x2 - xp, y2 - yp, z2 - zp];
    r2p = norm(r2p_vec) + eps;
    a2_planet = - planetMass * r2p_vec / r2p^3;
    a2_newt = a2_star1 + a2_planet;
    
    %% (b) --- Radiation Reaction (Damping) on the Stars ---
    % Damping coefficient based on instantaneous separation:
    D = (32/5) * (mu^2 * M_tot^2) / (r12^4);
    % Compute center-of-mass velocity of the binary:
    vcm = [(v1x*m1 + v2x*m2)/M_tot, (v1y*m1 + v2y*m2)/M_tot, (v1z*m1 + v2z*m2)/M_tot];
    % Radiation reaction acceleration (applied only to the orbital motion):
    a1_rad = - (D/m1) * ([v1x, v1y, v1z] - vcm);
    a2_rad = - (D/m2) * ([v2x, v2y, v2z] - vcm);
    
    % Total acceleration on each star:
    a1_total = a1_newt + a1_rad;
    a2_total = a2_newt + a2_rad;
    
    %% (c) --- Update Velocities and Positions (Symplectic Euler) ---
    % Update star 1:
    v1x = v1x + a1_total(1)*dt;
    v1y = v1y + a1_total(2)*dt;
    v1z = v1z + a1_total(3)*dt;
    x1 = x1 + v1x*dt;
    y1 = y1 + v1y*dt;
    z1 = z1 + v1z*dt;
    
    % Update star 2:
    v2x = v2x + a2_total(1)*dt;
    v2y = v2y + a2_total(2)*dt;
    v2z = v2z + a2_total(3)*dt;
    x2 = x2 + v2x*dt;
    y2 = y2 + v2y*dt;
    z2 = z2 + v2z*dt;
    
    % Compute gravitational accelerations for the planet (test particle):
    r_p1_vec = [xp - x1, yp - y1, zp - z1];
    r_p1 = norm(r_p1_vec) + eps;
    r_p2_vec = [xp - x2, yp - y2, zp - z2];
    r_p2 = norm(r_p2_vec) + eps;
    a_p1 = - m1 * r_p1_vec / r_p1^3;
    a_p2 = - m2 * r_p2_vec / r_p2^3;
    a_p_newt = a_p1 + a_p2;
    
    % Update planet velocity and position:
    vxp = vxp + a_p_newt(1)*dt;
    vyp = vyp + a_p_newt(2)*dt;
    vzp = vzp + a_p_newt(3)*dt;
    xp = xp + vxp*dt;
    yp = yp + vyp*dt;
    zp = zp + vzp*dt;
    
    %% (d) --- Save Positions and Compute Instantaneous Orbital Frequency ---
    x1save(k) = x1;  y1save(k) = y1;  z1save(k) = z1;
    x2save(k) = x2;  y2save(k) = y2;  z2save(k) = z2;
    xpSave(k) = xp;  ypSave(k) = yp;  zpSave(k) = zp;
    
    % Estimate orbital frequency omega from the binary's relative motion:
    v_rel = [v1x - v2x, v1y - v2y, v1z - v2z];
    L_vec = cross(r12_vec, v_rel);
    L_mag = norm(L_vec);
    omega_inst = L_mag / (mu * r12^2 + eps);
    omegaSave(k) = omega_inst;
end

%% =============== Step 4: Precompute TT Gauge Distortion for the Grid ===============
% Instead of a simple radial stretch, we apply the TT gauge transformation.
% For a wave propagating in the z direction, the TT gauge modifies only x and y.

A = 0.1;  % amplitude of gravitational waves

% Preallocate arrays for TT-transformed grid positions:
xGridTT = zeros(numPoints, numSteps); 
yGridTT = zeros(numPoints, numSteps);
zGridTT = zeros(numPoints, numSteps);

for k = 1:numSteps
    t = time(k);
    
    % Retrieve star positions at time step k:
    x1 = x1save(k);  y1 = y1save(k);  z1 = z1save(k);
    x2 = x2save(k);  y2 = y2save(k);  z2 = z2save(k);
    
    % Compute distances from each star to every grid point:
    d1 = sqrt((X - x1).^2 + (Y - y1).^2 + (Z - z1).^2);
    d2 = sqrt((X - x2).^2 + (Y - y2).^2 + (Z - z2).^2);
    
    % Retarded times:
    t1_ret = t - d1/c;
    t2_ret = t - d2/c;
    
    % Angle (theta) between star-to-grid vector and z-axis:
    theta1 = acos((Z - z1) ./ (d1 + eps));
    theta2 = acos((Z - z2) ./ (d2 + eps));
    
    % Quadrupole pattern factors:
    pattern1 = (3*cos(theta1).^2 - 1) / 2;
    pattern2 = (3*cos(theta2).^2 - 1) / 2;
    

    % Use precomputed orbital frequency:
    omega = omegaSave(k);
    
    % Compute polarization strains from each star:
    hplus1  = A * cos(2*omega*t1_ret) .* pattern1;
    hcross1 = A * sin(2*omega*t1_ret) .* pattern1;
    hplus2  = A * cos(2*omega*t2_ret) .* pattern2;
    hcross2 = A * sin(2*omega*t2_ret) .* pattern2;
    
    % Total strains (superposition):
    hplus_total = hplus1 + hplus2;
    hcross_total = hcross1 + hcross2;
    
    % Apply TT gauge transformation for each grid point:
    % For a point (x,y,z), the new coordinates are:
    %   x_new = x + 0.5*(hplus_total*x + hcross_total*y)
    %   y_new = y - 0.5*(hplus_total*y - hcross_total*x)
    %   z remains unchanged.
    X_new = X + 0.5 * (hplus_total .* X + hcross_total .* Y);
    Y_new = Y - 0.5 * (hplus_total .* Y - hcross_total .* X);
    Z_new = Z;
    
    % Save the TT-transformed grid positions:
    xGridTT(:, k) = X_new(:);
    yGridTT(:, k) = Y_new(:);
    zGridTT(:, k) = Z_new(:);
end

%% =============== Step 5: Precompute TT Gauge Distortion for the Planet ===============
% Apply the TT gauge transformation to the planet's background coordinates.
xpTT = zeros(1, numSteps);
ypTT = zeros(1, numSteps);
zpTT = zeros(1, numSteps);

for k = 1:numSteps
    t = time(k);
    
    % Planet's background position:
    xp_bg = xpSave(k);  yp_bg = ypSave(k);  zp_bg = zpSave(k);
    
    % For each star, compute the distance and retarded time from the planet:
    d1p = sqrt((xp_bg - x1save(k))^2 + (yp_bg - y1save(k))^2 + (zp_bg - z1save(k))^2);
    d2p = sqrt((xp_bg - x2save(k))^2 + (yp_bg - y2save(k))^2 + (zp_bg - z2save(k))^2);
    
    t1_ret = t - d1p/c;
    t2_ret = t - d2p/c;
    
    % Compute angles for the planet (with respect to each star):
    theta1p = acos((zp_bg - z1save(k)) / (d1p + eps));
    theta2p = acos((zp_bg - z2save(k)) / (d2p + eps));
    
    pattern1p = (3*cos(theta1p)^2 - 1)/2;
    pattern2p = (3*cos(theta2p)^2 - 1)/2;
    
    % Use the orbital frequency at step k:
    omega_k = omegaSave(k);
    
    % Compute strains at the planet's location from each star:
    hplus1 = A * cos(2*omega_k*t1_ret) * pattern1p;
    hcross1 = A * sin(2*omega_k*t1_ret) * pattern1p;
    hplus2 = A * cos(2*omega_k*t2_ret) * pattern2p;
    hcross2 = A * sin(2*omega_k*t2_ret) * pattern2p;
    
    hplus_total = hplus1 + hplus2;
    hcross_total = hcross1 + hcross2;
    
    % Apply the TT gauge transformation to the planet's coordinates:
    xp_new = xp_bg + 0.5*(hplus_total*xp_bg + hcross_total*yp_bg);
    yp_new = yp_bg - 0.5*(hplus_total*yp_bg - hcross_total*xp_bg);
    zp_new = zp_bg;  % unchanged
    
    xpTT(k) = xp_new;
    ypTT(k) = yp_new;
    zpTT(k) = zp_new;
end

%% =============== Step 6: Animate the Results Using the TT Distortion ===============
% Adjust simulation speed via simSpeed:
simSpeed = 0.00001;  % lower for faster animation, higher for slower

figure;
% Plot the TT-transformed grid (as a scatter3 plot)
hGrid = scatter3(xGridTT(:,1), yGridTT(:,1), zGridTT(:,1), 20, 'k', 'filled'); 
hold on;

% Plot the stars (using their background positions)
hStar1 = plot3(x1save(1), y1save(1), z1save(1), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);
hStar2 = plot3(x2save(1), y2save(1), z2save(1), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);

% Plot the planet (using the TT-transformed positions)
hPlanet = plot3(xpTT(1), ypTT(1), zpTT(1), 'go', 'MarkerFaceColor','g', 'MarkerSize',6);

axis([-L L -L L -L L]); 
axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Inspiraling Binary with TT Gauge Gravitational Waves Affecting the Planet');
rotate3d on;

for k = 1:numSteps
    % Update grid positions:
    set(hGrid, 'XData', xGridTT(:,k), 'YData', yGridTT(:,k), 'ZData', zGridTT(:,k));
    
    % Update star positions:
    set(hStar1, 'XData', x1save(k), 'YData', y1save(k), 'ZData', z1save(k));
    set(hStar2, 'XData', x2save(k), 'YData', y2save(k), 'ZData', z2save(k));
    
    % Update planet position:
    set(hPlanet, 'XData', xpTT(k), 'YData', ypTT(k), 'ZData', zpTT(k));
    
    drawnow;
    pause(simSpeed);
end
