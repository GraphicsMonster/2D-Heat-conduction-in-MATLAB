clc;
% Define parameters
Lx = 2; Ly = 1; % Define length only in x and y dimensions
Nx = 50; Ny = 25; % Adjust grid size for the rectangular rod
dx = Lx / (Nx - 1); dy = Ly / (Ny - 1); % Adjust grid spacing
dt = 0.01; T_final = 5;
alpha = 0.02; % diffusivity constant

% Initialize temperature field
T = zeros(Nx, Ny);

% Set boundary conditions
T(1, :) = 100; % Left boundary condition
T(end, :) = 100; % Right boundary condition

% Initialize storage for temperature at each time step
num_steps = ceil(T_final / dt);
T_history = zeros(Nx, Ny, num_steps);

% Time-stepping loop
t = 0;
step = 1;
while t < T_final
    % Store temperature at current time step
    T_history(:,:,step) = T;
    
    % Compute temperature at the next time step using finite differences
    T_new = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            T_new(i,j) = T(i,j) + alpha*dt*(...
                (T(i+1,j) - 2*T(i,j) + T(i-1,j))/dx^2 + ...
                (T(i,j+1) - 2*T(i,j) + T(i,j-1))/dy^2);
        end
    end
    T = T_new;
    
    % Update time and step
    t = t + dt;
    step = step + 1;
end

% Visualize the temporal evolution of the temperature field
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
figure;
for i = 1:num_steps
    % Plot the temperature distribution
    contourf(X, Y, T_history(:,:,i)', 'EdgeColor', 'none'); % Transpose T_history to match X and Y
    
    % Set aspect ratio to make it look more like a rod
    axis equal; % Set equal aspect ratio
    
    % Set limits of the axes to focus on the rod
    xlim([0, Lx]);
    ylim([0, Ly]);
    
    % Label axes and add colorbar
    xlabel('x');
    ylabel('y');
    colorbar;
    
    % Title with current time
    title(['Temperature Distribution at Time = ', num2str((i-1)*dt), 'seconds']);
    
    drawnow;
    pause(0.1); % Adjust the pause time as needed
end

% Extract temperature data for the last time step
T_surface = squeeze(T_history(:, :, end));

% Creating a mesh for the surface plot
[X_surf, Y_surf] = meshgrid(linspace(0, Lx, size(T_surface, 2)), linspace(0, Ly, size(T_surface, 1)));

% Creating the surface plot with colors representing temperature
figure;
surf(X_surf, Y_surf, T_surface, 'EdgeColor', 'none');
colorbar;

% Setting labels and title
xlabel('x');
ylabel('y');
zlabel('Temperature');
title('Temperature Distribution at the End of Simulation');
