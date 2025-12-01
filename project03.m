%% MAE 623 - CFD I 
% Project 03
% Jorge Luis Mares Zamora
% Due 12/01/2025

clear
clc 
close all

%% Input parameters

nu = 0.01; % kinematic viscosity
rho = 1; % density
Re = 100; % Reynolds number

L = 1; % Length of the square
N = 50; 

%% Implementation? 
[u, v, p] = incNavStokes(nu, rho, Re, L, N, N, 0.1);

%% Helper functions

function [u, v, P] = incNavStokes(nu, rho, Re, L, Nx, Ny, maxco, tol)
    % allowing users to NOT input max courant no. if they dont care??
    if nargin == 6
        maxco = 0.5;
        tol = 10 ^ (-6); 
    elseif nargin == 7
        tol = 10 ^ (-6); 
    end

    % initialize variables
    dx = L / (Nx - 1); 
    dy = L / (Ny - 1); 
    time = 0; 

    % Calculating dt assuming the max velocity is = 1 (since it is the velocity at the lid!)
    dt = maxco * dy * 1;    

    % Creating grid for u, v, and P
    % zero velocities as the initial condition
    u = zeros(Ny, Nx); 
    v = zeros(Ny, Nx); 
    ustar = zeros(Ny, Nx); 
    vstar = zeros(Ny, Nx); 
    P = zeros(size(u)); % Initial guess for P for the gauss seidel iterations!
    projector_algorithm = true; 

    % Iterations of the projection method algorithm. 
    while projector_algorithm 
        time = time + dt 

        % Step 01: Find u* and v* using (4) and (5)
        % --------------------------------------------------------
        % Important to note that u* and v* are only used for the calculation 
        % of the pressure, P, by using Poisson's Equation. This equation is
        % only used in the interior nodes! so u* and v* need only be calculated 
        % in the interior nodes of the domain and the boundary values can be 
        % safely ignored for now! :D

        for i = 2:(Nx-1) % going through all of the interior nodes!
            for j = 2:(Ny-1)
                % first find u*
                %convx = ((u(j, i+1))^2 - (u(j, i-1))^2)/(2*dx) + (u(j+1,i) * v(j+1,i) - u(j-1,i) * v(j-1, i))/(2*dy);  % convection term in x
                convx = u(j,i) * (u(j, i+1) - u(j,i-1))/(2 * dx) + v(j,i) * (u(j+1,i) - u(j-1,i))/(2*dy); 
                diffx = nu * ( (u(j,i+1) + u(j, i-1) - 2*u(j,i))/(dx^2) + (u(j+1, i) + u(j-1,i) - 2*u(j,i))/(dy^2)); % diffusion term in x
                ustar(j, i) = (- convx + diffx) * dt + u(j,i); 

                % then; we find v* 
                %convy = (u(j, i+1) * v(j, i+1) - u(j, i-1) * v(j, i-1))/(2*dx) + ((v(j+1, i))^2 - (v(j-1, i))^2)/(2*dy); 
                convy = u(j,i) * (v(j,i+1) - v(j,i-1))/(2*dx) + v(j,i) * (v(j+1, i) - v(j-1, i))/(2*dy); 
                diffy = nu * ( (v(j,i+1) + v(j, i-1) - 2*v(j,i))/(dx^2) + (v(j+1, i) + v(j-1, i) - 2*v(j, i))/(dy^2)); 
                vstar(j, i) = (-convy + diffy) * dt + v(j,i); 
            end 
        end

        % Step 02: Find P by solving Poisson's Equation (8)
        % ---------------------------------------------------------------------
        % P is not given as an IC, so a system has to be solved for this step
        % in this case, it was chosen to use Gauss Seidel with an initial guess
        % of P = 0 everywhere in the domain

        gauss = true; 
        gauss_iteration = 0; 

        while gauss
            P_old = P; 
            gauss_iteration = gauss_iteration + 1;

            % updating all of the interior nodes
            for i = 2:(Nx-1)
                for j = 2:(Ny-1) 
                    num1 = (P(j, i+1) + P(j, i-1))/(dx^2) + (P(j+1,i) + P(j-1, i))/(dy^2); 
                    num2 = (ustar(j, i+1) - ustar(j, i))/dx + (vstar(j+1, i) - vstar(j, i))/dy; 
                    den = -2/(dx^2) - 2/(dy^2); 
                    P(j, i) = (-num1 + (rho /dt) * num2) / den; 
                end 
            end

            % updating boundary values using linear extrapolation
            % error('worry about the interior corner nodes??? and corner ones as well lol')
            for i = 1:(Nx)
                % SOUTH boundary
                P(1, i) = 2 * P(2, i) - P(3, i); 
                % NORTH boundary
                P(Ny, i) = 2 * P(Ny-1, i) - P(Ny-2, i); 
            end
            for j = 1:(Ny)
                % EAST boundary
                P(j, 1) = 2 * P(j, 2) - P(j, 3); 
                % WEST boundary 
                P(j, Nx) = 2 * P(j, Nx-1) - P(j, Nx-2); 
            end
            P(2,2) = 0; % Fixing a node

            % boundary conditions on corner nodes! (8 nodes adjacent to corners)
            % here we apply 0th order extrapolation
            % SOUTH Boundary
            P(1, 2) = P(2, 2); 
            P(1, Nx-1) = P(2, Nx-1); 
            % EAST Boundary
            P(2, Nx) = P(2, Nx-1); 
            P(Ny-1, Nx) = P(Ny-1, Nx-1); 
            % WEST Boundary
            P(2, 1) = P(2, 2); 
            P(Ny-1, 1) = P(Ny-1, 2); 
            % NORTH Boundary
            P(Ny, 2) = P(Ny-1, 2); 
            P(Ny, Nx-1) = P(Ny-1, Nx-1); 

            % LASTLY, we can do the 4 corner nodes! (just extending value from corner)
            P(1,1) = P(2,2); % bottom left 
            P(1, Nx) = P(2, Nx-1); % bottom right 
            P(Ny, Nx) = P(Ny-1, Nx-1); % top right
            P(Ny, 1) = P(Ny-1, 2); % top left

            % Convergence criteria
            MAE = max(max(abs((P - P_old)))); % Magnitude approximate error
            if MAE < 10 ^ (-3) 
                break
            end
        end

        % STEP 03: Find U^(n+1) and V^(n+1)
        u_new = zeros(size(u)); 
        v_new = zeros(size(v)); 
        
        % Inner nodes first
        for i = 2:(Nx-1)
            for j = 2:(Ny-1)
                u_new(j,i) = - (dt / (rho * dx)) * (P(j,i) - P(j, i-1)) + ustar(j, i); 
                v_new(j,i) = - (dt / (rho * dy)) * (P(j,i) - P(j-1, i)) + vstar(j, i); 
            end 
        end
        % boundary nodes
        for i = 1:Nx
            % SOUTH boundary
            u_new(1, i) = 0; 
            v_new(1, i) = 0; 
            % NORTH boundary
            u_new(Ny, i) = 0; 
            v_new(Ny, i) = 0; 
        end
        for j = 1:Ny
            % EAST boundary
            u_new(j, 1) = 0; 
            v_new(j, 1) = 1; 
            % WEST boundary
            u_new(j, Nx) = 0; 
            v_new(j, Nx) = 0; 
        end

        % Calculating max courant numbers
        Co_x = max(max(u)) * dt / dx; 
        Co_y = max(max(v)) * dt / dy; 
        if max(Co_x, Co_y) > 0.5 
            fprintf('The max. Courant number is: Co_x = %.4f and Co_y = %.4f\n', Co_x, Co_y)
            disp('Warning! Courant number too high!!!')
        end
    
        MAE_x = max(max(abs(u_new - u))); 
        MAE_y = max(max(abs(v_new - v)));
        % Checking for steady state!
        if MAE_x < tol && MAE_y < tol
            disp('Converged to steady state!')
            break
        else
            u = u_new; 
            v = v_new; 
        end
    end
end