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

% Grid sizes 
% (2x each time for grid refinement study!! ^_^)
N_coarse = 20; 
N_medium = 100; 
N_fine = 200; 

%% Report Results!!
[u_coarse, v_coarse, p_coarse] = incNavStokes(nu, rho, Re, L, N_coarse, N_coarse, 0.1);
%[u_medium, v_medium, p_medium] = incNavStokes(nu, rho, Re, L, N_medium, N_medium, 0.1); 
%[u_fine, v_fine, p_fine] = incNavStokes(nu, rho, Re, L, N_fine, N_fine, 0.1); 

%reportIndividualPlotting(u_coarse, v_coarse, p_coarse); 

% vector field plot
% fig = figure('Units','inches','Position',[1 1 5 4]);
% x = linspace(0,1, N_coarse); 
% y = linspace(0,1, N_coarse); 
% [X, Y] = meshgrid(x, y); 
% q = quiver(X, Y, u_coarse * 2, v_coarse * 2,'k', 'MarkerSize',6);
% q.ShowArrowHead = 'off';
% q.Marker = '.';
% xlim([0 1])
% ylim([0 1])
% axis tight
% xlabel('x')
% ylabel('y')
% title('Velocity vector field')
% set(gca,'FontSize',12)
% set(fig,'Color','w') 
% exportgraphics(fig,'velocity_vector_field.png','Resolution',600)

% %% Grid Refinement study!!
% gridRefinement(u_coarse, u_medium, u_fine, 'u'); 
% gridRefinement(v_coarse, v_medium, v_fine, 'v'); 
% gridRefinement(p_coarse, p_medium, p_fine, 'p'); 

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
    target = 0; 

    % Calculating dt assuming the max velocity is = 1 (since it is the velocity at the lid!)
    dt = maxco * dy * 1;    

    % Creating grid for u, v, and P
    % zero velocities as the initial condition
    u = zeros(Ny, Nx); 
    v = zeros(Ny, Nx); 
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

        ustar = zeros(size(u)); 
        vstar = zeros(size(v)); 

        for i = 2:(Nx-1) % going through all of the interior nodes!
            for j = 2:(Ny-1)
                % first find u*
                convx = ((u(j, i+1))^2 - (u(j, i-1))^2)/(2*dx) + (u(j+1,i) * v(j+1,i) - u(j-1,i) * v(j-1, i))/(2*dy);  % convection term in x
                %convx = u(j,i) * (u(j, i+1) - u(j,i-1))/(2 * dx) + v(j,i) * (u(j+1,i) - u(j-1,i))/(2*dy); 
                diffx = nu * ( (u(j,i+1) + u(j, i-1) - 2*u(j,i))/(dx^2) + (u(j+1, i) + u(j-1,i) - 2*u(j,i))/(dy^2)); % diffusion term in x
                ustar(j, i) = (- convx + diffx) * dt + u(j,i); 

                % then; we find v* 
                convy = (u(j, i+1) * v(j, i+1) - u(j, i-1) * v(j, i-1))/(2*dx) + ((v(j+1, i))^2 - (v(j-1, i))^2)/(2*dy); 
                %convy = u(j,i) * (v(j,i+1) - v(j,i-1))/(2*dx) + v(j,i) * (v(j+1, i) - v(j-1, i))/(2*dy); 
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
        gauss_tol = 10 ^ (-3);

        while gauss
            P_old = P; 
            gauss_iteration = gauss_iteration + 1; 

            % updating all of the interior nodes
            for i = 2:(Nx-1)
                for j = 2:(Ny-1) 
                    num1 = (P(j, i+1) + P(j, i-1))/(dx^2) + (P(j+1,i) + P(j-1, i))/(dy^2); 
                    num2 = (ustar(j, i+1) - ustar(j, i))/dx + (vstar(j+1, i) - vstar(j, i))/dy; 
                    den = -2/(dx^2) - 2/(dy^2); 
                    P(j, i) = (-num1 + (rho/dt) * num2)/den; 
                end 
            end
            P(2,2) = 0; 

            % Updating boundaries! (corners are IGNORED)
            % North and South boundaries first
            for i = 2:(Nx-1)
                if i == 2 || i == (Nx-1) % Corner values
                    % SOUTH boundary
                    P(1, i) = P(2, i);
                    % NORTH boundary
                    P(Ny, i) = P(Ny-1, i); 
                end
                
                % Rest of SOUTH boundary
                %P(1, i) = 2 * P(2, i) - P(3, i); 
                P(1, i) = P(2,i); 
                % Rest of NORTH boundary
                %P(Ny, i) = 2 * P(Ny-1, i) - P(Ny-2, i); 
                P(Ny, i) = P(Ny-1, i); 
            end
            % Now East and West boundaries!
            for j = 2:(Ny-1)
                if j == 2 || j == (Ny-1) % corner values first
                    % EAST boundary
                    P(j, 1) = P(j, 2); 
                    % WEST boundary
                    P(j, Nx) = P(j, Nx-1); 
                end

                % EAST boundary again!
                %P(j, 1) = 2 * P(j, 2) - P(j, 3); 
                P(j, 1) = P(j, 2); 
                % WEST boundary again!
                %P(j, Nx) = 2 * P(j, Nx-1) - P(j, Nx-2); 
                P(j, Nx) = P(j, Nx-1); 
            end

           % Convergence criteria
           e = abs(P - P_old); % magnitude approximate error
            if max(max(e)) < gauss_tol
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
        Co_x = max(max(u_new)) * dt / dx; 
        Co_y = max(max(v_new)) * dt / dy; 
        if max(Co_x, Co_y) > 0.5 
            fprintf('The max. Courant number is: Co_x = %.4f and Co_y = %.4f\n', Co_x, Co_y)
            disp('Warning! Courant number too high!!!')
        end
    
        MAE_x = max(max(abs(u_new - u))); 
        MAE_y = max(max(abs(v_new - v)));
        % Checking for steady state!
        if MAE_x < tol && MAE_y < tol
            disp('Converged to steady state!')

            % before breaking we shift p! remember we only care about the 
            % shape of p and not about the absolute pressure of the fluid!
            shift = P(Ny/2,Nx/2) - target; 
            P = P - shift; 
            break
        else
            u = u_new; 
            v = v_new; 
        end
    end
end

function reportIndividualPlotting(u, v, p)
    n = size(u, 1); 

    % Create arrays for x and y
    y = 0:(1/(n-1)):1; 
    x = y; 

    uvals = u(n/2, :); 
    vvals = v(n/2, :); 
    pvals = p(n/2, :); 

    % Plot u, v, p vs x @ y = 0.5
    figure()
    subplot(3,1,1)
    plot(x, uvals, 'm')
    title('u vs x @ y = 0.5')
    xlabel('x values')
    ylabel('u values')
    grid on 

    subplot(3,1,2)
    sgtitle(['Plots for grid resolution of Nx = Ny = ', num2str(n)])
    plot(x, vvals, 'm')
    title('v vs x @ y = 0.5')
    xlabel('x values')
    ylabel('v values')
    grid on 
 
    subplot(3,1,3)
    plot(x, pvals, 'm')
    title('p vs x @ y = 0.5')
    xlabel('x values')
    ylabel('p values')
    grid on 

    % Plot u, v, p, vs y @ x = 0.5
    uvals = u(:, n/2); 
    vvals = v(:, n/2); 
    pvals = p(:, n/2); 

    figure()
    subplot(3,1,1)
    sgtitle(['Plots for grid resolution of Nx = Ny = ', num2str(n)])

    plot(uvals, y, 'm')
    title('u vs y @ x = 0.5')
    xlabel('u values')
    ylabel('y values')
    grid on

    subplot(3,1,2)
    plot(vvals, y, 'm')
    title('v vs y @ x = 0.5')
    xlabel('v values')
    ylabel('y values')
    grid on
    
    subplot(3,1,3)
    plot(pvals, y, 'm')
    title('p vs y @ x = 0.5')
    xlabel('p values')
    ylabel('y values')
    grid on
end

function gridRefinement(phi1, phi2, phi3, field)
    n1 = size(phi1, 1); 
    n2 = size(phi2, 1); 
    n3 = size(phi3, 1); 

    x1 = 0:(1/(n1-1)):1; 
    x2 = 0:(1/(n2-1)):1; 
    x3 = 0:(1/(n3-1)):1; 

    y1 = x1; 
    y2 = x2; 
    y3 = x3; 


    % first plotting against x @ y = 0.5
    phi1vals = phi1(n1/2, :); 
    phi2vals = phi2(n2/2, :); 
    phi3vals = phi3(n3/2, :); 

    figure()
    sgtitle(['Grid Refinement Study on ', field, ' values.'])
    subplot(2,1,1)
    plot(x1, phi1vals, 'k', x2, phi2vals, 'b', x3, phi3vals, 'm')
    title([field, ' vs x @y = 0.5 for different grid resolutions'])
    xlabel('x values')
    ylabel([field, ' values'])
    legend(['Nx = Ny = ', num2str(n1)], ['Nx = Ny = ', num2str(n2)], ['Nx = Ny = ', num2str(n3)])
    grid on

    % Then plotting against y @ x = 0.5
    phi1vals = phi1(:, n1/2); 
    phi2vals = phi2(:, n2/2); 
    phi3vals = phi3(:, n3/2); 

    subplot(2,1,2)
    plot(phi1vals, y1, 'k', phi2vals, y2, 'b', phi3vals, y3, 'm')
    title([field, ' vs y @ x = 0.5 for different grid resolutions'])
    ylabel('y values')
    xlabel([field, ' values'])
    legend(['Nx = Ny = ', num2str(n1)], ['Nx = Ny = ', num2str(n2)], ['Nx = Ny = ', num2str(n3)])
    grid on
end