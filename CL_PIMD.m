classdef CL_PIMD
    % continuous path path integral molecular dynamics 
    %   Underdamped Langevin dynamics is evolved numerically
    %   Sample the normal modes rather than ring polymer beads
    %   All vectors are understood as column vectors
    
    properties
        %% method parameters
        %    N      [1] number of normal modes (odd)
        %    D      [1] number of grid points
        %    beta_D [1] beta / D
        N; D; beta_D;
    end
    
    methods
        %% initialize method parameters
        
        function pimd = CL_PIMD(N, D)
            pimd.N = N;
            pimd.D = D;
            pimd.beta_D = para.beta / pimd.D;
        end
        
        % initialize xi [N] and eta [N]
        % Output: xi  [N] normal modes coordinates
        %         eta [N] fictional velocity variables
        function [xi, eta] = initialize(pimd)
            xi  = zeros(pimd.N,1); % [N]
            eta = zeros(pimd.N,1); % [N]
            for k = 1:pimd.N
                w_k = para.eigen_w(k);
                xi(k)  = randn() / sqrt(w_k^2+1);
                eta(k) = randn() / sqrt(w_k^2+1);
            end
        end
 
        %% truncated underdamped Langevin dynamics
        %    xi  [N] normal mode coordinates
        %    eta [N] fictional velocities
        %    dt  [1] time step for numerical simulation of Langevin dynamics
        
        % compute observable integral
        % Input:  xi [N] normal modes coordinates
        % Output: o  [5] observable integral
        function o = compute_o(pimd, xi)
            % compute the continuous loop x(tau)
            tau = (0:(pimd.D-1))' * pimd.beta_D; % [D] grid points in [0,beta]
            x = zeros(pimd.D,1); % [D] continuous loop on grid points
            for k = 1:pimd.N
                x = x + xi(k) * para.eigen_c(k, tau); % [D]
            end
            % compute observable integral
            o = 0;
            for j = 1:pimd.D
                o = o + para.O(x(j));
            end
            o = o / pimd.D;
        end
        
        % compute drift force
        % Input:  xi [N] normal modes coordinates
        % Output: f  [N] drift force on N normal modes
        function f = compute_f(pimd, xi)
            % compute the continuous loop x(tau)
            tau = (0:(pimd.D-1))' * pimd.beta_D; % [D] grid points in [0,beta]
            x = zeros(pimd.D,1); % [D] continuous loop on grid points
            for k = 1:pimd.N
                x = x + xi(k) * para.eigen_c(k, tau); % [D]
            end
            % compute grad_V on grid points 
            grad_v = zeros(pimd.D,1); % [D] grad_V on grid points
            for j = 1:pimd.D
                grad_v(j) = para.grad_V(x(j));
            end
            % compute gradient integral
            f = zeros(pimd.N,1);
            for k = 1:pimd.N
                w_k = para.eigen_w(k);
                f(k) = - pimd.beta_D * dot(para.eigen_c(k, tau), grad_v) / (w_k^2+1); 
            end
        end
        
        % Caley integrator for part A
        % free ring polymer evolution
        % Input:  xi  [N+1] normal modes coordinates
        %         eta [N+1] fictional velocity variables
        %         dt  [1] time step
        % Output: xi  [N+1] normal modes coordinates (updated)
        %         eta [N+1] fictional velocity variables (updated)
        function [xi, eta] = simulation_C(pimd, xi, eta, dt)
            A = [0 1;-1 0]; % [2*2]
            cay_A = (eye(2)-dt/2*A) \ (eye(2)+dt/2*A); % Cayley integrator
            for k = 1:pimd.N
                xi_eta_k = [xi(k); eta(k)]; % [2] initial xi(k) and eta(k)
                xi_eta_k = cay_A * xi_eta_k; % Cayley integrator
                xi(k) = xi_eta_k(1); eta(k) = xi_eta_k(2); % update xi(k) and eta(k)
            end
        end
        
        % numerical integrator for part B    
        % Input:  xi  [N] normal modes coordinates
        %         eta [N] fictional velocity variables
        %         dt  [1] time step
        % Output: eta [N] fictional velocity variables (updated)
        function eta = simulation_B(pimd, xi, eta, dt)
            f = pimd.compute_f(xi);
            eta = eta + dt*f; % update xi in major modes
        end
        
        % numerical integrator for part O
        % Input:  eta [N] fictional velocity variables
        %         dt  [1] time step
        % Output: eta [N] fictional velocity variables (updated)
        function eta = simulation_O(pimd, eta, dt)
            gamma = 1; % damping rate
            c = exp(-gamma*dt); % contraction rate
            for k = 1:pimd.N
                w_k = para.eigen_w(k);
                eta(k) = c*eta(k) + sqrt((1-c^2)/(w_k^2+1))*randn(); % update eta in all modes
            end
        end
        
        % ---- BAOAB numerical integrator ---- %
        % Input:  xi  [N+1] normal modes coordinates
        %         eta [N+1] fictional velocity variables
        %         dt  [1] time step
        % Output: xi  [N+1] normal modes coordinates (updated)
        %         eta [N+1] fictional velocity variables (updated)
        function [xi, eta] = BCOCB(pimd, xi, eta, dt)
            eta = pimd.simulation_B(xi, eta, dt/2);
            [xi, eta] = pimd.simulation_C(xi, eta, dt/2);
            eta = pimd.simulation_O(eta, dt);
            [xi, eta] = pimd.simulation_C(xi, eta, dt/2);
            eta = pimd.simulation_B(xi, eta, dt/2);
        end
    end
end