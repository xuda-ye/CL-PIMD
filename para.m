classdef para
% parameters of the system to compute the quantum thermal average
% Hamiltonian operator: H = H^a + V^a
%   H^a is the quantum harmonic oscillator
%   V^a is the reduced potential function
%   a = 1 is assumed

    properties (Constant)
        beta = 2; % inverse temperature
        o_refer = [-0.2364364, -0.2932499]'; 
    end
    
    methods (Static)
        %% model parameters
        %    V^a(q)      [1] potential function         
        %    grad_V^a(q) [1] gradient of potential function
        %    O(q)        [1] observable function        
        
        % potential function V^a(q)
        function v = V(q)
            % Input:  q [1] position of particle
            % Output: v [1] potential function
            %v = q*cos(q);
             v = 4*sin(2*pi*q);
        end
        
        % gradient of potential function V^a(q)
        function grad_v = grad_V(q)
            % Input:  q [1] position of particle
            % Output: grad_v [1] gradient of potential function
            %grad_v = cos(q) - q*sin(q);
            grad_v = 8*pi*cos(2*pi*q);
        end
        
        % observable function
        function o = O(q)
            % Input:  q [1] position of particle
            % Output: o [1] observable value
            o = sin(pi/2*q)';
        end
        
        %% eigenvalues and eigenfunctions
        
        % calculate the square root of the eigenvalues
        % Input:  k [1] mode number
        % Output: w [1] eigenvalue w_k^2
        function w = eigen_w(k)
            if k == 1 % trivial eigenvalue
                w = 0;
                return
            end
            k = floor(k/2);
            w = 2*k*pi / para.beta;
        end
        
        % calculate the values of the eigenfunction on grid points
        % Input:  k   [1] mode number
        %         tau [-] grid points on [0,beta]
        % Output: c   [1] eigenfunction
        function c = eigen_c(k, tau)
            if k == 1 % trivial eigenvalue
                c = sqrt(1/para.beta) * ones(size(tau)); % constant eigenfunction
                return
            end
            c = mod(k,2); % c is odd or even
            k = floor(k/2);
            if c == 0 % sin
                c = sqrt(2/para.beta) * sin(2*k*pi*tau/para.beta);
            else % cos
                c = sqrt(2/para.beta) * cos(2*k*pi*tau/para.beta);
            end
        end
    end
end

