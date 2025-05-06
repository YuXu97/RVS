function [X,J,LYN] = bsi_wiener(data, model, par, X, J)
% Conditional bootstrap particle filter for SISO Wiener system with full noise
% parametrisation

input = (size(model.ss.B,2) > 0);
conditioning = (nargin > 3);

T = par.finalTime;

% Extract the model
A = model.ss.A; % Dynamic equation
B = model.ss.B;
C = model.ss.C; % Measurement equation
h = model.ss.h;
Q = model.ss.Q;

Ri = inv(model.ss.R); % Measurement precision
Qi = inv(Q);

% Data
if(input), u = data.u; end;
y = data.y;
LYN = zeros(1,T); % Total log-likelihood

% Initialize the particles
x_pred = zeros(model.nx, par.Np, T);
x_pred(:,:,1) = mvnrnd(model.ss.X1, model.ss.P1, par.Np)'; % Gaussian initial state
if(conditioning)
    x_pred(:, J(1), 1) = X(:,1); % Set the initial particle according to the conditioning
end

% Allocate memory for weights and ancestor indices
w = zeros(1,par.Np,T);
ancestors = zeros(1,par.Np,T);
extrc = 0; % Extrapolation counter

for(t = 1:T)
    % ===   PF TU - Predict new particles   ===
    % If we have our initial predicition, we don't have to predict again.
    if(t ~= 1)
        if(input)
            x_pred(:,:,t) = A*x_pred(:,ind,t-1) + repmat(B*u(t-1),1,par.Np) + mvnrnd(zeros(1,model.nx), Q, par.Np)';
        else
            x_pred(:,:,t) = A*x_pred(:,ind,t-1) + mvnrnd(zeros(1,model.nx), Q, par.Np)';
        end
        
        if(conditioning)
            x_pred(:, J(t),t) = X(:,t); % Set the t:th particle according to the conditioning
        end
    end

    % ===   PF MU - Compute importance weights and resample  ===
    [pred, extr] = h(C*x_pred(:,:,t));
    extrc = extrc + extr;
    innovation = y(t) - pred;
    % Compute weights - using the log-trick
    lq = loggausspdf2(innovation, Ri);    
    const = max(lq); % Subtract the maximum value
    q = exp(lq-const);
    if(const == -Inf), error('PF: Weights = 0!'); end
    LYN(t) = const + log(sum(q)) - log(par.Np);
    q = q/sum(q);
    w(:,:,t) = q; % Save the weights
    
    ind = resampling(q, par.resampling);
    ind = ind(randperm(par.Np));
    if(conditioning && t < T)
        % The J(t+1):th particle at t+1 originates from A_t^{J(t+1)) -
        % sample this using backward simulation
        if(input)
            p = gausspdf(coladd(-A*x_pred(:,:,t), X(:,t+1)-B*u(t)), Qi);
        else
            p = gausspdf(coladd(-A*x_pred(:,:,t), X(:,t+1)), Qi);
        end
            
        w_i = w(:,:,t).*p;
        w_i = w_i/sum(w_i);
        ind(J(t+1)) = catrnd(w_i);
    end
    
    % Save the ancestral line - x_t(i) originates from x_{t-1}(ancestors_t(i))
    ancestors(:,:,t) = ind;
end

% BACKWARD SIMULATION

% Initialize the smoothed particles and weights
X = zeros(model.nx, T);  
J = zeros(1, T);            % Index list for the sampled backward trajectory

% Initialise the backward trajectory by sampling from the FF particles at time T
J(T) = catrnd(w(:,:,T));
X(:,T) = x_pred(:,J(T),T);

for(t = (T-1) : (-1) : 1)
    xi_t = x_pred(:,:,t);    % x^i_t  - All particles from the forward pass at time t    
    x_t1 = X(:,t+1);    % x*_t+1 - Backward trajectory at time t+1
        
    % Linear Gaussian model
    if(input)
        f_xi_t = A*xi_t + repmat(B*u(t),1,par.Np);
    else
        f_xi_t = A*xi_t;
    end
        
    % Find the new weights (P(x*_t = xi_t | x*_t1))
    % Eq. (5) in Godsill et al. (2004)
    p = gausspdf(coladd(-f_xi_t, x_t1), Qi);
    w_i = w(:,:,t).*p;
    w_i = w_i/sum(w_i);
    
    J(t) = catrnd(w_i);
    X(:,t) = x_pred(:,J(t),t);
end
