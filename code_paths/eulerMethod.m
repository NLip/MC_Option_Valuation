%% eulerMethod
% Implementation of Euler-Maruyama method for numerical integration of SODE dSt = a(t,St)dt + b(t,St)dWt
%
%   Input arguments:
%		a: function handle for drift term
%       b: function handle for diffusion term
%		X0: initial value
%		M: number of samples to produce in temporal dimension
%		T: final time
%		N: number of samples to produce
%		dt: timestep
%   Output arguments:
%       X: [N,M] matrix representing N samples of the process S (each M equidistant samples on [0,T])
function X = eulerMethod(a,b,X0,M,T,N,dt)
    X = zeros(N,M);
    X(:,1) = X0;
    t = 0;
    for k = 1:(M-1)
        Xk = X(:,k);
        for i = 1:((T/M)/dt) 
            dW = sqrt(dt)*randn(N,1);
            t = t + dt;
            Xk = Xk + a(t,Xk)*dt + b(t,Xk).*dW;
        end
        X(:,k+1) = Xk;
    end
end