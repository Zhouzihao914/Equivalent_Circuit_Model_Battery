% A customized aging function for discharge capacity
%
% Inputs:  func_ind - function id, 1->linear, 2->exponential
%          a - parameter a for calendary aging
%          b - parameter b for dynamic aging
%          N - cycle number (linear-N=0)
%          t0 - begin time step
%          Ltime - time length
%          Q - init capacity


% Outputs: delta_Q - corresponding drop volume of discharge capacity 

function delta_Q = AgingFuncDisQ(func_ind, a, b, N, Ltime)
format long
if func_ind == 1,
    a_delta = a*[1:1:Ltime];
    b_delta = b*[1:1:Ltime];
elseif func_ind == 2,
    if b > 0,
        dyn_drop = (exp((N-60)/b) - exp((N-1-60)/b));
        b_delta = linspace(0,dyn_drop,Ltime);
    else
        b_delta = 0;
    end
    a_delta = a*[1:1:Ltime];
else
    error('ERROR: The Function type is not included!');
end
    delta_Q = a_delta + b_delta;
