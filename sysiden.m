%% System Identification - to recover h, given U (= U_0 + W_i) and y (= y_0 + w_o)
% Refer to line 31 for verification of this algorithm

clc; clear; close all
h = [-0.3 -0.9 0.8]';
snr1 = 10;
snr2 = 5;
u0 = randi([0 1], 10, 1);
u0(u0==0) = -1;
% u0 = [-1 1 -1 1 -1 1 1 -1 1 -1]';
u = awgn(u0,snr1);
U0 = toeplitz(u0,zeros(1,length(h)));
y0 = U0*h;
U = toeplitz(u,zeros(1,length(h)));
y = awgn(y0,snr2);

%% Generating Ui and yi
p = size(U,1);
n = size(U,2);
Ui = cell(1,p);
yi = eye(p);

for i = 1:p
    Ui{i} = toeplitz(yi(i,:),zeros(1,n));
end

%% Algorithm 4
rho = 0;
k = 1;
fak = ones(p,2)*diag([2 1]);

while norm(fak(:,k+1),2) < norm(fak(:,k),2)
    
    % this algorithm can be tested by replacing U with U0 and y with y0
    % expected result should be xk ~ h
    [xk fak1 alphak1] = sysidenalgo3(rho, U, y, Ui, yi);
%     [xk fak1 alphak1] = sysidenalgo3(rho, U0, y0, Ui, yi);
    
    alphak = alphak1(:,end);
    
    k = k + 1;
    fak(:,k+1) = fak1(:,end);
    
    rho = rho + 1e-3;
end

fprintf('Done!\n');
disp(xk)