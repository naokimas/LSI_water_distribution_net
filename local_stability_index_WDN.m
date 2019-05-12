% Calculate the local stability index of water distribution networks
% At minimum, one should provide the following as input:
%     - adjacency matrix, A (or incidence matrix, C), of the network
%     - diameter of each pipe, diam
%     - length of each pipe, ell
%     - number of reservoirs (i.e. nodes whose total energy is fixed), N0
%     - water demand at each node, tildeQ
%
cd '/Users/masuda/Google Drive/social/water'
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

if_pipe_removal = 0; % =1 if we carry out a pipe-removal experiment

load('AllInfo80Nets.mat')
% In ‘Property_Matrix’, the number of tanks refers to the total number of tanks
% and reservoirs. Because these networks have no tanks, the number
% refers to the number of reservoirs. 
Nnet = size(Adjacency_matrix_all,1); % number of networks

LSI_all = []; % local stability index for all networks
for ind=1:Nnet

    A = Adjacency_matrix_all{ind};
    tildeQ = NodeDemandAll{ind}; % demand at each of the N nodes
    N = size(NodeDemandAll{ind},2); % number of nodes
    M = size(LinkinfoAll{ind},1); % number of edges
    node1 = LinkinfoAll{ind}(:,1); % node on one end of the edge
    node2 = LinkinfoAll{ind}(:,2); % node on the other end of the edge
    diam = LinkinfoAll{ind}(:,3) * 1e-3; % pipe diameter [m]. originally provided in [mm]
    ell = LinkinfoAll{ind}(:,4); % pipe length [m]
    NodeType = NodeTypeAll{ind}; % 0: junction/demand node, 1: reservoir, 2: tank
    N0 = size(find(NodeType==1),2); % number of nodes whose energy is fixed (i.e. reservoirs)
    if (sum(find(NodeType==1) < N-N0+0.01) > 0.5)
        % There is at least one reservoir whose node index <= N-N0.
        % The rest of this code assumes that the indices for the
        % reservoir nodes are: N-N0+1, ..., N.
        error('Indices of reservoir nodes should be N-N0+1, ..., N') 
    end
    if size(find(NodeType==2),2) >= 1 % We do not treat tanks in this paper
        error('There is at least one tank; there should be no tank')
    end
    tildeQ = tildeQ(find(NodeType==0)); % demand at each of the N-N0 nodes whose energy is not fixed. In other words, exclude the reservoir nodes.

    tildeQ = tildeQ'; % now N-N0 \times 1 matrix
    tildeQ = tildeQ / 3600; % original data: [m^3/h]. Now transformed to [m^3/s]

    % validation of the number of nodes and edges
    N_tmp = size(A,1);
    M_tmp = nnz(A)/2; % number of edges = number of nonzero elements / 2
    if N ~= N_tmp
        [N N_tmp]
        error ('Number of nodes is inconsistent')
    end
    if M_tmp ~= M
        [M M_tmp]
        error('Number of edges is inconsistent');
    end

    C = zeros(N,M); % incidence matrix
    for i=1:M
        C(node1(i),i) = 1;
        C(node2(i),i) = -1;
    end
    if (if_connected(C) == false)
        error('The network is not a connected network');
    end

    % initialization
    h1 = 100 * ones(N-N0,1); % an arbitrary initial condition for the Newton-Raphson iteration
    h2 = 65 * ones(N0,1); % height at reservoirs. 65 [m] is the default value we used in our paper

    [h1, Q] = Newton_Raphson_water(C,tildeQ,h2,ell,diam); % steady state
    % sectional_area = pi * diam .^ 2 / 4; % cross-sectional area [m^2]
    % Re = abs(Q) .* diam ./ sectional_area ./ v; % Reynolds number
    % pic = ell ./ (g * sectional_area); % pipe inertia constant
    LSI = local_stab_index(C,N0,Q,ell,diam); % local stability index

    % local stability index after removing a pipe (not used in the paper)
    if if_pipe_removal == 1
        results_pipe_removal = [];
        for i=1:M
            C_tmp = C;
            C_tmp(:,i) = []; % remove the i-th edge
            if (if_connected(C_tmp) == true)
                ell_tmp = ell;
                ell_tmp(i) = [];
                diam_tmp = diam;
                diam_tmp(i) = [];
                [h1, Q_tmp] = Newton_Raphson_water(C_tmp,tildeQ,h2,ell_tmp,diam_tmp);
                LSI = local_stab_index(C_tmp,N0,Q_tmp,ell_tmp,diam_tmp);
                results_pipe_removal = [results_pipe_removal; LSI i ell(i) Q(i)/sectional_area(i) Re(i)/10^4 pic(i)];
            end
        end
        results_pipe_removal(:,[1 4:end])
    end % pipe removal experiment done

    fprintf('%d %d %d %d %f\n', ind, N, N0, M, LSI)
    LSI_all = [LSI_all; LSI];
end % all networks done

results_Meng2018 = xlsread('Results_80Nets');
LSI_vs_Meng2018 = []; % corr(local stability index, strain index) (i=1:6)
% or corr(local stability index, network property) (i=7:14)
for i=1:14 % 14 indices measured in Meng et al., Water Research (2018)
    [rho,pval] = corr(LSI_all, results_Meng2018(:,i));
    LSI_vs_Meng2018 = [LSI_vs_Meng2018; rho pval];
end
LSI_vs_Meng2018
    
[rho, pval] = corr(results_Meng2018(:,3), results_Meng2018(:,5));
fprintf('Corr(failure magnitude, recovery rate) = %f, (p = %f)\n', rho, pval)
    
[rho, pval] = corr(results_Meng2018(:,4), results_Meng2018(:,5));
fprintf('Corr(failure rate, recovery rate) = %f, (p = %f)\n', rho, pval)

NetProp_vs_RecoveryRate = [];
for i=1:8
    % results_Meng2018(:,5): recovery rate of each network
    % results_Meng2018(:,6+i): 8 structural properties for each network (i=1:8)
    [rho,pval] = corr(results_Meng2018(:,5), results_Meng2018(:,6+i));
    NetProp_vs_RecoveryRate = [NetProp_vs_RecoveryRate; rho pval];
end
NetProp_vs_RecoveryRate

% save('local_stab_wds_all.txt', 'LSI_all', '-ascii')
clear
    

function [h1, Q] = Newton_Raphson_water(C,tildeQ,h2,ell,diam)
%
% Calculate the steady state
%
% Input
%   C: incidence matrix
%   tildeQ: demand at each node whose energy is not fixed (N-N0\times 1 matrix)
%   h2: total energy at the nodes whose energy is fixed (N0\times 1 matrix)
%   ell: length of each pipe (M\times 1 matrix)
%   diam: diameter of each pipe (M\times 1 matrix)
%
% Output
%   h1: energy at N-N0 nodes in the steady state (N-N0\times 1 matrix)
%   Q: flor rate for the M pipes in the steady state (M\times 1 matrix)
%
g = 9.80665; % gravitational acceleration [m/s^2]

N = size(C,1); % number of nodes
M = size(C,2); % number of edges
N0 = size(h2,1); % number of nodes whose energy is fixed

% initialization for the Newton-Raphson iteration
h1 = 100 * ones(N-N0,1);
Q = 1.0 * ones(M,1);
error = 1e+2;

% tolerance for the Newton-Raphson iteration
tolerance = 0.01;

while error > tolerance
    Rm = 8 * Bellos(Q,diam) .* ell ./ (g * pi^2 * (diam .^ 5)); % pipe coefficient
    Fnode = C(1:N-N0, :) * Q + tildeQ;
    h = [h1; h2]; % total energy at each node
    Fedge = C' * h - Rm .* Q .* abs(Q);
    J = [zeros(N-N0, N-N0) C(1:N-N0, :) ; C(1:N-N0, :)' -2 * diag(Rm .* abs(Q))];
    v = [h1; Q] - inv(J) * [Fnode; Fedge]; % Newton-Raphson
    h1_new = v(1:N-N0);
    Q_new = v(N-N0+1:end);
    
    error = sum((h1_new - h1) .^ 2) + sum((Q_new - Q) .^ 2);
    h1 = h1_new;
    Q = Q_new;
end
end


function LSI = local_stab_index(C,N0,Q,ell,diam)
%
% Local stability index
%
% Input
%   C: incidence matrix (N\times M matrix)
%   N0: number of nodes whose energy is fixed (i.e. reservoirs)
%   Q: steady-state flow rate in each pipe (M\times 1 matrix)
%   ell: length of each pipe (M\times 1 matrix)
%   diam: diameter of each pipe (M\times 1 matrix)
%
% Output
%   LSI: local stability index
%
N = size(C,1); % number of nodes
M = size(C,2); % number of edges (pipes)
g = 9.80665; % gravitational acceleration [m/s^2]
sectional_area = pi * diam .^ 2 / 4; % cross-sectional area [m^2]
pic = ell ./ (g * sectional_area); % pipe inertia constant
D = diag(1 ./ pic); % diagonal matrix
Rm = 8 * Bellos(Q,diam) .* ell ./ (g * pi^2 * (diam .^ 5)); % pipe coefficient
dRm_dQ = 2 * Rm .* abs(Q) + ...
    8 * Bellos_derivative(Q,diam) .* ell ./ (g * pi^2 * diam.^5) .* Q .* abs(Q);
Jdyn = (D * C(1:N-N0,:)' * inv(C(1:N-N0,:) * D * C(1:N-N0,:)') * C(1:N-N0,:) * D - D) ...
    * 2 * diag(dRm_dQ);
eig_real = real(eig(Jdyn));
eig_real = sort(eig_real); % real part of the eigenvalues in ascending order

num_zero_eigs = sum(eig_real < 1e-8 & eig_real > -1e-8); % number of zero eigenvalues with the threshold of 1e-8
if num_zero_eigs ~= N-N0
    error('J^{dyn} should have N-N0 zero eigenvalues')
end

eig_real = eig_real(1:M-N+N0); % zero eigenvalues removed
LSI = - max(eig_real); % spectral gap (as a positive value)
end


function f = Bellos(Q,diam)
%
% Friction factor for a pipe.
% The formula proposed in Bellos et al. 2018 is used.
%
% Input
%   Q: flow rate through the pipe
%   diam: diameter of the pipe
%
% Output
%   f: friction factor
%
e = 2.591 * 1e-4; % roughness coefficient of cast iron [m]
v = 1.007 * 1e-6; % kinetic viscosity of water at 20 degree celcius [m^2/s]
sectional_area = pi * diam .^ 2 / 4; % cross-sectional area [m^2]
Re = abs(Q) .* diam ./ sectional_area ./ v; % Reynolds number
a = (1+(Re/2712).^(8.4)).^(-1);
b = (1 + (Re*e/150 ./ diam).^(1.8)).^(-1);

% Friction factor, f = f1 * f2 * f3
f1 = (64 ./ Re).^a;
f2 = (0.75*log(Re/5.37)).^(2*(a-1) .* b);
f3 = (0.88*log(6.82 * diam ./ e)).^(2*(a-1) .* (1-b));
f = f1 .* f2 .* f3;
end


function df_dQ = Bellos_derivative(Q,diam)
%
% Derivative of the friction factor, f, with respect to the flow rate, Q
% A single pipe is assumed.
% For f(Q), the formula proposed in Bellos et al. 2018 is used.
%
% Input
%   Q: flow rate through the pipe
%   diam: diameter of the pipe
%
% Output
%   df_dQ: df/dQ
%
e = 2.591 * 1e-4; % roughness coefficient of cast iron [m]
v = 1.007 * 1e-6; % kinetic viscosity of water at 20 degree celcius [m^2/s]
sectional_area = pi * diam .^ 2 / 4; % cross-sectional area [m^2]
Re = abs(Q) .* diam ./ sectional_area ./ v; % Reynolds number

a = (1+(Re/2712).^(8.4)).^(-1);
b = (1 + (Re*e/150 ./ diam).^(1.8)).^(-1);
dRe_dQ = sign(Q) .* diam ./ sectional_area ./ v; % dRe/dQ. Not continuous at Q=0 though
da_dRe = -8.4 * (Re/2712).^(7.4) / 2712 .* (a .^ 2); % da/dRe
db_dRe = -1.8 * (Re*e/150 ./ diam).^(0.8) .* (e/150 ./ diam) .* (b .^ 2); % db/dRe

% Friction factor, f = f1 * f2 * f3
f1 = (64 ./ Re).^a;
f2 = (0.75*log(Re/5.37)).^(2*(a-1) .* b);
f3 = (0.88*log(6.82 * diam ./ e)).^(2*(a-1) .* (1-b));

df1_dRe = f1 .* (da_dRe .* log(64 ./ Re) - a ./ Re); % df1/dRe
df2_dRe = f2 .* (2*(da_dRe .* b + (a-1) .* db_dRe) .* log(0.75*log(Re/5.37)) ...
    + 2*(a-1) .* b ./ log(Re/5.37) ./ Re); % df2/dRe
df3_dRe = f3 .* 2 .* (da_dRe .* (1-b) - (a-1) .* db_dRe) .* log(0.88*log(6.82*diam/e)); % df3/dRe
df_dQ = (df1_dRe .* f2 .* f3 + f1 .* df2_dRe .* f3 + f1 .* f2 .* df3_dRe) .* dRe_dQ; % df/dQ
end


function out = if_connected(C)
%
% Check if the network is a connected network.
%
% Input
%   C: incidence matrix
%
% Output
%   out = true if connected. out = falst otherwise
%
L = C * C'; % Laplacian matrix
eig_ascending = sort(eig(L)); % eigenvalue of L in ascending order
if eig_ascending(2) > 1e-8 % 2nd smallest eigenvalue is positive
    out = true;
else
    out = false;
end
end