clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script identifies a state space system.

% We assume an observable, canonical form of 
% * x(k+1) = A*x(k) + B*u(k)                                                
% * y(k) = C*x(k) + D*u(k)                             
% where C = I and D = 0.

% The state-space model is defined in terms of deviation variables around 
% a nominal operating condition.

% * current input follows Chan et. al. [10.23919/ACC55779.2023.10156650]   
% * where nx = 2, nu = 2. eventually, need to modify nx = 5*?

% INPUT CSV DATA REQUIREMENTS:
% * states MUST be listed first prior to inputs

% Input:                                                                  
% - file name of CSV                                                      
% - identification of number of states nx and inputs nu                   
 
% Output:                                                                 
% - A, B, C matrices corresponding to process form with state feedback    
% - nominal steady state values for centering controller                  
% - maximum and minimum errors between linear model and data              
% - performance metrics of identified model                               
% - copy of information on data used                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open loop data from experimental studies
filename = '/data/ExperimentalData/OL_data_0_inputOutputData.csv';
filedate = '2022_09_22_17h28m06s';
nx = 2; % number of states
nu = 2; % number of control inputs
data_direction = 0; % 0 for column-wise, 1 for row-wise data

modelOrder = 2;
est_function = 'n4sid'; 

plot_fit = 1; % 1 for yes, 0 for no; plot comparison of data/identified model
center_data = 1; % 1 for yes, 0 for no
num_pts2center = 60; % number of points to use to center the data

saveModel = 0; % 1 for yes, 0 for no
out_filename = ['subspace_id_', filedate, '.mat'];

% input specific to K. Chan's paper system parameters
I_bkgd_collect = 0; % background optical intensity (if collected)
sub_I_bkgd = 0; % option to subtract background optical intensity
Ts = 0.5; % Sampling time
norm_intensity = 1; % 1 for yes, 0 for no
filter_temp = 0;
I_norm_factor = 0.5e5; % intensity normalization factor
T_col = 1; 
I_col = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT:                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD, FORMAT, and CLEAN DATA %
data = readtable(filename); % load in data
if data_direction % transpose data if row-wise
    data = data'; 
end

% determine indices for data file for state and input data
numCols = width(data); % get total number of states and inputs
% obtain state and input labels
headers = data.Properties.VariableNames; % cell array of character vectors

x_idxs = zeros(1, nx); 
u_idxs = zeros(1, nu);
x_labels = {};
u_labels = {};

for i = 1:numCols
    if i < nx
        x_idxs(n) = n; % populate state indices
        x_labels{n} = headers{n}; % populate state labels
    else
        u_idxs(n) = n; % populate input indices
        u_labels{n} = headers{n}; % populate input labels
    end
end

% CLEAN UP DATA HERE. Majority of this is data/system specific.
if I_bkgd_collect
    I_bkgd = data(1:60, I_col);
    I_bkgd_sum = mean(I_bkgd);
end

data = data(120:end, :); % remove startup data

if sub_I_bkgd
    data(:, I_col) = data(:, I_col) - I_bkgd_sum;
end

if norm_intensity
    if isempty(I_col)
        warning('Normalization not possible. Intensity row/col not specified.')
    else
        data(:, I_col) = data (:, I_col) ./ I_norm_factor;
    end
end

if filter_temp
    data(2:end,T_col) = 0.3*data(1:end-1,T_col) + 0.7*data(2:end,T_col);
end

% split data into input and output data
udata = data(:, u_idxs);
xdata = data(:, x_idxs);

subIDdata = iddata(xdata, udata, Ts);
Ndata = subIDdata.N; % amount of data collected, based on frequency Ts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL VERIFICATION: Plot & identify model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plot data to visualize.')
figure(1)

for i = 1:nx % state subplot
    subplot(nx, nu, i)
    plot(data(:, x_idxs(i)), 'linewidth', 2)
    ylabel(x_labels{i})
    set(gca, 'fontsize', 15)
end
xlabel('Time step')
set(gcf, 'color', 'w');

% define estimator system
sys = n4sid(subIDdata, modelOrder, 'DisturbanceModel', 'none', 'Form', 'canonical', 'Ts', Ts); 
% get matrices from est.
A = sys.A;
B = sys.B;
C = sys.C;

% verify model graphically
if plot_fit
    disp('Verifying model graphically.')
    simTime = 0:Ts:Ts*(Ndata - 1);
    % Plot simulated time response of dynamic system to arbitrary inputs
    xCompare = lsim(sys, udata, simTime);
    
    % Create defult options for comparing Name, Value pairs
    opt = compareOptions('InitialCondition', zeros(modelOrder,1));
    
    figure(3)
    compare(subIDdata, sys, opt)
    xlabel('Time / s')
    legend('Experimental data', 'Linear model')
    title('Trained Model')
    set(gcf, 'color', 'w')
end

% validate model
wmaxTrain = max(ydata-xCompare);
wminTrain = min(ydata-xCompare);
wmaxValid = zeros(1,nx);
wminValid = zeros(1,nx);

% determine max and min errors
maxErr = max([wmaxTrain; wmaxValid], [], 1);
minErr = min([wminTrain; wminValid], [], 1);
disp(['Maximum Output Errors: ', num2str(maxErrors)])
disp(['Minimum Output Errors: ', num2str(minErrors)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE DATA: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
