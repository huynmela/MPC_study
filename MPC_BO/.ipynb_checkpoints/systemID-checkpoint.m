clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script identifies a state space system.

% We assume an observable, canonical form of 
% * x(k+1) = A*x(k) + B*u(k)                                                
% * y(k) = C*x(k) + D*u(k)                             
% where C = I and D = 0.

% The state-space model is defined in terms of deviation variables around 
% a nominal operating condition.

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
filename = '.csv';

y_idxs = [1,2]; % row/column indices in the data file corresponding to the output data
u_idxs = [3,4]; % row/column indices in the data file corresponding to the input data
y_labels = {'T (^\circC)', 'I (arb. units.)'}; % outputs
u_labels = {'P (W)', 'q (SLM)'}; % inputs

Ts = 0.5; % sampling time

modelOrder = 5;
est_function = 'n4sid'; 

norm_intensity = 1; % 1 for yes, 0 for no
I_norm_factor = 0.5e5; % intensity normalization factor
T_col = 1; 
I_col = 2;

plot_fit = 1; % 1 for yes, 0 for no; plot comparison of data/identified model
center_data = 1; % 1 for yes, 0 for no

num_pts2center = 60; % number of points to use to center the data

saveModel = 0; % 1 for yes, 0 for no
out_filename = ['/data/mat/subspace_id.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT:                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD, FORMAT, and CLEAN DATA %
data = readmatrix(filename); % load in data

% CLEAN UP DATA HERE. Majority of this is data/system specific.
data = data(120:end, :); % remove startup data

if norm_intensity
    if isempty(I_col)
        warning('Normalization not possible. Intensity row/col not specified.')
    else
        data(:, I_col) = data (:, I_col) ./ I_norm_factor;
    end
end

% split data into input and output data
udata = data(:, u_idxs);
ydata = data(:, y_idxs);

subIDdata = iddata(ydata, udata, Ts);
Ndata = subIDdata.N; % amount of data collected, based on frequency Ts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL VERIFICATION: Plot & identify model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define estimator system
sys = n4sid(subIDdata, modelOrder, 'DisturbanceModel', 'none', 'Form', 'canonical', 'Ts', Ts); 
% get matrices from est.
A = sys.A;
B = sys.B;
C = sys.C;

% verify model graphically
disp('Verifying model graphically.')
simTime = 0:Ts:Ts*(Ndata - 1);
% Plot simulated time response of dynamic system to arbitrary inputs
yCompare = lsim(sys, udata, simTime);

% Create defult options for comparing Name, Value pairs
opt = compareOptions('InitialCondition', zeros(modelOrder,1));

figure(1)
compare(subIDdata, sys, opt)
xlabel('Time / s')
legend('Experimental data', 'Linear model')
title('Trained Model')
set(gcf, 'color', 'w')


% validate model
wmaxTrain = max(ydata-yCompare);
wminTrain = min(ydata-yCompare);
wmaxValid = zeros(1,ny);
wminValid = zeros(1,ny);

% determine max and min errors
maxErrors = max([wmaxTrain; wmaxValid], [], 1);
minErrors = min([wminTrain; wminValid], [], 1);
disp(['Maximum Output Errors: ', num2str(maxErrors)])
disp(['Minimum Output Errors: ', num2str(minErrors)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE DATA: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataInfo.ylabels = y_labels;
dataInfo.uLabels = u_labels;
dataInfo.fileName = filename;
dataInfo.ydata = ydata;
dataInfo.udata = udata;
dataInfo.sys = sys;
dataInfo.ypred = yCompare;

% specific to K. Chan system
dataInfo.samplingTime = Ts;

if saveModel
    if isempty(out_filename)
        out_filename = ['empty_', filedate, '.mat'];
    end

    if isfile(out_filename)
        overwrite_file = input('Warning: File already exists in current path! Do you want to overwrite? [1 for yes, 0 for no]: ');
        if overwrite_file
            save(out_filename, 'A', 'B', 'C', 'yss', 'uss', 'maxErrors', 'minErrors', 'dataInfo')
        else
            out_filename = input('Input a new filename (ensure .mat is included in your filename) or press Enter if you no longer want to save the identified model: \n', 's');
            if isempty(out_filename)
                disp('Identified system not saved.')
            else
                save(out_filename, 'A', 'B', 'C', 'yss', 'uss', 'maxErrors', 'minErrors', 'dataInfo')
            end
        end
    else
        save(out_filename, 'A', 'B', 'C', 'yss', 'uss', 'maxErrors', 'minErrors', 'dataInfo')
    end
end
