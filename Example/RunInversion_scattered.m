  %Run the Test6 Ensemble Kalman Inversion.

%I am running this using EnKF_Combined_Code12.
%The data being fit are horizon points, fault points, and restored-state

options_file_name = 'Options_Test6_Scattered_Forward.m';
% options_file_name = 'Options_Test6_Scattered_Restoration.m';

run(options_file_name);
[params_initial_raw,params_final_raw,info] = RunEKI(opt); %Run the inversion.
save('N200_Forward_scattered.mat','params_initial_raw','params_final_raw','info','opt');
% save('N200_Restoration_scattered.mat','params_initial_raw','params_final_raw','info','opt');