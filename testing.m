%% MAIN SCRIPTS AGENT TESTING
clc
clearvars
close all
% Initial settings
addpath(genpath("data/"))
addpath(genpath("src"))
run plot_settings.m
do_training = false;
% Load configuration
config = jsondecode(fileread("src/config/config.json")); % Load JSON
if ~ismember(config.agent_type, ["SAC","DQN"])
    error('Agent not valid. Please select between "SAC", "DQN"')
elseif ~ismember(config.env_type, ["Matlab","Simulink"])
    error('Environment not valid. Please select between "Matlab", "Simulink"')
end
agent_type = config.agent_type;
env_type = config.env_type;
results_folder = "results/"+agent_type+env_type+"/";
figure_folder = "figures/"+agent_type+env_type+"/";
config_add.results_folder = results_folder;
% Agent loading
ds = fileDatastore(fullfile(results_folder,"saved_agents/"),"ReadFcn",@load,"FileExtensions",".mat");
agent_names = cell(length(ds.Files),1);
for i = 1:length(ds.Files)
    name = strsplit(ds.Files{i},"\");
    agent_names{i} = name{end};
end
[sel_agent,~]=listdlg('PromptString','SELECT THE AGENT TO BE TESTED','ListString',agent_names,"ListSize",[400 600],"Name",'Test agent selection');
load(ds.Files{sel_agent})
% Driving cycle selection
driving_cycles = load("DrivingCycleObjects.mat");
drv_cyc = fieldnames(driving_cycles);
[sel_cyc,~]=listdlg('PromptString','SELECT THE TEST DRIVING CYCLE','ListString',drv_cyc,"ListSize",[400 600],"Name",'Test cycle selection');
rm_drv_cyc = fieldnames(driving_cycles);
rm_drv_cyc(rm_drv_cyc == convertCharsToStrings(drv_cyc{sel_cyc})) = [];
driving_cycles = rmfield(driving_cycles,rm_drv_cyc);
VehicleEnv.driving_cycles = driving_cycles.(drv_cyc{sel_cyc});
if VehicleEnv.driving_cycles.Ts ~= VehicleEnv.Ts
    VehicleEnv.driving_cycles.resample(VehicleEnv.Ts);
end
if strcmp(env_type,'Matlab')
    env.driving_cycles = VehicleEnv.driving_cycles;
end
% Simulation
rng(123,"twister")
cycleTime = VehicleEnv.driving_cycles.time_simu(end);
stepSize = VehicleEnv.Ts;
maxSteps = cycleTime/stepSize;
saved_agent.UseExplorationPolicy = false;
simOpts = rlSimulationOptions;
simOpts.MaxSteps = maxSteps;
simOpts.UseParallel = false;
simOpts.SimulationStorageType = "memory";
if strcmp(env_type,'Matlab')
    sim(env,saved_agent,simOpts);
    VehicleEnv.outputPostprocess();
else
    agentSampleTime = saved_agent.AgentOptions.SampleTime;
    sim_out = sim(env,saved_agent,simOpts);
    VehicleEnv.outputPostprocess(sim_out);
end
% Plot the results
figure_folder = figure_folder+agent_names(sel_agent)+'/'+convertStringsToChars(driving_cycles.(drv_cyc{sel_cyc}).name)+'/';
figure_folder = erase(figure_folder,".mat");
if ~isfolder(figure_folder)
    mkdir(figure_folder)
end
VehicleEnv.plotTestOutput(figure_folder);
