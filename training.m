%% MAIN SCRIPTS FOR RL TRAINING
clc
clearvars
close all

% Initial settings
addpath(genpath("data/"))
addpath(genpath("src"))
run plot_settings.m
do_training = true;
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
config_add.results_folder = results_folder;
% Set the agent and the environment
[RLagent, agent_name, fileLogger] = getAgent(agent_type, config_add);
obs_info = RLagent.obs_info;
act_info = RLagent.act_info;
agentSampleTime = RLagent.agent.AgentOptions.SampleTime;
[env, VehicleEnv] = getEnv(env_type, obs_info, act_info);
env.reset;
% Perform the training
save_folder = fullfile(results_folder,'saved_agents');
if ~isfolder(save_folder)
    mkdir(save_folder)
end 
figure_folder = "figures/"+agent_type+env_type+"/"+agent_name+"/";
if ~isfolder(figure_folder)
    mkdir(figure_folder)
end
training_stats = train(RLagent.agent,env,RLagent.trainingOpts,Logger=fileLogger);
if training_stats.EpisodeIndex(end) == RLagent.trainingOpts.MaxEpisodes
    saved_agent = RLagent.agent;
    saved_agents_results = processTrainingOut(training_stats, agent_type, env_type, VehicleEnv.Ts, agent_name);
    save(fullfile(save_folder,agent_name),'saved_agent','saved_agents_results','env','VehicleEnv');
    RLagent.plotTrainingData(saved_agents_results, results_folder, figure_folder, agent_name, VehicleEnv)
end
