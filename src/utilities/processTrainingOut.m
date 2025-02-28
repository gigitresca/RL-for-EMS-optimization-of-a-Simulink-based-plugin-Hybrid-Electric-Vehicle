function saved_agents_results = processTrainingOut(training_stats, agent_type, env_type, env_step_size, agent_name)
    % Funtion to process the output of the RL training in a format compliant to
    % the plot and testing of the codebase
    saved_agents_results.EpisodeIndex = training_stats.EpisodeIndex;
    saved_agents_results.EpisodeReward = training_stats.EpisodeReward;
    saved_agents_results.EpisodeSteps = training_stats.EpisodeSteps;
    saved_agents_results.AverageReward = training_stats.AverageReward;
    saved_agents_results.TotalAgentSteps = training_stats.TotalAgentSteps;
    saved_agents_results.AverageSteps = training_stats.AverageSteps;
    saved_agents_results.EpisodeQ0 = training_stats.EpisodeQ0;
    saved_agents_results.SimulationOutput = cell(training_stats.SimulationInfo.NumSimulations,1);
    saved_agents_results.envStepSize = env_step_size;
    saved_agents_results.TrainingOptions = training_stats.TrainingOptions;
    if env_type == "Simulink"
        for i = 1 :training_stats.SimulationInfo.NumSimulations
            saved_agents_results.SimulationOutput{i} = timetable( ...
                seconds(training_stats.SimulationInfo(i).batt_soc.Time), ...
                training_stats.SimulationInfo(i).batt_soc.Data, ...
                training_stats.SimulationInfo(i).ice_fuel_cum.Data, ...
                VariableNames=["batt_soc" "ice_fuel_cum"] ...
            );
        end
    else
        ds = fileDatastore(fullfile("results/"+agent_type+env_type+"/rl_logs/",agent_name),"ReadFcn",@load,"FileExtensions",".mat");
        while ds.hasdata
            data = read(ds);
            if contains("info",fieldnames(data))
                continue
            end
            time = 1:1:length(data.episodeData.Action{1});
            saved_agents_results.SimulationOutput{data.EpisodeCount} = timetable( ...
                seconds(time'), ...
                data.episodeData.SOC{1}', ...
                cumtrapz(time',data.episodeData.Action{1}'), ...
                VariableNames=["batt_soc" "ice_fuel_cum"] ...
            );
        end
    end
end