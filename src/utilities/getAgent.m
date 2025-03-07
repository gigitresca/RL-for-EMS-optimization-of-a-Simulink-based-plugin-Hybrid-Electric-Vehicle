function [agent, agent_name, fileLogger] = getAgent(agent_type,config_add)
    % Function to extract the selected agent, given the type
    switch agent_type
    case 'SAC'
        if exist("src/config/SAC.json", 'file') ~= 2
            error('Config file does not exist: src/config/SAC.json! Please add');
        end
        config = jsondecode(fileread("src/config/SAC.json")); % Load JSON
        agent = SAC( ...
            n_nodes_base_act=config.parameters.n_nodes_base_act, ...
            n_nodes_base_cr=config.parameters.n_nodes_base_cr, ...
            lr_act=config.parameters.lr_act, ...
            lr_cr1=config.parameters.lr_cr1, ...
            lr_cr2=config.parameters.lr_cr2, ...
            ent_wgt=config.parameters.ent_wgt, ...
            ent_trg=config.parameters.ent_trg, ...
            ent_lr=config.parameters.ent_lr, ...
            max_episodes=config.max_episodes, ...
            sample_time=config.sample_time, ...
            cycle_time=config.cycle_time, ...
            empty_obs_name=true ...
        );
    case 'DQN'
        if exist("src/config/DQN.json", 'file') ~= 2
            error('Config file does not exist: src/config/DQN.json! Please add');
        end
        config = jsondecode(fileread("src/config/DQN.json")); % Load JSON
        agent = DQN( ...
            n_nodes_base_cr=config.parameters.n_nodes_base_cr, ...
            lr_cr=config.parameters.lr_cr, ...
            max_episodes=config.max_episodes, ...
            sample_time=config.sample_time, ...
            cycle_time=config.cycle_time, ...
            empty_obs_name=true ...
        );
    otherwise
        error('Unsupported agent type: %s! Please add the agent in the getAgent.m file and add the .json config file', agent_type);
    end
    agent_name = set_agent_name(config,agent.trainingOpts.MaxEpisodes,config_add.results_folder);
    fileLogger = agent.setLogger(fullfile(config_add.results_folder,"rl_logs/"),agent_name);
end