function agent_name = set_agent_name(config,name_feature,results_folder)
    date = yyyymmdd(datetime('today'));
    agent_name = date+"_Agent"+name_feature;
    if config.agent_name
        agent_name = config.agent_name;
    else
        while isfile(results_folder+"saved_agents/"+agent_name+".mat")
            agent_name_split = strsplit(agent_name,"_");
            if isnan(str2double(agent_name_split(end)))
                agent_name_suff = "_1";
                agent_name = agent_name+agent_name_suff;
            else
                agent_name_suff = num2str(str2double(agent_name_split(end))+1);
                agent_name = "";
                for i = 1:length(agent_name_split)-1
                    agent_name = strcat(agent_name, agent_name_split(i));
                    agent_name = agent_name+"_";
                end
                agent_name = agent_name+agent_name_suff;
            end
        end
    end
    if isfolder(results_folder+"rl_logs/"+agent_name)
        error("The agent name was already selected! Please select a different name")
    end
end