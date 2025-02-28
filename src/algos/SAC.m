classdef SAC < matlab.mixin.Copyable
    properties
        obs_info
        act_info
        actor = ActorGaussianNet
        critic1 = CriticNet
        critic2 = CriticNet
        agent
        trainingOpts = rlTrainingOptions
    end

    methods
        function this = SAC(param)
            % Constructor method: Set obs_info and act_info of the agent
            arguments (Input)
                param.n_nodes_base_act double = 16
                param.n_nodes_base_cr double = 16
                param.lr_act double = 0.01
                param.lr_cr1 double = 0.01
                param.lr_cr2 double = 0.01
                param.ent_wgt double = 2
                param.ent_trg double = -1
                param.ent_lr double = 0.01
                param.max_episodes double = 500
                param.sample_time double = 1
                param.cycle_time double = 1800
                param.empty_obs_name {mustBeNumericOrLogical} = false
            end
            % Initialize Observation settings
            this.obs_info=rlNumericSpec([5 1]);
            this.obs_info.Name = [
                "veh_spd"
                "veh_acc"
                "em_spd"
                "batt_soc"
                "veh_dist_perc"
            ];
            % Initialize Action settings   
            this.act_info=rlNumericSpec([1 1]);
            this.act_info.Name="ice_trq_max_perc";
            this.act_info.UpperLimit = 1;
            this.act_info.LowerLimit = 0;
            this.setSAC( ...
                param.n_nodes_base_act, ...
                param.n_nodes_base_cr, ...
                param.lr_act, ...
                param.lr_cr1, ...
                param.lr_cr2, ...
                param.ent_wgt, ...
                param.ent_trg, ...
                param.ent_lr, ...
                param.max_episodes, ...
                param.sample_time, ...
                param.cycle_time/param.sample_time, ...
                param.empty_obs_name ...
            );                
        end

        function this = setSAC(this,n_nodes_base_act,n_nodes_base_cr,lr_act,lr_cr1,lr_cr2,ent_wgt,ent_trg,ent_lr,max_episodes,sample_time,max_steps,empty_obs_name)
            % Set the actor and critic nets and the SAC agent
            obs_dim = this.obs_info.Dimension(1);
            act_dim = this.act_info.Dimension(1);
            % Set the DNN
            actorNet = ActorGaussianNet(n_nodes_base_act,obs_dim,act_dim);
            criticNet = CriticNet(n_nodes_base_cr,obs_dim,act_dim);
            criticNet1 = initialize(criticNet.net);
            criticNet2 = initialize(criticNet.net);
            % Set the actor and the critic
            if empty_obs_name
                obs_info_name = this.obs_info.Name;
                this.obs_info.Name = "";
                for i = 1:length(obs_info_name)
                    this.obs_info.Name = this.obs_info.Name+obs_info_name(i)+"/";
                end
            end
            this.actor=rlContinuousGaussianActor(actorNet.net,this.obs_info,this.act_info,"ObservationInputNames","obs","ActionMeanOutputNames","mean_out","ActionStandardDeviationOutputNames","std_out");
            this.critic1 = rlQValueFunction(criticNet1,this.obs_info,this.act_info,"ActionInputNames","act","ObservationInputNames","obs");
            this.critic2 = rlQValueFunction(criticNet2,this.obs_info,this.act_info,"ActionInputNames","act","ObservationInputNames","obs");
            % Set the agent
            this.agent=rlSACAgent(this.actor,[this.critic1,this.critic2]);
            if empty_obs_name
                this.obs_info.Name = obs_info_name;
            end
            % Set agent options
            this.setActorOpts(lr_act);
            this.setCriticsOpts(lr_cr1,lr_cr2);
            this.setEntropyOpts(ent_wgt,ent_trg,ent_lr);
            this.setAgentOpts(sample_time,max_steps);
            % Set training options
            this.setTrainingOpts(max_episodes,max_steps);
        end

        function this = setActorOpts(this,lr)
            % Set actor options
            this.agent.AgentOptions.ActorOptimizerOptions.LearnRate = lr;
            this.agent.AgentOptions.ActorOptimizerOptions.L2RegularizationFactor = 1e-4;
            this.agent.AgentOptions.ActorOptimizerOptions.GradientThreshold =1;
            this.agent.AgentOptions.ActorOptimizerOptions.Algorithm="adam";
            this.agent.AgentOptions.InfoToSave.Optimizer = false;
            this.agent.AgentOptions.InfoToSave.PolicyState = false;
            this.agent.AgentOptions.InfoToSave.PolicyState = false;
        end

        function this = setCriticsOpts(this,lr1,lr2)
            % Set critics options
            this.agent.AgentOptions.CriticOptimizerOptions(1,1).LearnRate = lr1;
            this.agent.AgentOptions.CriticOptimizerOptions(1,1).L2RegularizationFactor = 1e-4;
            this.agent.AgentOptions.CriticOptimizerOptions(1,1).GradientThreshold =1;
            this.agent.AgentOptions.CriticOptimizerOptions(1,1).Algorithm="adam";
            this.agent.AgentOptions.CriticOptimizerOptions(1,2).LearnRate = lr2;
            this.agent.AgentOptions.CriticOptimizerOptions(1,2).L2RegularizationFactor = 1e-4;
            this.agent.AgentOptions.CriticOptimizerOptions(1,2).GradientThreshold =1;
            this.agent.AgentOptions.CriticOptimizerOptions(1,2).Algorithm="adam";
        end

        function this = setEntropyOpts(this,weigth,trg,lr)
            % Set entropy
            this.agent.AgentOptions.EntropyWeightOptions.EntropyWeight = weigth;
            this.agent.AgentOptions.EntropyWeightOptions.TargetEntropy = trg;
            this.agent.AgentOptions.EntropyWeightOptions.LearnRate = lr;
            this.agent.AgentOptions.EntropyWeightOptions.GradientThreshold = Inf;
            this.agent.AgentOptions.EntropyWeightOptions.Algorithm = "adam";
        end

        function this = setAgentOpts(this,sample_time,max_steps)
            % Set main agent options
            this.agent.AgentOptions.SequenceLength = 1;
            this.agent.AgentOptions.SampleTime = sample_time;
            this.agent.AgentOptions.DiscountFactor = 0.99;
            this.agent.AgentOptions.MiniBatchSize = 64;  %%%%%%%%%%%
            % this.agent.ExperienceBuffer = rlPrioritizedReplayMemory(obs_info,act_info);
            this.agent.AgentOptions.ExperienceBufferLength = 100*max_steps;
            this.agent.AgentOptions.PolicyUpdateFrequency = 1;
            this.agent.AgentOptions.CriticUpdateFrequency = 1;
            this.agent.AgentOptions.TargetSmoothFactor = 1;
            this.agent.AgentOptions.TargetUpdateFrequency = 500;
        end

        function this = setTrainingOpts(this,max_episodes,max_steps)
            % Set training options
            this.trainingOpts.MaxEpisodes=max_episodes;
            this.trainingOpts.MaxStepsPerEpisode=max_steps;
            this.trainingOpts.ScoreAveragingWindowLength=5;
            this.trainingOpts.Verbose=true;
            this.trainingOpts.Plots="none";
            this.trainingOpts.StopTrainingCriteria="EpisodeCount";
            this.trainingOpts.StopTrainingValue=max_episodes;
            this.trainingOpts.StopOnError="on";
            this.trainingOpts.SaveAgentCriteria="None";
            % this.trainingOpts.SaveAgentValue=maxEpisodes;
        end

        function fileLogger = setLogger(this,path,agent_name)
            % Function to set the filelogger to be used during the training
            fileLogger = rlDataLogger();
            fileLogger.LoggingOptions.LoggingDirectory = fullfile(path,agent_name);
            fileLogger.LoggingOptions.FileNameRule = "episode<id>";
            fileLogger.LoggingOptions.DataWriteFrequency = 10;
            fileLogger.EpisodeFinishedFcn = @this.EpisodeLoggingFcn;
            fileLogger.AgentStepFinishedFcn = @this.AgentStepLoggingFcn;
            fileLogger.AgentLearnFinishedFcn = @this.AgentLearnLoggingFcn;
        end

        function plotTrainingData(this,saved_agents_results, results_folder, figure_folder_name, agent_name, Mercedes)
            % Funtion to plot the reward trend and the losses
            run plot_settings.m
            avg_window = saved_agents_results.TrainingOptions.ScoreAveragingWindowLength;
            reward_mean = movmean(saved_agents_results.EpisodeReward,avg_window);
            reward_std = movstd(saved_agents_results.EpisodeReward,avg_window);
            % Load actor and critic log
            ds = fileDatastore(fullfile(results_folder,"rl_logs/",agent_name),"ReadFcn",@load,"FileExtensions",".mat");
            episode = 1:25:(length(ds.Files)-1);
            % Color map definition
            colormap default
            mycolormap = colormap;
            plot_multiplier = floor(length(mycolormap)/length(episode));
            % Reward plot
            h = figure;
            figure_name = 'Reward-training';
            plot(saved_agents_results.EpisodeIndex,saved_agents_results.EpisodeReward);
            hold on
            plot(saved_agents_results.EpisodeIndex,saved_agents_results.AverageReward);
            plot(saved_agents_results.EpisodeIndex,saved_agents_results.EpisodeQ0);
            hold off
            grid on
            xlabel("Episode")
            ylabel("Reward [-]")
            title("Cumulative reward")
            xlim([0 length(saved_agents_results.EpisodeIndex)])
            legend('Reward','Average reward','Q0',"Location","northoutside","Orientation","horizontal")
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Shaded reward plot
            h = figure;
            figure_name = 'Reward-training_shaded';
            upper = reward_mean + reward_std;
            bottom = reward_mean - reward_std;
            x_axis = [saved_agents_results.EpisodeIndex', fliplr(saved_agents_results.EpisodeIndex')];
            inBetween = [bottom', fliplr(upper')];
            hold on;
            fill(x_axis, inBetween,[0.6350 0.0780 0.1840],'FaceAlpha',0.3);
            plot(saved_agents_results.EpisodeIndex,reward_mean,'Color',[0.6350 0.0780 0.1840]);
            hold off
            grid on;
            xlabel("Episode")
            ylabel("Reward [-]")
            title("Cumulative reward")
            xlim([0 length(saved_agents_results.EpisodeIndex)])
            legend("Std","Mean","Location","northoutside","Orientation","horizontal")
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % SOC plot
            h = figure;
            figure_name = 'SOC-training';
            for i = 1:length(episode)
                time = seconds(saved_agents_results.SimulationOutput{episode(i)}.Time);
                soc = saved_agents_results.SimulationOutput{episode(i)}.batt_soc;
                if i == 1
                    plot(time,soc,'color',mycolormap(i,:));
                else
                    plot(time,soc,'color',mycolormap(plot_multiplier*(i-1),:));
                end
                hold on
            end
            grid on
            colorbar;
            clim([0 saved_agents_results.TrainingOptions.MaxEpisodes]); % Set the colorbar limits
            t = colorbar;
            title(t, 'Ep. \#','interpreter','latex');
            t.TickLabelInterpreter = 'latex';
            t.Label.Interpreter = 'latex';
            xlabel("Time [s]")
            ylabel("SOC [-]")
            title("State of Charge")
            xlim([0 seconds(saved_agents_results.SimulationOutput{end}.Time(end))])
            % legend("Episode = "+strings(length(episode),1)+num2str(episode'),Location="northwest")
            colororder(jet(length(episode)))
            hold off
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Fuel consumption plot
            h = figure;
            figure_name = 'Fuel consumption-training';
            for i = 1:length(episode)
                time = seconds(saved_agents_results.SimulationOutput{episode(i)}.Time);
                fc = saved_agents_results.SimulationOutput{episode(i)}.ice_fuel_cum;
                if i == 1
                    plot(time,fc,'color',mycolormap(i,:));
                else
                    plot(time,fc,'color',mycolormap(plot_multiplier*(i-1),:));
                end
                hold on
            end
            grid on
            colorbar;
            clim([0 saved_agents_results.TrainingOptions.MaxEpisodes]); % Set the colorbar limits
            t = colorbar;
            title(t, 'Ep. \#','interpreter','latex');
            t.TickLabelInterpreter = 'latex';
            t.Label.Interpreter = 'latex';
            xlabel("Time [s]")
            ylabel("$m_f$ [g]")
            title("Fuel consumption")
            xlim([0 seconds(saved_agents_results.SimulationOutput{end}.Time(end))])
            % legend("Episode = "+strings(length(episode),1)+num2str(episode'),Location="northwest")
            colororder(jet(length(episode)))
            hold off
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Critic loss plot
            h = figure;
            figure_name = 'Critic loss';
            step_init = 1;
            for i = 1:length(episode)
                load(ds.Files{i});
                if ~isempty(agentLearnData)
                    step_end = step_init+length(cell2mat(agentLearnData.CriticLoss))-1;
                    if i == 1
                        plot(step_init:1:(step_end),cell2mat(agentLearnData.CriticLoss),'color',mycolormap(i,:));
                    else
                        plot(step_init:1:(step_end),cell2mat(agentLearnData.CriticLoss),'color',mycolormap(plot_multiplier*(i-1),:));
                    end
                    hold on
                    step_init = step_end+1;
                end
            end
            grid on
            colorbar;
            clim([0 saved_agents_results.TrainingOptions.MaxEpisodes]); % Set the colorbar limits
            t = colorbar;
            title(t, 'Ep. \#','interpreter','latex');
            t.TickLabelInterpreter = 'latex';
            t.Label.Interpreter = 'latex';
            xlabel("Agent step [-]")
            xlim([0 step_end])
            ylabel("Loss [-]")
            title("Critic loss")
            %legend("Episode = "+strings(length(episode),1)+num2str(episode'),Location="northwest")
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Actor loss plot
            h = figure;
            figure_name = 'Actor loss';
            step_init = 1;
            for i = 1:length(episode)
                load(ds.Files{i});
                if ~isempty(agentLearnData)
                    step_end = step_init+length(cell2mat(agentLearnData.ActorLoss))-1;
                    if i == 1
                        plot(step_init:1:(step_end),cell2mat(agentLearnData.ActorLoss),'color',mycolormap(i,:));
                    else
                        plot(step_init:1:(step_end),cell2mat(agentLearnData.ActorLoss),'color',mycolormap(plot_multiplier*(i-1),:));
                    end
                    hold on
                    step_init = step_end+1;
                end
            end
            grid on
            colorbar;
            clim([0 saved_agents_results.TrainingOptions.MaxEpisodes]); % Set the colorbar limits
            t = colorbar;
            title(t, 'Ep. \#','interpreter','latex');
            t.TickLabelInterpreter = 'latex';
            t.Label.Interpreter = 'latex';
            xlabel("Agent step [-]")
            xlim([0 step_end])
            ylabel("Loss [-]")
            title("Actor loss")
            %legend("Episode = "+strings(length(episode),1)+num2str(episode'),Location="northwest")
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Entropy weight plot
            h = figure;
            figure_name = 'Entropy weight';
            step_init = 1;
            for i = 1:length(episode)
                load(ds.Files{i});
                if ~isempty(agentLearnData)
                    step_end = step_init+length(cell2mat(agentLearnData.AgentLearnCount))-1;
                    doubleEntropyWeights = cellfun(@double, agentLearnData.EntropyWeight, 'UniformOutput', false);
                    matrixEntropyWeights = cell2mat(doubleEntropyWeights);
                    if i == 1
                        plot(step_init:1:step_end,matrixEntropyWeights,'color',mycolormap(i,:));
                    else
                        plot(step_init:1:step_end,matrixEntropyWeights,'color',mycolormap(plot_multiplier*(i-1),:));
                    end
                    hold on
                    step_init = step_end+1;
                end
            end
            grid on
            colorbar;
            clim([0 saved_agents_results.TrainingOptions.MaxEpisodes]); % Set the colorbar limits
            t = colorbar;
            title(t, 'Ep. \#','interpreter','latex');
            t.TickLabelInterpreter = 'latex';
            t.Label.Interpreter = 'latex';
            xlabel("Agent step [-]")
            xlim([0 step_end])
            ylabel("Weight factor [-]")
            title("Entropy")
            %legend("Episode = "+strings(length(episode),1)+num2str(episode'),Location="northwest")
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Policy and critic analysis settings
            n_el = 30;
            n_act = 30;
            n_level = 3;
            critic = this.agent.getCritic;
            actor = this.agent.getActor;
            x_axis = "veh_spd";
            sweep_var = ["veh_acc" "batt_soc" "veh_dist_perc"];
            [x_axis_swp,action_swp,batch_obs_tot,batch_act,states_comb,states_comb_norm] = Mercedes.statesBatchCreation(x_axis,sweep_var,n_el,n_act,n_level);
            [XAxisSwp,ActionSwp] = meshgrid(x_axis_swp,action_swp);
            [QValue,QValueMax,QvalueMaxIdx] = computeQvalue(critic,states_comb_norm,batch_obs_tot,batch_act,n_act,n_el,this.obs_info);
            [~,action_mean,action_std] = computeActionDistribution(actor,batch_obs_tot,states_comb_norm,this.obs_info,n_el);
            % Q-value function plots
            sweep_var = strrep(sweep_var,"veh_","");
            sweep_var = strrep(sweep_var,"batt_","");
            sweep_var = strrep(sweep_var,"_perc","");
            % Fix vehicle distance
            idx_init = 1;
            idx_end = n_level*n_level;
            for j = 1:n_level
                h = figure;
                figure_name = "QValue_speed_SoC_sweep_dist"+round(states_comb{idx_init,3},2)*100;
                figure_name = char(figure_name);
                tld = tiledlayout(n_level,n_level,"TileSpacing","compact");
                for i = idx_init:idx_end
                    nexttile
                    contourf(XAxisSwp,ActionSwp,QValue(:,:,i))
                    hold on
                    plot(x_axis_swp,action_swp(QvalueMaxIdx(i,:)),"Color",[0 0 0])
                    hold off
                    colorbar
                    title(sweep_var(1)+"="+round(states_comb{i,1},2)+", "+sweep_var(2)+"="+round(states_comb{i,2},2)+", "+sweep_var(3)+"="+round(states_comb{i,3},2),"Interpreter","none")
                    set(gca,'Fontsize',print_fontsize*2/3);
                end
                xlabel(tld,"Vehicle speed [km/h]","Fontsize",print_fontsize,"Fontname",print_font,"Interpreter","latex")
                ylabel(tld,"Action [-]","Fontsize",print_fontsize,"Fontname",print_font,"Interpreter","latex")
                title(tld,"Q-value with fixed vehicle distance = "+round(states_comb{idx_init,3},2),"Fontsize",print_fontsize,"Fontname",print_font,"Interpreter","latex")
                idx_init = i+1;
                idx_end = idx_end+n_level*n_level;
                print_figure(h,print_size.*[3,3],fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
                savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            end
            % Action gaussian plot
            upper = action_mean + action_std;
            bottom = action_mean - action_std;
            x_axis_shape = [x_axis_swp, fliplr(x_axis_swp)];
            % Fix vehicle distance
            idx_init = 1;
            idx_end = n_level*n_level;
            for j = 1:n_level
                h = figure;
                figure_name = "Action_distribution_speed_SoC_sweep_dist"+round(states_comb{idx_init,3},2)*100;
                figure_name = char(figure_name);
                tld = tiledlayout(n_level,n_level,"TileSpacing","compact");
                for i = idx_init:idx_end
                    nexttile
                    inBetween = [bottom(i,:), fliplr(upper(i,:))];
                    hold on;
                    fill(x_axis_shape, inBetween,[0.6350 0.0780 0.1840],'FaceAlpha',0.3);
                    plot(x_axis_swp,action_mean(i,:),'Color',[0.6350 0.0780 0.1840]);
                    hold off
                    grid on;
                    xlim([x_axis_swp(1) x_axis_swp(end)])
                    ylim([0 1]);
                    title(sweep_var(1)+"="+round(states_comb{i,1},2)+", "+sweep_var(2)+"="+round(states_comb{i,2},2)+", "+sweep_var(3)+"="+round(states_comb{i,3},2),"Interpreter","none")
                    set(gca,'Fontsize',print_fontsize*2/3);
                end
                xlabel(tld,"Vehicle speed [km/h]","Fontsize",print_fontsize,"Fontname",print_font,"Interpreter","latex")
                ylabel(tld,"Action [-]","Fontsize",print_fontsize,"Fontname",print_font,"Interpreter","latex")
                title(tld,"Action distribution with fixed vehicle distance = "+round(states_comb{idx_init,3},2),"Fontsize",print_fontsize,"Fontname",print_font,"Interpreter","latex")
                idx_init = i+1;
                idx_end = idx_end+n_level*n_level;
                print_figure(h,print_size.*[3,3],fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
                savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            end
        end
    end

    methods (Static)
        % Logging functions for the agent
        function dataToLog = EpisodeLoggingFcn(data)
            % data is a structure that contains the following fields:
            % EpisodeCount: The current episode number.
            % Environment: Environment object.
            % Agent: Agent object.
            % Experience: A structure containing the experiences from the current episode.
            % EpisodeInfo: A structure containing the fields CumulativeReward, StepsTaken, and InitialObservation.
            % SimulationInfo: A Simulink.SimulationOutput object containing logged signals in Simulink environments.
            %
            % dataToLog is a structure containing the data to be logged to disk.           
            for i = 1:length(data.Experience)
                obsInfo = data.Environment.getObservationInfo;
                obsInfo = obsInfo.Name;
                dataToLog.SOC(i) = data.Experience(i).Observation{1}(strcmp(obsInfo,"batt_soc"));
                dataToLog.Action(i) = data.Experience(i).Action{1}(1);
            end
        end

        function dataToLog = AgentStepLoggingFcn(data)
            % data is a structure that contains the following fields:
            % EpisodeCount: The current episode number.
            % AgentStepCount: The cumulative number of steps taken by the agent.
            % SimulationTime: The current simulation time in the environment.
            % Agent: Agent object.
            %
            % dataToLog is a structure containing the data to be logged to disk.
            Epsilon = getState(getExplorationPolicy(data.Agent));
            dataToLog.noiseState = Epsilon;      
        end

        function dataToLog = AgentLearnLoggingFcn(data)
            % data is a structure that contains the following fields:
            % EpisodeCount: The current episode number.
            % AgentStepCount: The cumulative number of steps taken by the agent.
            % AgentLearnCount: The cumulative number of learning steps taken by the agent.
            % EnvModelTrainingInfo: A structure containing the fields TransitionFcnLoss, RewardFcnLos, IsDoneFcnLoss. This is applicable for model-based agent training.
            % Agent: Agent object.
            % ActorLoss: Training loss of actor function.
            % Agent: Training loss of critic function.
            %
            % dataToLog is a structure containing the data to be logged to disk.
            dataToLog.CriticLoss = min(data.CriticLoss);
            dataToLog.ActorLoss = data.ActorLoss;
            dataToLog.AgentLearnCount = data.AgentLearnCount;
            dataToLog.EntropyWeight = data.Agent.getEntropyWeight;
        end
    end
end