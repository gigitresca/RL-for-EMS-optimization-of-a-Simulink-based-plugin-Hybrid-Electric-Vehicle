classdef DQN < matlab.mixin.Copyable
    properties
        obs_info
        act_info
        critic = CriticNet
        agent
        trainingOpts = rlTrainingOptions
    end

    methods
        function this = DQN(param)
            % Constructor method: Set obs_info and act_info of the agent
            arguments (Input)
                param.n_nodes_base_cr double = 16
                param.lr_cr double = 0.01
                param.max_episodes double = 6000
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
            act_discr = 0.05;
            this.act_info=rlFiniteSetSpec(0:act_discr:1);
            this.act_info.Name="ice_trq_max_perc";
            this.setDQN( ...
                param.n_nodes_base_cr, ...
                param.lr_cr, ...
                param.max_episodes, ...
                param.sample_time, ...
                param.cycle_time/param.sample_time, ...
                param.empty_obs_name ...
            );   
        end

        function this = setDQN(this,n_nodes_base_cr,lr_cr,max_episodes,sample_time,max_steps,empty_obs_name)
            % Set the critic net and the DQN agent
            % Initialize Observation settings
            obs_dim = this.obs_info.Dimension(1);
            act_dim = numel(this.act_info.Elements);
            % Set the DNN
            criticNet = CriticDQNNet(n_nodes_base_cr,obs_dim,act_dim);
            if empty_obs_name
                obs_info_name = this.obs_info.Name;
                this.obs_info.Name = "";
                for i = 1:length(obs_info_name)
                    this.obs_info.Name = this.obs_info.Name+obs_info_name(i)+"/";
                end
            end
            this.critic=rlVectorQValueFunction(criticNet.net,this.obs_info,this.act_info);
            this.agent=rlDQNAgent(this.critic);
            if empty_obs_name
                this.obs_info.Name = obs_info_name;
            end
            % Set agent options
            this.setCriticsOpts(lr_cr);
            this.setExplorationOpts(max_episodes,max_steps);
            this.setAgentOpts(sample_time,max_steps);
            % Set training options
            this.setTrainingOpts(max_episodes,max_steps);
        end

        function this = setCriticsOpts(this,lr)
            % Set critics options
            this.agent.AgentOptions.CriticOptimizerOptions.LearnRate = lr;
            this.agent.AgentOptions.CriticOptimizerOptions.L2RegularizationFactor = 1e-4;
            this.agent.AgentOptions.CriticOptimizerOptions.GradientThreshold=1;
            this.agent.AgentOptions.CriticOptimizerOptions.Algorithm="adam";
            this.agent.AgentOptions.InfoToSave.ExperienceBuffer = false;
            this.agent.AgentOptions.InfoToSave.Optimizer = false;
            this.agent.AgentOptions.InfoToSave.PolicyState = false;
            this.agent.AgentOptions.InfoToSave.PolicyState = false;
        end

        function this = setExplorationOpts(this,max_episode,max_steps)
            % Set epsilon-greedy exploration strategy options
            this.agent.AgentOptions.EpsilonGreedyExploration.Epsilon = 1;
            this.agent.AgentOptions.EpsilonGreedyExploration.EpsilonMin = 0.01;
            min_epsilon_episode = 1/4*max_episode;
            this.agent.AgentOptions.EpsilonGreedyExploration.EpsilonDecay = ...
                1-(this.agent.AgentOptions.EpsilonGreedyExploration.EpsilonMin/this.agent.AgentOptions.EpsilonGreedyExploration.Epsilon)^(1/(min_epsilon_episode*max_steps-1));
        end

        function this = setAgentOpts(this,sample_time,max_steps)
            % Set main agent options
            this.agent.AgentOptions.SequenceLength = 1;
            this.agent.AgentOptions.UseDoubleDQN = true;
            this.agent.AgentOptions.SampleTime = sample_time;
            this.agent.AgentOptions.DiscountFactor = 0.99;
            this.agent.AgentOptions.MiniBatchSize = 64;  %%%%%%%%%%%
            % agent.ExperienceBuffer = rlPrioritizedReplayMemory(obs_info,act_info);
            this.agent.AgentOptions.ExperienceBufferLength = 100*max_steps;
            this.agent.AgentOptions.TargetSmoothFactor = 1;
            this.agent.AgentOptions.TargetUpdateFrequency = 1*max_steps/4;
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

        function this = resetObsName(this)
            this.obs_info.Name = "";
        end

        function this = epsilonPlot(this)
            % Plot exploration strategy
            agent_step_axes = (0:1:maxSteps*maxEpisodes)';
            Epsilon = this.agent.AgentOptions.EpsilonGreedyExploration.Epsilon * ones(size(agent_step_axes));
            for i = 1:length(Epsilon)-1
                Epsilon(i+1) = Epsilon(i)*(1-this.agent.AgentOptions.EpsilonGreedyExploration.EpsilonDecay);
            end
            Epsilon(Epsilon<this.agent.AgentOptions.EpsilonGreedyExploration.EpsilonMin) = this.agent.AgentOptions.EpsilonGreedyExploration.EpsilonMin;
            figure
            plot(agent_step_axes/maxSteps,Epsilon,"LineWidth",2);
            grid on
            xlabel("Episodes [-]")
            ylabel("\epsilon")
            title("Epsilon Greedy Trend")
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
            dataToLog.CriticLoss = data.CriticLoss;
            dataToLog.AgentLearnCount = data.AgentLearnCount;
        end

        function plotTrainingData(saved_agents_results, results_folder, figure_folder_name, agent_name, Mercedes)
            % Funtion to plot the reward trend and the losses
            run plot_settings.m
            avg_window = saved_agents_results.TrainingOptions.ScoreAveragingWindowLength;
            reward_mean = movmean(saved_agents_results.EpisodeReward,avg_window);
            reward_std = movstd(saved_agents_results.EpisodeReward,avg_window);
            % Load critic log
            ds = fileDatastore(fullfile(results_folder,"rl_logs/",agent_name),"ReadFcn",@load,"FileExtensions",".mat");
            episode = 1:500:(length(ds.Files)-1);
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
        end
    end
end