classdef CriticDQNNet < matlab.mixin.Copyable
    properties
        n_nodes_base
        net
    end

    methods
        function this = CriticDQNNet(n_nodes_base_act,obs_dim,act_dim)
            % Constructor method: Intializes the values of the reward object
            % (the maps are set with the proper method)
            if nargin>0
                this.setCriticNet(n_nodes_base_act,obs_dim,act_dim);           
            end
        end
        function this = setCriticNet(this,n_nodes_base_cr,obs_dim,act_dim)
            % Method to initialize the actor network
            this.n_nodes_base = n_nodes_base_cr;
            criticNet = [
                featureInputLayer(obs_dim,'Normalization','none','Name','obs')
                fullyConnectedLayer(n_nodes_base_cr*2,'Name','fc1','WeightsInitializer','he')
                reluLayer('Name','relu1')
                fullyConnectedLayer(n_nodes_base_cr,'Name','fc2','WeightsInitializer','he')
                reluLayer('Name','relu2')
                fullyConnectedLayer(act_dim,'Name','fc3')
            ];
            criticNet = dlnetwork(criticNet);
            this.net = criticNet;
        end
        function netSummary(this)
            % Method to display net parameters
            summary(this.net)
            plot(layerGraph(this.net));
        end
    end
end