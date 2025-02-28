classdef CriticNet < matlab.mixin.Copyable
    properties
        n_nodes_base
        net
    end

    methods
        function this = CriticNet(n_nodes_base_act,obs_dim,act_dim)
            % Constructor method: Intializes the values of the reward object
            % (the maps are set with the proper method)
            if nargin>0
                this.setCriticNet(n_nodes_base_act,obs_dim,act_dim);           
            end
        end
        function this = setCriticNet(this,n_nodes_base_cr,obs_dim,act_dim)
            % Method to initialize the actor network
            this.n_nodes_base = n_nodes_base_cr;
            obsPath = [
                featureInputLayer(obs_dim(1),'Normalization','none','Name','obs')
                fullyConnectedLayer(n_nodes_base_cr*4,'Name','fc_obs','WeightsInitializer','he')
            ];
            actPath = [
                featureInputLayer(act_dim(1),'Normalization','none','Name','act')
                fullyConnectedLayer(n_nodes_base_cr*4,'Name','fc_act','WeightsInitializer','he')
            ];
            commonPath = [
                concatenationLayer(1,2,'Name','concat')
                reluLayer('Name','relu_com')
                fullyConnectedLayer(n_nodes_base_cr*2,'Name','fc_com1','WeightsInitializer','he')
                reluLayer('Name','relu_com1')
                fullyConnectedLayer(n_nodes_base_cr,'Name','fc_com2','WeightsInitializer','he')
                reluLayer('Name','relu_com2')
                fullyConnectedLayer(1,'Name','critic_out','WeightsInitializer','he')
            ];
            % Add layers to layergraph object
            criticNet = layerGraph;
            criticNet = addLayers(criticNet,obsPath);
            criticNet = addLayers(criticNet,actPath);
            criticNet = addLayers(criticNet,commonPath);
            % Connect layers
            criticNet = connectLayers(criticNet,"fc_obs","concat/in1");
            criticNet = connectLayers(criticNet,"fc_act","concat/in2");
            criticNet = dlnetwork(criticNet,"Initialize",false);
            this.net = criticNet;
        end
        function netSummary(this)
            % Method to display net parameters
            summary(this.net)
            plot(layerGraph(this.net));
        end
    end
end