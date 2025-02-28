classdef ActorGaussianNet < matlab.mixin.Copyable
    properties
        n_nodes_base
        net
    end

    methods
        function this = ActorGaussianNet(n_nodes_base_act,obs_dim,act_dim)
            % Constructor method: Intializes the values of the reward object
            % (the maps are set with the proper method)
            if nargin>0
                this.setActorGaussianNet(n_nodes_base_act,obs_dim,act_dim);                
            end
        end
        function this = setActorGaussianNet(this,n_nodes_base_act,obs_dim,act_dim)
            % Method to initialize the actor network
            this.n_nodes_base = n_nodes_base_act;
            commonPath = [
                featureInputLayer(obs_dim(1),'Normalization','none','Name','obs')
                fullyConnectedLayer(n_nodes_base_act*4,'Name','fc_com1','WeightsInitializer','he')
                reluLayer('Name','relu_com1')
                fullyConnectedLayer(n_nodes_base_act*2,'Name','fc_com2','WeightsInitializer','he')
                reluLayer('Name','relu_com2')
            ];
            meanPath = [
                fullyConnectedLayer(n_nodes_base_act,'Name','fc_mean','WeightsInitializer','he')
                reluLayer('Name','relu_mean')
                fullyConnectedLayer(prod(act_dim),'Name','mean_out')
            ];
            stdPath = [
                fullyConnectedLayer(prod(act_dim),'Name','fc_std')
                reluLayer('Name','relu_std')
                softplusLayer(Name="std_out")
            ];
            % Add layers to layerGraph object 
            actorNet = layerGraph(commonPath);
            actorNet = addLayers(actorNet,meanPath);
            actorNet = addLayers(actorNet,stdPath);
            % Connect layers
            actorNet = connectLayers(actorNet,"relu_com2","fc_mean/in");
            actorNet = connectLayers(actorNet,"relu_com2","fc_std/in");
            % Convert to dlnetwork and display the number of weights.
            actorNet = dlnetwork(actorNet);
            this.net = actorNet;
        end
        function netSummary(this)
            % Method to display net parameters
            summary(this.net)
            plot(layerGraph(this.net));
        end
    end
end