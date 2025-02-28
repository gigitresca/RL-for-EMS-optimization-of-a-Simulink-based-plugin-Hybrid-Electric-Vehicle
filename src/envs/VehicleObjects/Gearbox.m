classdef Gearbox < matlab.mixin.Copyable
    % object that define the parameters of the four wheels of the vehicle

    properties
        inertia_in
        inertia_out
        gear_idx
        gear_ratio
    end

    properties (SetAccess = private)
        EffTrqInMap = Lookup2D
        EffTrqOutMap = Lookup2D
    end

    methods
        function this = Gearbox(tab_data)
            % Constructor method: Intializes the scalar value of the
            % gearbox (the efficiency maps are set with the proper method)
            if nargin>0
                this.inertia_in = tab_data.("Inertia in")(1);
                this.inertia_out = tab_data.("Inertia out")(1);
                this.gear_idx = [0;tab_data.("N gear")];
                this.gear_ratio = [0;tab_data.("Gear ratios")];
            end
        end

        function setGearboxParameters(this,tab_data)
            % Method to set the gearbox parameters (the efficiency maps are set with the proper method)
            this.inertia_in = tab_data.("Inertia in")(1);
            this.inertia_out = tab_data.("Inertia out")(1);
            this.gear_idx = [0;tab_data.("N gear")];
            this.gear_ratio = [0;tab_data.("Gear ratios")];
        end

        function initializeEfficiencyMap(this,tab_data,var_info)
            % Method to initialize the efficiency map of the gearbox
            var_name = var_info.Properties.VariableNames;
            for i = 1:length(this.gear_idx)
               this.EffTrqInMap(i,1) = Lookup2D; 
               this.EffTrqOutMap(i,1) = Lookup2D; 
            end
            Map = ones(size(tab_data)-[1 1]);
            var_name{1} = 'Speed in';
            var_name{2} = 'Torque in';
            this.EffTrqInMap(1).setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),Map,var_name,var_info.Variables);
            var_name{1} = 'Speed out';
            var_name{2} = 'Torque out';
            this.EffTrqOutMap(1).setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),Map,var_name,var_info.Variables);
        end

        function setEfficiencyMap(this,tab_data,var_info,gear)
            % Method to initialize the efficiency map of the gearbox
            var_name = var_info.Properties.VariableNames;
            if strcmp(var_name(1),"Speed in")
                this.EffTrqInMap(gear+1).setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                [spd_idx_out,trq_idx_out,EffMap_out] = this.fromInToOutMap(this.EffTrqInMap(gear+1).brkp1,this.EffTrqInMap(gear+1).brkp2,this.EffTrqInMap(gear+1).tab_data,this.gear_ratio(gear+1));
                var_name{1} = 'Speed out';
                var_name{2} = 'Torque out';
                this.EffTrqOutMap(gear+1).setLookup2DParameters(spd_idx_out,trq_idx_out,EffMap_out,var_name,var_info.Variables);
            elseif strcmp(var_name(1),"Speed out")
                this.EffTrqOutMap(gear+1).setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
            else
                error("Wrong axis selection")
            end
        end

        function plotEffMap(this,gear,gb_side,ContourLevels,LineWidth,Fontsize,FontName)
            % Method to plot the gearbox efficiency map. NOTE: gb_side is
            % the side on which to plot the map. Can be "Input" or "Output"
            if strcmp(gb_side,"Input")
                spd_idx = this.EffTrqInMap(gear+1).brkp1*30/pi;
                y_idx = this.EffTrqInMap(gear+1).brkp2;
                Map = this.EffTrqInMap(gear+1).tab_data;
            elseif strcmp(gb_side,"Output")
                spd_idx = this.EffTrqOutMap(gear+1).brkp1*30/pi;
                y_idx = this.EffTrqOutMap(gear+1).brkp2;
                Map = this.EffTrqOutMap(gear+1).tab_data;
            else
                error("Wrong side definition. Can be 'Input' or 'Output'")
            end
            contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{1},'ShowText','on');
            hold on
            for i = 2:length(ContourLevels)
                contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{i},'ShowText','on',"HandleVisibility","off");
            end
            colormap(flipud(gray(256)));
            clim([0.7 1.2]);
            hold off
            xlim([min(spd_idx) max(spd_idx)])
            ylim([min(y_idx) max(y_idx)])
            xlabel(gb_side+" Speed [rpm]");
            ylabel(gb_side+" Torque [Nm]");
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title("Efficiency map-"+gear+" gear",'fontsize',Fontsize);
        end
    end

    methods (Static)
        function [spd_idx_out,y_idx_out,Map_out] = fromInToOutMap(spd_idx_in,y_idx_in,Map_in,gear_ratio)
            % Method to convert the efficiency map from input size to
            % output size of the gearbox
            spd_idx_out = spd_idx_in/gear_ratio;
            y_idx_out = linspace(0,y_idx_in(end)*gear_ratio,length(y_idx_in))';
            [SpdIn,YIn] = meshgrid(spd_idx_in,y_idx_in);
            SpdOut = SpdIn/gear_ratio;
            YOut = YIn*gear_ratio.*Map_in;
            F = scatteredInterpolant(SpdOut(:),YOut(:),Map_in(:),'linear','linear');
            [SpdOut,YOut] = meshgrid(spd_idx_out,y_idx_out);
            Map_out = F(SpdOut,YOut);
        end

        function [spd_idx_out,y_idx_out,Map_out] = fromOutToInMap(spd_idx_out,y_idx_out,Map_out,gear_ratio)
            % Method to convert the efficiency map from input size to
            % output size of the gearbox
            spd_idx_out = spd_idx_out*gear_ratio;
            y_idx_out = linspace(0,y_idx_out(end)/gear_ratio,length(y_idx_out))';
            [SpdIn,YIn] = meshgrid(spd_idx_out,y_idx_out);
            SpdOut = SpdIn*gear_ratio;
            YOut = YIn/gear_ratio./Map_out;
            F = scatteredInterpolant(SpdOut(:),YOut(:),Map_out(:),'linear','linear');
            [SpdOut,YOut] = meshgrid(spd_idx_out,y_idx_out);
            Map_out = F(SpdOut,YOut);
        end

        function [spd_idx,pwr_idx,PwrMap] = fromTrqToPwrMap(spd_idx,trq_idx,pwr_idx,TrqMap)
            % Method to convert a torque based map into a power based map
            [Spd,Trq] = meshgrid(spd_idx,trq_idx);
            Pwr = Spd.*Trq;
            F = scatteredInterpolant(Spd(:),Pwr(:),TrqMap(:),'linear','linear');
            [Spd,Pwr] = meshgrid(spd_idx,pwr_idx);
            PwrMap = F(Spd,Pwr);
        end
    end
end