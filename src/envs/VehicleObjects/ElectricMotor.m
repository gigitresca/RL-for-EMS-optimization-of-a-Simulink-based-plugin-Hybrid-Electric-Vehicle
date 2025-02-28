classdef ElectricMotor < matlab.mixin.Copyable
    % object that define a generic Internal Combustion Engine described by
    % fuel consumption map and bsfc map

    properties
        mass
        inertia
        gear_ratio 
    end

    properties (SetAccess = private)
        trq_max = Lookup1D
        pwr_max = Lookup1D
        trq_min = Lookup1D
        pwr_min = Lookup1D
        EffPwrMap = Lookup2D
        RegMap = Lookup2D
        pwr_maxeff = Lookup1D
    end

    methods
        function this = ElectricMotor(tab_data)
            % Constructor method: Intializes the scalar value of the
            % electric motor
            % (the maps are set with the proper method)
            if nargin>0
                this.mass = tab_data.("Mass");
                this.inertia = tab_data.("Inertia");
                this.gear_ratio = tab_data.("Gear ratio");
            end
        end

        function setElectricMotorParameteres(this,tab_data)
            % Method to set the scalar value of the electric motor
            % (the maps are set with the proper method)
            this.mass = tab_data.("Mass");
            this.inertia = tab_data.("Inertia");
            this.gear_ratio = tab_data.("Gear ratio");
        end

        function setLimits(this,tab_data,var_info)
            % Method to initialize the maximum and minimum torque and power
            % limits
            this.trq_max.setLookup1DParameters(tab_data.Speed,tab_data.("Max Torque"),{'Speed','Torque'},[var_info.Speed var_info.('Max Torque')]);
            this.pwr_max.setLookup1DParameters(tab_data.Speed,tab_data.("Max Power"),{"Speed",'Power'},[var_info.Speed var_info.('Max Power')]);
            this.trq_min.setLookup1DParameters(tab_data.Speed,tab_data.("Min Torque"),{'Speed','Torque'},[var_info.Speed var_info.('Min Torque')]);
            this.pwr_min.setLookup1DParameters(tab_data.Speed,tab_data.("Min Power"),{'Speed','Power'},[var_info.Speed var_info.('Min Power')]);
        end

        function setMap(this,tab_data,var_info)
            % Method to initialize the electric motor map based on the
            % input
            var_name = var_info.Properties.VariableNames;
            if strcmp(var_name(3),"Efficiency")
                if strcmp(var_name(2),"Power")
                    this.EffPwrMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                else
                    error("Wrong y axis selection");
                end
            elseif strcmp(var_name(3),"Regeneration ratio")
                if strcmp(var_name(1),"Acceleration") && strcmp(var_name(2),"Speed")
                    this.RegMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                else
                    error("Wrong axis selection");
                end
            else
                error("Wrong map")
            end
        end
      
        function maxEffComputation(this)
            % Method to compute the Optimal Operating Line (OOL).
            spd_idx = this.EffPwrMap.brkp1;
            y_idx = this.EffPwrMap.brkp2;
            fl_spd_idx = this.pwr_max.brkp1;
            y_fl = this.pwr_max.tab_data;
            Map = this.EffPwrMap.tab_data;
            var_names = this.pwr_max.var_names;
            var_units = this.pwr_max.var_units;
            y_fl_idx = interp1_lim(fl_spd_idx,y_fl,spd_idx');
            [SpdMtx,YMtx] = meshgrid(spd_idx,y_idx);
            Map(YMtx>y_fl_idx') = NaN;
            Map(YMtx<-y_fl_idx') = NaN;
            Map(Map==1) = NaN;
            maxEff = max(Map);
            maxEff = repmat(maxEff,[length(y_idx) 1]);
            maxEff( isnan(maxEff) ) = 0;
            spd_maxeff = SpdMtx(Map==maxEff);
            y_maxeff = YMtx(Map==maxEff);
            this.pwr_maxeff.setLookup1DParameters(spd_maxeff,y_maxeff,var_names,var_units);
        end

        function plotEffMap(this,ContourLevels,LineWidth,Fontsize,FontName)
            % Method to plot the EM efficiency map
            spd_idx = this.EffPwrMap.brkp1*30/pi;
            y_idx = this.EffPwrMap.brkp2/1000;
            spd_fl_trac_idx = this.pwr_max.brkp1*30/pi;
            y_fl_trac = this.pwr_max.tab_data/1000;
            spd_fl_brak_idx = this.pwr_min.brkp1*30/pi;
            y_fl_brak = this.pwr_min.tab_data/1000;
            spd_ool = this.pwr_maxeff.brkp1*30/pi;
            y_ool = this.pwr_maxeff.tab_data/1000;
            Map = this.EffPwrMap.tab_data;
            y_fl_trac_idx = interp1_lim(spd_fl_trac_idx,y_fl_trac,spd_idx');
            y_fl_brak_idx = interp1_lim(spd_fl_brak_idx,y_fl_brak,spd_idx');
            [~,YMtx] = meshgrid(spd_idx,y_idx);
            Map(YMtx>y_fl_trac_idx') = NaN;
            Map(YMtx<y_fl_brak_idx') = NaN;
            contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{1},'ShowText','on');
            hold on
            for i = 2:length(ContourLevels)
                contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{i},'ShowText','on',"HandleVisibility","off");
            end
            plot(spd_fl_trac_idx,y_fl_trac,'ro-.','LineWidth',LineWidth);
            plot(spd_fl_brak_idx,y_fl_brak,'ro-.','LineWidth',LineWidth,'HandleVisibility','off');
            plot(spd_ool,y_ool,'bo-.','LineWidth',LineWidth);
            plot(spd_ool,-y_ool,'bo-.','LineWidth',LineWidth,'HandleVisibility','off');
            colormap(flipud(gray(256)));
            clim([0.7 1.2]);
            hold off
            xlim([500 5000])
            ylim([min(y_fl_brak) max(y_fl_trac)])
            xlabel('EM Speed [rpm]');
            ylabel('Power [kW]');
            legend('$\eta$','Full load','OOL','FontSize',Fontsize);
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title('Efficiency map','fontsize',Fontsize);
        end

        function plotRegMap(this,ContourLevels,LineWidth,Fontsize,FontName)
            % Method to plot the EM regeneration map
            acc_idx = this.RegMap.brkp1;
            spd_idx = this.RegMap.brkp2*3.6;
            Map = this.RegMap.tab_data;
            contourf(acc_idx,spd_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{1},'ShowText','on');
            hold on
            for i = 2:length(ContourLevels)
                contourf(acc_idx,spd_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{i},'ShowText','on',"HandleVisibility","off");
            end
            colormap(flipud(gray(256)));
            clim([0 2]);
            hold off
            xlim([min(acc_idx) 0])
            ylim([0 max(spd_idx)])
            xlabel('Acceleration [$m/s^2$]');
            ylabel('Speed [km/h]');
            legend('$P_{EM}$/$P{gb}$','FontSize',Fontsize);
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title('Regeneration map','fontsize',Fontsize);
        end
    end

    methods (Static)
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