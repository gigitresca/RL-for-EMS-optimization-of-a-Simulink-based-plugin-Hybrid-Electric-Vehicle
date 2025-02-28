classdef Engine < matlab.mixin.Copyable
    % object that define a generic Internal Combustion Engine described by
    % fuel consumption map and bsfc map

    properties
        fuel_heating_val
        fuel_density_val
        spd_idle
        spd_max
        bore
        stroke
        n_cyl
        inertia
        disp
        co2_mol_mass
        fuel_mol_mass
        cut_off_lim
    end

    properties (SetAccess = private)
        trq_max = Lookup1D
        pwr_max = Lookup1D
        trq_min = Lookup1D
        pwr_min = Lookup1D
        FuelPwrMap = Lookup2D
        BsfcPwrMap = Lookup2D
        FuelMap = Lookup2D
        BsfcMap = Lookup2D
        EffMap = Lookup2D
        trq_minbsfc = Lookup1D
        pwr_minbsfc = Lookup1D
    end

    methods
        function this = Engine(tab_data)
            % Constructor method: Intializes the scalar value of the engine
            % (the maps are set with the proper method)
            if nargin>0
                this.fuel_heating_val = tab_data.("Fuel LHV");
                this.fuel_density_val = tab_data.("Fuel density value");
                this.spd_idle = tab_data.("Idle speed");
                this.spd_max = tab_data.("Max speed");
                this.bore = tab_data.("Bore");
                this.stroke = tab_data.("Stroke");
                this.n_cyl = tab_data.("N cylinder");
                this.inertia = tab_data.("Inertia");
                if isnan(tab_data.("Displecement"))
                    this.disp = this.bore^2*pi/4*this.stroke*this.n_cyl/1000;
                else
                    this.disp = tab_data.("Displecement");
                end
                this.co2_mol_mass = tab_data.("CO2 molar mass");
                this.fuel_mol_mass = tab_data.("Fuel molar mass");
                this.cut_off_lim = tab_data.("Cut-off limit");
            end
        end

        function setEngineParameteres(this,tab_data)
            % Method to set the scalar value of the engine
            % (the maps are set with the proper method)
            this.fuel_heating_val = tab_data.("Fuel LHV");
                this.fuel_density_val = tab_data.("Fuel density value");
                this.spd_idle = tab_data.("Idle speed");
                this.spd_max = tab_data.("Max speed");
                this.bore = tab_data.("Bore");
                this.stroke = tab_data.("Stroke");
                this.n_cyl = tab_data.("N cylinder");
                this.inertia = tab_data.("Inertia");
                if isnan(tab_data.("Displecement"))
                    this.disp = this.bore^2*pi/4*this.stroke*this.n_cyl/1000;
                else
                    this.disp = tab_data.("Displecement");
                end
                this.co2_mol_mass = tab_data.("CO2 molar mass");
                this.fuel_mol_mass = tab_data.("Fuel molar mass");
                this.cut_off_lim = tab_data.("Cut-off limit");
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
            % Method to intialize the engine maps depending on the input
            var_name = var_info.Properties.VariableNames;
            if strcmp(var_name(3),"Fuel rate")
                if strcmp(var_name(2),"Torque")
                    this.FuelMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                elseif strcmp(var_name(2),"Power")
                    this.FuelPwrMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                else
                    error("Wrong y axis selection")
                end
            elseif strcmp(var_name(3),"BSFC")
                if strcmp(var_name(2),"Torque")
                    this.BsfcMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                elseif strcmp(var_name(2),"Power")
                    
                    this.BsfcPwrMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                else
                    error("Wrong y axis selection")
                end
            elseif strcmp(var_name(3),"Efficiency")
                if strcmp(var_name(2),"Torque")
                    this.EffMap.setLookup2DParameters(tab_data(1,2:end),tab_data(2:end,1),tab_data(2:end,2:end),var_name,var_info.Variables);
                else
                    error("Wrong y axis selection")
                end
            else
                error("Invalid map")
            end
        end

        function bmep = bmepFromTrq(this,trq,i)
            % Method to compute the ice bmep [bar] from torque [Nm] value
            bmep = (trq*2*pi*i)/(this.disp*10^(-6))*10^(-5);
        end

        function bmep = bmepFromPwr(this,spd,pwr,i)
            % Method to compute the ice bmep [bar] from speed [rad/s] and
            % power [W] values
            bmep = (60*i*pwr)./((spd*30/pi)*(this.disp*10^(-6)))*10^(-5);
        end

        function pwr = pwrFromBmep(this,spd,bmep,i)
            % Method to compute the ice pwr [W] from speed [rad/s] and
            % bmep [bar] values
            pwr = ((bmep*10^5)*spd*(30/pi)*(this.disp*10^-6))/(60*i);
        end

        function minBsfcComputation(this,y_axis)
            % Method to compute the Optimal Operating Line (OOL). NOTE: y_axis input requires a
            % string that can be "Torque","Bmep","Power"
            if strcmp(y_axis,"Torque")
                spd_idx = this.BsfcMap.brkp1;
                y_idx = this.BsfcMap.brkp2;
                fl_spd_idx = this.trq_max.brkp1;
                y_fl = this.trq_max.tab_data;
                Map = this.BsfcMap.tab_data;
                var_names = this.trq_max.var_names;
                var_units = this.trq_max.var_units;
            elseif strcmp(y_axis,"Power")
                spd_idx = this.BsfcPwrMap.brkp1;
                y_idx = this.BsfcPwrMap.brkp2;
                fl_spd_idx = this.pwr_max.brkp1;
                y_fl = this.pwr_max.tab_data;
                Map = this.BsfcPwrMap.tab_data;
                var_names = this.pwr_max.var_names;
                var_units = this.pwr_max.var_units;
            else
                error('Invalid x_axis setting: choose between "Torque", "Bmep", "Power"')
            end
            y_fl_idx = interp1_lim(fl_spd_idx,y_fl,spd_idx');
            [SpdMtx,YMtx] = meshgrid(spd_idx,y_idx);
            Map(YMtx>y_fl_idx') = NaN;
            minBSFC = min(Map);
            minBSFC = repmat(minBSFC,[length(y_idx) 1]);
            minBSFC( isnan(minBSFC) ) = 10^6;
            spd_minbsfc = SpdMtx(Map==minBSFC);
            if strcmp(y_axis,"Torque")
                y_minbsfc = YMtx(Map==minBSFC);
                this.trq_minbsfc.setLookup1DParameters(spd_minbsfc,y_minbsfc,var_names,var_units);
            elseif strcmp(y_axis,"Power")
                y_minbsfc = YMtx(Map==minBSFC);
                this.pwr_minbsfc.setLookup1DParameters(spd_minbsfc,y_minbsfc,var_names,var_units);
            end
        end

        function [YIdxTD,SpdIdxTD,TimeDistribution] = timeDistribution(this,Map,y_axis,time,ice_y,ice_spd,ice_state)
            % Method to compute the time distribution starting from ice
            % trq, spd and state signal. It is possible to select Bsfc or
            % Fuel map on which compute the time distribution and different
            % y axis between Bmep, Torque and Power
            if strcmp(Map,"Bsfc")
                if strcmp(y_axis,"Torque")
                    spd_idx = this.BsfcMap.brkp1;
                    y_idx = this.BsfcMap.brkp2;
                    Map = this.BsfcMap.tab_data;
                elseif strcmp(y_axis,"Bmep")
                    spd_idx = this.BsfcMap.brkp1;
                    y_idx = this.bmepFromTrq(this.BsfcMap.brkp2,2);
                elseif strcmp(y_axis,"Power")
                    spd_idx = this.BsfcPwrMap.brkp1;
                    y_idx = this.BsfcPwrMap.brkp2;
                else
                    error('Invalid x_axis setting: choose between "Torque", "Bmep", "Power"')
                end
            elseif strcmp(Map,"Fuel")
                if strcmp(y_axis,"Torque")
                    spd_idx = this.FuelMap.brkp1;
                    y_idx = this.FuelMap.brkp2;
                    Map = this.FuelMap.tab_data;
                elseif strcmp(y_axis,"Bmep")
                    spd_idx = this.FuelMap.brkp1;
                    y_idx = this.bmepFromTrq(this.FuelMap.brkp2,2);
                elseif strcmp(y_axis,"Power")
                    spd_idx = this.FuelPwrMap.brkp1;
                    y_idx = this.FuelPwrMap.brkp2;
                else
                    error('Invalid x_axis setting: choose between "Torque", "Bmep", "Power"')
                end
            else
                error('Invalid map selection: choose between "Bsfc", "Fuel"')
            end
            if iscategorical(ice_state)
                ice_state = double(ice_state);
                ice_state(ice_state==1) = 0;
                ice_state(ice_state==2) = 1;
            end
            n_spd = 26;
            n_y = 26;
            spd_idx_td = linspace(0,spd_idx(end),n_spd);
            y_idx_td = linspace(y_idx(1),y_idx(end),n_y);
            [SpdIdxTD,YIdxTD] = meshgrid(spd_idx_td,y_idx_td);
            TimeDistribution = zeros(size(YIdxTD));
                    
            for i = 1:length(time)
                for j = 1:size(YIdxTD,1)-1
                    if ice_y(i)==YIdxTD(1,1)
                        row = j;
                        break
                    elseif ice_y(i)>YIdxTD(j,1) && ice_y(i)<=YIdxTD(j+1,1)
                        row = j+1;
                        break
                    end
                end
                for j = 1:size(SpdIdxTD,2)-1
                    if ice_spd(i)==SpdIdxTD(1,1)
                        col = j;
                        break
                    elseif ice_spd(i)>SpdIdxTD(1,j) && ice_spd(i)<=SpdIdxTD(1,j+1)
                        col = j+1;
                        break
                    end
                end
                if i == 1
                    TimeDistribution(row,col) = TimeDistribution(row,col)+time(i);
                else
                    TimeDistribution(row,col) = TimeDistribution(row,col)+(time(i)-time(i-1))*ice_state(i);
                end
            end
            TimeDistribution = TimeDistribution/time(end)*100;
        end

        function plotBsfcMap(this,y_axis,ContourLevels,LineWidth,Fontsize,FontName)
            % Method to plot the BSFC map. NOTE: y_axis input requires a
            % string that can be "Torque","Bmep","Power"
            if strcmp(y_axis,"Torque")
                spd_idx = this.BsfcMap.brkp1*30/pi;
                y_idx = this.BsfcMap.brkp2;
                fl_spd_idx = this.trq_max.brkp1*30/pi;
                y_fl = this.trq_max.tab_data;
                spd_ool = this.trq_minbsfc.brkp1*30/pi;
                y_ool = this.trq_minbsfc.tab_data;
                Map = this.BsfcMap.tab_data;
            elseif strcmp(y_axis,"Bmep")
                spd_idx = this.BsfcMap.brkp1*30/pi;
                y_idx = this.bmepFromTrq(this.BsfcMap.brkp2,2);
                fl_spd_idx = this.trq_max.brkp1*30/pi;
                y_fl = this.bmepFromTrq(this.trq_max.tab_data,2);
                spd_ool = this.trq_minbsfc.brkp1*30/pi;
                y_ool = this.bmepFromTrq(this.trq_minbsfc.tab_data,2);
                Map = this.BsfcMap.tab_data;
            elseif strcmp(y_axis,"Power")
                spd_idx = this.BsfcPwrMap.brkp1*30/pi;
                y_idx = this.BsfcPwrMap.brkp2;
                fl_spd_idx = this.pwr_max.brkp1*30/pi;
                y_fl = this.pwr_max.tab_data;
                spd_ool = this.pwr_minbsfc.brkp1*30/pi;
                y_ool = this.pwr_minbsfc.tab_data;
                Map = this.BsfcPwrMap.tab_data;
            else
                error('Invalid x_axis setting: choose between "Torque", "Bmep", "Power"')
            end
            y_fl_idx = interp1_lim(fl_spd_idx,y_fl,spd_idx');
            [~,YMtx] = meshgrid(spd_idx,y_idx);
            Map(YMtx>y_fl_idx') = NaN;
            contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{1},'ShowText','on');
            hold on
            for i = 2:length(ContourLevels)
                contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{i},'ShowText','on',"HandleVisibility","off");
            end
            plot(fl_spd_idx,y_fl,'ro-.','LineWidth',LineWidth);
            plot(spd_ool,y_ool,'bo-.','LineWidth',LineWidth);
            colormap(flipud(gray(256)));
            clim([200 500]);
            hold off
            xlim([1000 5000])
            ylim([0 max(y_fl_idx)])
            xlabel('Engine Speed [rpm]');
            if strcmp(y_axis,"Torque")
                ylabel('Torque [Nm]');
            elseif strcmp(y_axis,"Bmep")
                ylabel('Bmep [bar]');
            elseif strcmp(y_axis,"Power")
                ylabel('Power [kW]');
            end
            legend('BSFC [g/kWh]','Full load','OOL','FontSize',Fontsize);
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title('BSFC','fontsize',Fontsize);
        end

        function plotFuelMap(this,y_axis,ContourLevels,LineWidth,Fontsize,FontName)
            % Method to plot the fuel consumption map. NOTE: y_axis input requires a
            % string that can be "Torque","Bmep","Power"
            if strcmp(y_axis,"Torque")
                spd_idx = this.FuelMap.brkp1*30/pi;
                y_idx = this.FuelMap.brkp2;
                fl_spd_idx = this.trq_max.brkp1*30/pi;
                y_fl = this.trq_max.tab_data;
                spd_ool = this.trq_minbsfc.brkp1*30/pi;
                y_ool = this.trq_minbsfc.tab_data;
                Map = this.FuelMap.tab_data*3600;
            elseif strcmp(y_axis,"Bmep")
                spd_idx = this.FuelMap.brkp1*30/pi;
                y_idx = this.bmepFromTrq(this.FuelMap.brkp2,2);
                fl_spd_idx = this.trq_max.brkp1*30/pi;
                y_fl = this.bmepFromTrq(this.trq_max.tab_data,2);
                spd_ool = this.trq_minbsfc.brkp1*30/pi;
                y_ool = this.bmepFromTrq(this.trq_minbsfc.tab_data,2);
                Map = this.FuelMap.tab_data*3600;
            elseif strcmp(y_axis,"Power")
                spd_idx = this.FuelPwrMap.brkp1*30/pi;
                y_idx = this.FuelPwrMap.brkp2;
                fl_spd_idx = this.pwr_max.brkp1*30/pi;
                y_fl = this.pwr_max.tab_data;
                spd_ool = this.pwr_minbsfc.brkp1*30/pi;
                y_ool = this.pwr_minbsfc.tab_data;
                Map = this.FuelPwrMap.tab_data*3600;
            else
                error('Invalid x_axis setting: choose between "Torque", "Bmep", "Power"')
            end
            y_fl_idx = interp1_lim(fl_spd_idx,y_fl,spd_idx');
            [~,YMtx] = meshgrid(spd_idx,y_idx);
            Map(YMtx>y_fl_idx') = NaN;
            contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{1},'ShowText','on');
            hold on
            for i = 2:length(ContourLevels)
                contourf(spd_idx,y_idx,Map,'LineWidth',LineWidth,'LevelList',ContourLevels{i},'ShowText','on','HandleVisibility','off');
            end
            plot(fl_spd_idx,y_fl,'ro-.','LineWidth',LineWidth);
            plot(spd_ool,y_ool,'bo-.','LineWidth',LineWidth);
            colormap(flipud(gray(256)));
            clim([0 40]);
            hold off
            xlim([1000 5000])
            ylim([0 max(y_fl_idx)])
            xlabel('Engine Speed [rpm]');
            if strcmp(y_axis,"Torque")
                ylabel('Torque [Nm]');
            elseif strcmp(y_axis,"Bmep")
                ylabel('Bmep [bar]');
            elseif strcmp(y_axis,"Power")
                ylabel('Power [kW]');
            end
            legend('Fuel consumption [kg/h]','Full Load','OOL','FontSize',Fontsize);
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title('Fuel consumption','fontsize',Fontsize);
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