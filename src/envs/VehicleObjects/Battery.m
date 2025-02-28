classdef Battery < matlab.mixin.Copyable
    % object that define a generic Internal Combustion Engine described by
    % fuel consumption map and bsfc map

    properties
        num_cell_series
        num_module_parallel
        num_cell
        volt_nom
        volt_max
        volt_min
        mass
        mass_cell
        soc_min
        soc_max
        soc_high
        soc_low
        soc_init
        soc_trg
        cell_cap
        cell_curr_max_chg
        cell_curr_max_dis
        accelec
    end

    properties (SetAccess = private)
        cell_voc = Lookup1D
        cell_rint_chg = Lookup1D
        cell_rint_dis = Lookup1D
    end

    methods
        function this = Battery(tab_data)
            % Constructor method: Intializes the scalar value of the
            % battery
            % (the maps are set with the proper method)
            if nargin>0
                this.num_cell_series = tab_data.("N cell series");
                this.num_module_parallel = tab_data.("N module parallel");
                this.num_cell = this.num_cell_series*this.num_module_parallel;
                this.volt_nom = tab_data.("Nominal voltage");
                this.volt_max = tab_data.("Max voltage");
                this.volt_min = tab_data.("Min voltage");
                this.mass = tab_data.("Mass");
                this.mass_cell = tab_data.("Mass cell");
                this.soc_min = tab_data.("Min SoC");
                this.soc_max = tab_data.("Max SoC");
                this.soc_high = tab_data.("High SoC");
                this.soc_low = tab_data.("Low SoC");
                this.cell_cap = tab_data.("Cell capacity");
                this.cell_curr_max_chg = tab_data.("Cell max current charge");
                this.cell_curr_max_dis = tab_data.("Cell max current discharge");
                this.accelec = tab_data.("Electric accessory");
            end
        end

        function setBatteryParameteres(this,tab_data)
            % Method to set the scalar value of the battery
            % (the maps are set with the proper method)
            this.num_cell_series = tab_data.("N cell series");
            this.num_module_parallel = tab_data.("N module parallel");
            this.num_cell = this.num_cell_series*this.num_module_parallel;
            this.volt_nom = tab_data.("Nominal voltage");
            this.volt_max = tab_data.("Max voltage");
            this.volt_min = tab_data.("Min voltage");
            this.mass = tab_data.("Mass");
            this.mass_cell = tab_data.("Mass cell");
            this.soc_min = tab_data.("Min SoC");
            this.soc_max = tab_data.("Max SoC");
            this.soc_high = tab_data.("High SoC");
            this.soc_low = tab_data.("Low SoC");
            this.cell_cap = tab_data.("Cell capacity");
            this.cell_curr_max_chg = tab_data.("Cell max current charge");
            this.cell_curr_max_dis = tab_data.("Cell max current discharge");
            this.accelec = tab_data.("Electric accessory");
        end

        function setBatteryStrategy(this,soc_init,soc_trg)
            % Method to set the intial and target value of battery state of
            % charge
            this.soc_init = soc_init;
            this.soc_trg = soc_trg;
        end

        function setCellData(this,tab_data,var_info)
            % Method to initialize cell open circuit voltage [V]and internal
            % resistance [Ohm]
            var_name = var_info.Properties.VariableNames;
            this.cell_voc.setLookup1DParameters(tab_data.("State of charge"),tab_data.("Cell open circuit voltage"),{var_name{1},var_name{2}},[var_info.("State of charge") var_info.("Cell open circuit voltage")]);
            this.cell_rint_chg.setLookup1DParameters(tab_data.("State of charge"),tab_data.("Cell internal resistance charge"),{var_name{1},var_name{3}},[var_info.("State of charge") var_info.("Cell internal resistance charge")]);
            this.cell_rint_dis.setLookup1DParameters(tab_data.("State of charge"),tab_data.("Cell internal resistance discharge"),{var_name{1},var_name{4}},[var_info.("State of charge") var_info.("Cell internal resistance discharge")]);
        end

        function plotCellData(this,LineWidth,Fontsize,FontName)
            % Method to plot the open circuit voltage and internal
            % resistance of the cell
            yyaxis left
            plot(this.cell_voc.brkp1,this.cell_voc.tab_data,"LineWidth",LineWidth);
            ylabel('Open Circuit Voltage [V]')
            yyaxis right
            plot(this.cell_rint_chg.brkp1,this.cell_rint_chg.tab_data*1000,"LineWidth",LineWidth);
            hold on
            plot(this.cell_rint_dis.brkp1,this.cell_rint_dis.tab_data*1000,"LineWidth",LineWidth);
            hold off
            ylabel('Internal Resistance [$m\Omega$]')
            xlabel('State of Charge [-]')
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title('Battery Cell Feature')
            legend('$V_{oc}$','$R_{int,chg}$','$R_{int,dis}$')
            grid on
        end
    end
end