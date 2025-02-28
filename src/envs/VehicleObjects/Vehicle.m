classdef Vehicle < matlab.mixin.Copyable
    % object that define the longitudinal dynamics of a generic vehicle
    % through coast-down approach

    properties
        axal_base
        cg_height
        cargo_mass
        body_mass
        mass
        F0
        F1
        F2
    end

    methods
        function this = Vehicle(tab_data)
            % Constructor method: Intializes the scalar value of the
            % vehicle
            if nargin>0
                this.axal_base = tab_data.("Axal base");
                this.cg_height = tab_data.("CG height");
                this.cargo_mass = tab_data.("Cargo mass");
                this.body_mass = tab_data.("Body mass");
                if isnan(tab_data.("Mass"))
                    this.mass = this.body_mass+this.cargo_mass;
                else
                    this.mass = tab_data.("Mass");
                end
                this.F0 = tab_data.("F0");
                this.F1 = tab_data.("F1");
                this.F2 = tab_data.("F2");
            end
        end

        function setVehicleParameters(this,tab_data)
            % Method to set the vehicle parameters
            this.axal_base = tab_data.("Axal base");
            this.cg_height = tab_data.("CG height");
            this.cargo_mass = tab_data.("Cargo mass");
            this.body_mass = tab_data.("Body mass");
            if isnan(tab_data.("Mass"))
                this.mass = this.body_mass+this.cargo_mass;
            else
                this.mass = tab_data.("Mass");
            end
            this.F0 = tab_data.("F0");
            this.F1 = tab_data.("F1");
            this.F2 = tab_data.("F2");
        end

        function plotRoadReistance(this,LineWidth,Fontsize,FontName)
            % Method to plot road resistance based on coast-down
            % coefficient approach
            spd_idx = linspace(0,150,50);
            force_idx = this.F0+this.F1*spd_idx+this.F2*spd_idx.^2;
            pwr_idx = force_idx.*(spd_idx/3.6);
            yyaxis left
            plot(spd_idx,force_idx,"LineWidth",LineWidth);
            ylabel('Force [N]')
            yyaxis right
            plot(spd_idx,pwr_idx/1000,"LineWidth",LineWidth);
            ylabel('Power [kW]')
            xlabel('Speed [km/h]')
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title('Road resistance')
            legend('Force','Power','Location','northwest')
            grid on
        end
    end
end