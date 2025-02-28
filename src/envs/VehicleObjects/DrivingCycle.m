classdef DrivingCycle < matlab.mixin.Copyable
    properties
        name
        time_simu
        veh_spd
        veh_acc
        veh_gear
        veh_dist
        grade
        Ts
    end

    methods
        function this = DrivingCycle(tab_data,name)
            % Constructor method: initializes the parameters of the driving
            % cycle object
            if nargin>0
                if iscell(name) || ischar(name)
                    convertCharsToStrings(name)
                end
                this.name = name;
                this.time_simu = tab_data.Time;
                this.time_simu(isnan(this.time_simu)) = [];
                this.veh_spd = tab_data.("Vehicle Speed")/3.6;
                this.veh_spd(isnan(this.veh_spd)) = [];
                this.veh_gear = tab_data.("Gear Number");
                this.veh_gear(isnan(this.veh_gear)) = [];
                if nnz(strcmp("Grade",tab_data.Properties.VariableNames))>0
                    this.grade = tab_data.("Grade");
                    this.grade(isnan(this.grade)) = [];
                else
                    this.grade = zeros(size(this.time_simu));
                end
                this.veh_acc = [0;diff(this.veh_spd)./diff(this.time_simu)];
                this.veh_dist = cumtrapz(this.time_simu,this.veh_spd);
                this.Ts = this.time_simu(2)-this.time_simu(1);
            end
        end

        function setDrivingCycleParameters(this,tab_data)
            % Method to initializes the parameters of the driving
            % cycle object
            if iscell(name) || ischar(name)
                convertCharsToStrings(name)
            end
            this.name = name;
            this.time_simu = tab_data.Time;
            this.time_simu(isnan(this.time_simu)) = [];
            this.veh_spd = tab_data.("Vehicle Speed")/3.6;
            this.veh_spd(isnan(this.veh_spd)) = [];
            this.veh_gear = tab_data.("Gear Number");
            this.veh_gear(isnan(this.veh_gear)) = [];
            if nnz(strcmp("Grade",tab_data.Properties.VariableNames))>0
                this.grade = tab_data.("Grade");
                this.grade(isnan(this.grade)) = [];
            else
                this.grade = zeros(size(this.time_simu));
            end
            this.veh_acc = [0;diff(this.veh_spd)./diff(this.time_simu)];
            this.veh_dist = cumtrapz(this.time_simu,this.veh_spd);
            this.Ts = this.time_simu(2)-this.time_simu(1);
        end

        function resample(this,Ts)
            % Method to resample the driving cycle
            time_simu_new = (this.time_simu(1):Ts:this.time_simu(end))';
            veh_spd_new = interp1_lim(this.time_simu,this.veh_spd,time_simu_new);
            veh_gear_new = round(interp1_lim(this.time_simu,this.veh_gear,time_simu_new));
            grade_new = interp1_lim(this.time_simu,this.grade,time_simu_new);
            veh_acc_new = [0;diff(veh_spd_new)./diff(time_simu_new)];
            veh_dist_new = cumtrapz(time_simu_new,veh_spd_new);
            % Assign the new interpolated vector to object properties
            this.Ts = Ts;
            this.time_simu = time_simu_new;
            this.veh_spd = veh_spd_new;
            this.veh_gear = veh_gear_new;
            this.grade = grade_new;
            this.veh_acc = veh_acc_new;
            this.veh_dist = veh_dist_new;
        end

        function plotSpeedProfile(this,LineWidth,Fontsize,FontName)
            plot(this.time_simu,this.veh_spd*3.6,"LineWidth",LineWidth);
            xlabel('Time [s]')
            ylabel('Vehicle Speed [km/h]')
            xlim([0 this.time_simu(end)])
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title(this.name+" speed profile")
            grid on
        end

        function plotAccelerationProfile(this,LineWidth,Fontsize,FontName)
            plot(this.time_simu,this.veh_acc,"LineWidth",LineWidth);
            xlabel('Time [s]')
            ylabel('Vehicle Acceleration [$m/s^2$]')
            xlim([0 this.time_simu(end)])
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title(this.name+" acceleration profile")
            grid on
        end

        function plotGearProfile(this,LineWidth,Fontsize,FontName)
            plot(this.time_simu,this.veh_gear,"LineWidth",LineWidth);
            xlabel('Time [s]')
            ylabel('Gear [-]')
            xlim([0 this.time_simu(end)])
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title(this.name+" gear profile")
            grid on
        end

        function plotGradeProfile(this,LineWidth,Fontsize,FontName)
            plot(this.time_simu,this.grade,"LineWidth",LineWidth);
            xlabel('Time [s]')
            ylabel('Grade [-]')
            xlim([0 this.time_simu(end)])
            set(gca,'Fontsize',Fontsize,'FontName',FontName);
            title(this.name+" grade profile")
            grid on
        end
    end
end