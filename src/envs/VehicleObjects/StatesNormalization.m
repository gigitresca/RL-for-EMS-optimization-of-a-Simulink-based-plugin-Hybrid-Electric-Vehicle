classdef StatesNormalization < matlab.mixin.Copyable
    properties
        states_max
        states_min
        upper
        lower
    end

    methods
        function this = StatesNormalization(states_max,states_min,up_value,low_value)
            % Constructor method: Intializes the normalization object with
            % the proper values
            if nargin>0
                this.states_max = states_max;
                this.states_min = states_min;
                this.upper = up_value;
                this.lower = low_value;
            end
        end

        function setNormalizationParameters(this,states_max,states_min,up_value,low_value)
            % Method to initializethe normalization object with
            % the proper values
            if iscolumn(states_max)
                this.states_max = states_max;
            else
                this.states_max = states_max';
            end
            if iscolumn(states_min)
                this.states_min = states_min;
            else
                this.states_min = states_min';
            end
            this.upper = up_value;
            this.lower = low_value;
        end

        function setNormalizationFromStates(this,env_obj,up_value,low_value)
            % Method to compute the max and min directly from states trend
            % input
            this.upper = up_value;
            this.lower = low_value;
            veh_spd = vertcat(env_obj.driving_cycles.veh_spd);
            veh_acc = vertcat(env_obj.driving_cycles.veh_acc);
            veh_gear = vertcat(env_obj.driving_cycles.veh_gear);
            grade = vertcat(env_obj.driving_cycles.grade);
            % Power demad computation
            wh_spd = veh_spd/env_obj.wh.radius;
            fd_spd = wh_spd*env_obj.fd.ratio;
            gb_in_spd = fd_spd.*env_obj.gb.gear_ratio(veh_gear+1);
            em_spd = gb_in_spd*env_obj.em.gear_ratio;
            ice_state = 1;
            veh_eq_mass = env_obj.veh.mass+(2*env_obj.wh.inertiaR+2*env_obj.wh.inertiaR)*1/(env_obj.wh.radius^2)+env_obj.ice.inertia*env_obj.gb.gear_ratio...
                (veh_gear+1).^2*env_obj.fd.ratio^2/(env_obj.wh.radius^2)*ice_state+... %Conventional Vehicle Inertia
                env_obj.em.inertia*env_obj.gb.gear_ratio(veh_gear+1).^2*env_obj.fd.ratio^2*env_obj.em.gear_ratio^2/(env_obj.wh.radius^2); %Additional BAS Inertia
            F_acc = veh_eq_mass.*veh_acc; %[N]
            F_res = (env_obj.veh.F0+env_obj.veh.mass*9.81*sin(grade*pi/180)+env_obj.veh.F1*(veh_spd*3.6)+env_obj.veh.F2*(veh_spd*3.6).^2); %[N]
            F_tot = F_res+F_acc; %[N]
            fd_trq = F_tot*(env_obj.wh.radius/env_obj.fd.ratio); %[Nm] Torque Request at the inlet of the Final Drive
            fd_pwr = fd_trq.*fd_spd; %[W] Power Request at the inlet of the Final Drive
            gb_eff = zeros(size(fd_trq));
            for i = 1:length(gb_eff)
                gb_eff(i) = interp2_lim(env_obj.gb.EffTrqOutMap(veh_gear(i)+1).brkp1,env_obj.gb.EffTrqOutMap(veh_gear(i)+1).brkp2,env_obj.gb.EffTrqOutMap(veh_gear(i)+1).tab_data,gb_in_spd(i),abs(fd_trq(i))); %[-] gb efficiency calculated depending on the output torque 
            end
            gb_in_pwr = fd_pwr.*gb_eff.^sign(-fd_pwr);
            pwr_dmd = gb_in_pwr;
            % Compute the maximum and minimum
            states_info = env_obj.getObservationInfo;
            states_info = states_info.Name;
            states = zeros(length(veh_spd),length(states_info));
            if nnz(strcmp(states_info,"veh_spd"))>0
                states(:,strcmp(states_info,"veh_spd")) = veh_spd;
            end
            if nnz(strcmp(states_info,"veh_acc"))>0
                states(:,strcmp(states_info,"veh_acc")) = veh_acc;
            end
            if nnz(strcmp(states_info,"em_spd"))>0
                states(:,strcmp(states_info,"em_spd")) = em_spd;
            end
            if nnz(strcmp(states_info,"gb_in_spd"))>0
                states(:,strcmp(states_info,"gb_in_spd")) = gb_in_spd;
            end
            if nnz(strcmp(states_info,"gb_in_pwr"))>0
                states(:,strcmp(states_info,"gb_in_pwr")) = gb_in_pwr;
            end
            if nnz(strcmp(states_info,"pwr_dmd"))>0
                states(:,strcmp(states_info,"pwr_dmd")) = pwr_dmd;
            end
            this.states_max = max(states,[],1)';
            this.states_min = min(states,[],1)';
            if nnz(strcmp(states_info,"batt_soc"))>0
                this.states_max(strcmp(states_info,"batt_soc")) = env_obj.ess.soc_high;
                this.states_min(strcmp(states_info,"batt_soc")) = env_obj.ess.soc_low;
            end
            if nnz(strcmp(states_info,"veh_dist_perc"))>0
                this.states_max(strcmp(states_info,"veh_dist_perc")) = 1;
                this.states_min(strcmp(states_info,"veh_dist_perc")) = 0;
            end
            if nnz(strcmp(states_info,"time_perc"))>0
                this.states_max(strcmp(states_info,"time_perc")) = 1;
                this.states_min(strcmp(states_info,"time_perc")) = 0;
            end
            if nnz(strcmp(states_info,"ice_state"))>0
                this.states_max(strcmp(states_info,"ice_state")) = 1;
                this.states_min(strcmp(states_info,"ice_state")) = 0;
            end
        end

        function obs = normalize(this,states)
            % Method to normalize an input vector and obtain the relative
            % normalized observation
            obs = (states-this.states_min)./(this.states_max-this.states_min)*(this.upper-this.lower)+this.lower;       
        end

        function states = statesFromObs(this,obs)
            % Method to obtain the unormalized states from the normalized
            % observation
            states = (obs-this.lower)/(this.upper-this.lower).*(this.states_max-this.states_min)+this.states_min;
        end
    end
end