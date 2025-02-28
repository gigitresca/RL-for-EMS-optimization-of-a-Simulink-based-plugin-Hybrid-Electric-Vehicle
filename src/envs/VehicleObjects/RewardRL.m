classdef RewardRL < matlab.mixin.Copyable
    properties
        reward_type
        k1
        k2
        k3
        eff_dis
        eff_chg
        soc_gain
        soc_trg
        soc_bc_gain
        unfeaspenalty
        soc_bc = Lookup1D
    end

    methods
        function this = RewardRL(param)
            % Constructor method: Intializes the values of the reward object
            % (the maps are set with the proper method)
            arguments (Input)
                param.k1 double = 16
                param.k2 double = 0.01
                param.k3 double = 6000
                param.soc_bc_gain double = []
                param.unfeaspenalty double = 0
                param.reward_type string = "Dense"
                param.ice_fuel_heating_val double = []
                param.ess_cell_cap double = [];
                param.ess_volt_nom double = [];
                param.ess_num_cell_series double = [];
                param.ess_soc_trg double = [];
            end
            ice_info.fuel_heating_val = param.ice_fuel_heating_val;
            ess_info.cell_cap = param.ess_cell_cap;
            ess_info.volt_nom = param.ess_volt_nom;
            ess_info.num_cell_series = param.ess_num_cell_series;
            ess_info.soc_trg = param.ess_soc_trg;
            this.setRewardParameters( ...
                param.k1, ...
                param.k2, ...
                param.k3, ...
                param.soc_bc_gain, ...
                param.unfeaspenalty, ...
                param.reward_type, ...
                ice_info,ess_info ...
            );
        end
        function this = setRewardParameters(this,k1,k2,k3,soc_bc_gain,unfeaspenalty,reward_type,ice_info,ess_info)
            % Method to set the reward parameteres (sclaing factors ans
            % reward equation
            this.reward_type = reward_type;
            this.eff_dis = 0.95;
            this.eff_chg = 0.95;
            this.soc_gain = (1/this.eff_chg*ess_info.cell_cap*3600*ess_info.volt_nom*ess_info.num_cell_series/ice_info.fuel_heating_val); %[kg di combustibile equivalenti]
            this.soc_trg = ess_info.soc_trg;
            this.k1 = k1;
            this.k2 = k2;
            this.k3 = k3;
            this.soc_bc_gain = soc_bc_gain;
            this.unfeaspenalty = unfeaspenalty;
            if ~(strcmp(reward_type,"Dense") || strcmp(reward_type,"Sparse") || strcmp(reward_type,"ECMS") || strcmp(reward_type,"SOC_BC") || strcmp(reward_type,"Dense_multicycle"))
                error("The selected reward type is unavailable")
            end
        end
        function reward = getReward(this,ice_fuel_rate,batt_soc,batt_pwr,time,veh_dist_tot,veh_dist_perc,cycle_name,Ts,isdone,unfeasible)
            % Method to compute the reward from ice fuel rate and soc (or
            % battery power)
            fuel = 0;
            eq_fuel = 0;
            if strcmp(this.reward_type,"Dense")
                fuel = this.k2*ice_fuel_rate*(Ts);
                eq_fuel = this.k3*(this.soc_gain)*(batt_soc-this.soc_trg).^2;
            elseif strcmp(this.reward_type,"Sparse")
                fuel = this.k2*ice_fuel_rate*(Ts);
                eq_fuel = this.k3*(this.soc_gain)*(batt_soc-this.soc_trg).^2.*isdone;
            elseif strcmp(this.reward_type,"ECMS")
                fuel = this.k2*ice_fuel_rate*(Ts);
                eq_fuel = this.k3*batt_pwr.*(1/this.eff_dis*(batt_pwr>=0)+1*this.eff_chg*(batt_pwr<0))/(ice_info.fuel_heating_val)*Ts;
            elseif strcmp(this.reward_type,"SOC_BC")
                soc_up = interp1_lim(this.soc_bc.brkp1,this.soc_bc.tab_data,time);
                soc_low = this.soc_trg-(soc_up-this.soc_trg);
                out_range = (batt_soc>=soc_up)||(batt_soc<=soc_low);
                fuel = this.k2*ice_fuel_rate*(Ts);
                eq_fuel = this.k3*(this.soc_gain)*(batt_soc-this.soc_trg).^2.*(this.soc_bc_gain.*(~out_range)+out_range);
            elseif strcmp(this.reward_type,"Dense_multicycle")
                fuel = (this.k2*10)*ice_fuel_rate*(Ts)/(veh_dist_tot/(1000));
                eq_fuel = (this.k3*10)*(this.soc_gain)*(batt_soc-this.soc_trg).^2.*(0.3+veh_dist_perc*(1-0.3)).^2/(veh_dist_tot/(1000));
            end
            reward = (this.k1+fuel+eq_fuel).*(~unfeasible)+this.unfeaspenalty.*unfeasible;
        end
    end
end