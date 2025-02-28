classdef VehicleBackKinEnv < rl.env.MATLABEnvironment  
    properties
        % Specify and vehicle characteristics   
        veh = Vehicle
        wh = Wheels
        fd = FinalDrive
        gb = Gearbox
        em = ElectricMotor
        ess = Battery
        ice = Engine
        reward = RewardRL
        driving_cycles = DrivingCycle
        normalization = StatesNormalization
        Ts = 1
        sim_input
        sim_output
        states
    end
    
    properties(Access = protected)
        % Initialize internal flag to indicate episode termination
        isdone = false        
    end

    methods              
        function this = VehicleBackKinEnv(obs_info,act_info,veh,wh,fd,gb,em,ess,ice,driving_cycles_vec,reward_type)
            % Contructor method creates an instance of the environment
            % Change class name and constructor name accordingly
            
            if nargin<2
                error('Specify obs_info and act_info')
            end
            % The following line implements built-in functions of RL env
            this = this@rl.env.MATLABEnvironment(obs_info,act_info);
            if nargin>2
                % Set the vehicle paramteres objects and the available
                % driving cycles for the trainig
                updateVehicleParameteresAndDrivingCycles(this,veh,wh,fd,gb,em,ess,ice,driving_cycles_vec,reward_type);
                this.normalization.setNormalizationFromStates(this,1,0)
            end
        end

        function updateVehicleParameteresAndDrivingCycles(this,veh,wh,fd,gb,em,ess,ice,driving_cycles_vec,reward_type)
            % Method to assign the components objects to the relative env
            % properties
            this.veh = veh;
            this.wh = wh;
            this.fd = fd;
            this.gb = gb;
            this.em = em;
            this.ess = ess;
            this.ice = ice;
            % Driving cycles loading
            drv_cyc = fieldnames(driving_cycles_vec);
            for i = 1:length(drv_cyc)
                this.driving_cycles(i,1) = driving_cycles_vec.(drv_cyc{i});
                if this.driving_cycles(i,1).Ts ~= this.Ts
                    this.driving_cycles(i,1).resample(this.Ts);
                end
            end
            % Reward setting
            this.reward = RewardRL(k1=0, ...
                k2=-100, ...
                k3=-147, ...
                soc_bc_gain=1/10, ...
                reward_type=reward_type, ...
                unfeaspenalty=0, ...
                ice_fuel_heating_val=this.ice.fuel_heating_val, ...
                ess_cell_cap=this.ess.cell_cap, ...
                ess_num_cell_series=this.ess.num_cell_series, ...
                ess_volt_nom = this.ess.volt_nom, ...
                ess_soc_trg=this.ess.soc_trg ...
            );
        end

        function states = updateStates(this)
            % Function to update the states of the system given the name of
            % the states
            states_name = this.ObservationInfo.Name;
            states = zeros(length(states_name),1);

            veh_spd = this.sim_input.veh_spd;
            veh_acc = this.sim_input.veh_acc;
            veh_dist_perc = this.sim_input.veh_dist_perc;
            time_perc = this.sim_input.time_perc;
            gb_in_spd = this.sim_output.gb_in_spd(end);
            gb_in_pwr = this.sim_output.gb_in_pwr(end);
            pwr_dmd = this.sim_output.pwr_dmd(end);
            em_spd = this.sim_output.em_spd(end);
            ice_spd = this.sim_output.ice_spd(end);
            ice_state = this.sim_output.ice_state(end);
            batt_soc = this.sim_output.batt_soc(end);
            for i = 1:length(states)
                eval("states(i)"+" = "+states_name(i)+";");
            end
        end
        
        function InitialObservation = reset(this)
            % Method to reset environment to initial state and output initial observation
            
            % Reset output
            this.sim_output = [];            
            % Driving cycle selection
            if length(this.driving_cycles)>1
                cycle_num = randi(length(this.driving_cycles));
                cycle_num(cycle_num==0)=1;
            else
                cycle_num = 1;
            end
            this.sim_input.cycle = this.driving_cycles(cycle_num).name;
            this.sim_input.cycle_num = cycle_num;
            % Initialize the time
            time = 0;
            % Input initialization
            this.setSimulationInput(time);
            this.sim_output.batt_soc = this.ess.soc_init;
            % Power demand with ice off initialization
            [~,fd_spd,gb_in_spd] = this.computeSpeeds(this.sim_input.veh_spd,this.sim_input.veh_gear);
            [~,~,~,F_tot] = this.computeRoadForce(this.sim_input.veh_spd,this.sim_input.veh_acc,this.sim_input.veh_gear,this.sim_input.grade,0);
            [~,fd_pwr,~,gb_in_pwr] = this.computePowerDemand(this.sim_input.veh_gear,F_tot,fd_spd,gb_in_spd);
            this.sim_output.fd_spd = fd_spd;
            this.sim_output.fd_pwr = fd_pwr;
            this.sim_output.gb_in_spd = gb_in_spd;
            this.sim_output.gb_in_pwr = gb_in_pwr;
            this.sim_output.pwr_dmd = gb_in_pwr;
            this.sim_output.ice_state = 0;
            % State initialization
            this.states = zeros(length(this.ObservationInfo.Name),1);
            this.states(strcmp(this.ObservationInfo.Name,"batt_soc")) = this.ess.soc_init;
            InitialObservation = this.normalization.normalize(this.states);

            if strcmp(this.reward.reward_type,"SOC_BC")
                dsoc_init = 0.05;
                dsoc_end = 0.02;
                time_vec = this.driving_cycles(cycle_num).time_simu;
                time_soc_decay_init = 3/4*time_vec(end);
                time_soc_decay_end = time_vec(end)-1;
                soc_decay = (dsoc_init-dsoc_end)/(time_soc_decay_end-time_soc_decay_init);
                soc_bc_fun = @(time) -(soc_decay)*(time-time_soc_decay_init)+(this.ess.soc_trg+dsoc_init);
                soc_bc = (this.ess.soc_trg+dsoc_init)*ones(size(time_vec));
                soc_bc(time_vec>=time_soc_decay_init) = soc_bc_fun(time_vec(time_vec>=time_soc_decay_init));
                this.reward.soc_bc.setLookup1DParameters(time_vec,soc_bc,[],[]);
            end
        end

        function [NextObservation,Reward,isdone] = step(this,action)
            % Method that apply system dynamics and simulates the environment with the 
            % given action for one step.
            % Variables Definition
            time = this.sim_input.time;
            veh_spd = this.sim_input.veh_spd;
            veh_acc = this.sim_input.veh_acc;
            veh_gear = this.sim_input.veh_gear;
            grade = this.sim_input.grade;
            batt_soc = this.sim_output.batt_soc(end);
            fd_spd = this.sim_output.fd_spd(end);
            gb_in_spd = this.sim_output.gb_in_spd(end);
            gb_in_pwr = this.sim_output.gb_in_pwr(end);
            ice_state_0 = this.sim_output.ice_state(end);
            % Engine status definition
            ice_act_perc = action;
            % if time==1030
            %     keyboard
            % end
            if veh_spd==0 || gb_in_pwr<=0
                ice_state = 0;
                ice_act_perc = 0;
                if ice_state_0~=ice_state
                    [~,~,~,F_tot] = this.computeRoadForce(veh_spd,veh_acc,veh_gear,grade,ice_state);
                    [~,fd_pwr,~,gb_in_pwr] = this.computePowerDemand(veh_gear,F_tot,fd_spd,gb_in_spd);
                end
            elseif action>0
                ice_state = 1;
                if ice_state_0~=ice_state
                    [~,~,~,F_tot] = this.computeRoadForce(veh_spd,veh_acc,veh_gear,grade,ice_state);
                    [~,fd_pwr,~,gb_in_pwr] = this.computePowerDemand(veh_gear,F_tot,fd_spd,gb_in_spd);
                end
            elseif action == 0
                ice_state = 0;
                if ice_state_0~=ice_state
                    [~,~,~,F_tot] = this.computeRoadForce(veh_spd,veh_acc,veh_gear,grade,ice_state);
                    [~,fd_pwr,~,gb_in_pwr] = this.computePowerDemand(veh_gear,F_tot,fd_spd,gb_in_spd);
                end
            end
            % Powertrain and battery step with output computation
            [Reward,isdone] = this.powertrainStep(ice_act_perc,fd_spd,gb_in_spd,gb_in_pwr,ice_state,batt_soc);
            % Update simulation input
            time = time+this.Ts;    % Time update --> increase the time variable of 1 time step
            this.setSimulationInput(time);
            % Power demand with ice off computation
            [~,fd_spd,gb_in_spd] = this.computeSpeeds(this.sim_input.veh_spd,this.sim_input.veh_gear);
            [~,~,~,F_tot] = this.computeRoadForce(this.sim_input.veh_spd,this.sim_input.veh_acc,this.sim_input.veh_gear,this.sim_input.grade,ice_state);
            [~,fd_pwr,~,gb_in_pwr] = this.computePowerDemand(this.sim_input.veh_gear,F_tot,fd_spd,gb_in_spd);
            this.sim_output.fd_spd(time/this.Ts+1) = fd_spd;
            this.sim_output.fd_pwr(time/this.Ts+1) = fd_pwr;
            this.sim_output.gb_in_spd(time/this.Ts+1) = gb_in_spd;
            this.sim_output.gb_in_pwr(time/this.Ts+1) = gb_in_pwr;
            % States update and normalization
            this.states = this.updateStates;
            NextObservation = this.normalization.normalize(this.states);           
        end

        function setSingleCycleTraining(this, cycle_name)
            cycle_names = strings(size(this.driving_cycles));
            for i = 1:length(this.driving_cycles)
                cycle_names(i) = this.driving_cycles(i).name;
            end
            idx = cycle_names == cycle_name;
            this.driving_cycles = this.driving_cycles(idx);
        end

        function removeTestDrivingCycles(this)
            test_cycles = [
                "RDE_Poli_wGPS"...
                "RDE_Poli_wGPS_Urban"...
                "RDE_BRA_01"...
                "RDE_BRA_02"...
                "EPA5_Cycle"...
                "Artemis_Motorway"...
                "Artemis_Road"...
                "Artemis_Urban"
            ]';
            cycle_names = strings(size(this.driving_cycles));
            for i = 1:length(this.driving_cycles)
                cycle_names(i) = this.driving_cycles(i).name;
            end
            idx = find(ismember(cycle_names,test_cycles));
            this.driving_cycles(idx) = [];
        end

        function setSimulationInput(this,time)
            this.sim_input.time = time;
            this.sim_input.time_perc = time/this.driving_cycles(this.sim_input.cycle_num).time_simu(end);
            this.sim_input.veh_spd = this.driving_cycles(this.sim_input.cycle_num).veh_spd(time/this.Ts+1);
            this.sim_input.veh_acc = this.driving_cycles(this.sim_input.cycle_num).veh_acc(time/this.Ts+1);
            this.sim_input.veh_dist_perc = this.driving_cycles(this.sim_input.cycle_num).veh_dist(time/this.Ts+1)/this.driving_cycles(this.sim_input.cycle_num).veh_dist(end);
            this.sim_input.veh_gear = this.driving_cycles(this.sim_input.cycle_num).veh_gear(time/this.Ts+1);
            this.sim_input.grade = this.driving_cycles(this.sim_input.cycle_num).grade(time/this.Ts+1);
        end
        
        function [wh_spd,fd_spd,gb_in_spd] = computeSpeeds(this,veh_spd,veh_gear)
            wh_spd = veh_spd/this.wh.radius;
            fd_spd = wh_spd*this.fd.ratio;
            gb_in_spd = fd_spd*this.gb.gear_ratio(veh_gear+1);
        end

        function [veh_eq_mass,F_acc,F_res,F_tot] = computeRoadForce(this,veh_spd,veh_acc,veh_gear,grade,ice_state)
            veh_eq_mass = this.veh.mass+(2*this.wh.inertiaR+2*this.wh.inertiaR)*1/(this.wh.radius^2)+this.ice.inertia*this.gb.gear_ratio...
                (veh_gear+1)^2*this.fd.ratio^2/(this.wh.radius^2)*ice_state+... %Conventional Vehicle Inertia
                this.em.inertia*this.gb.gear_ratio(veh_gear+1)^2*this.fd.ratio^2*this.em.gear_ratio^2/(this.wh.radius^2); %Additional BAS Inertia
            F_acc = veh_eq_mass*veh_acc; %[N]
            F_res = (this.veh.F0+this.veh.mass*9.81*sin(grade*pi/180)+this.veh.F1*(veh_spd*3.6)+this.veh.F2*(veh_spd*3.6)^2); %[N]
            F_tot = F_res+F_acc; %[N]
        end

        function [fd_trq,fd_pwr,gb_eff,gb_in_pwr] = computePowerDemand(this,veh_gear,F_tot,fd_spd,gb_in_spd)
            fd_trq = F_tot*(this.wh.radius/this.fd.ratio); %[Nm] Torque Request at the inlet of the Final Drive
            fd_pwr = fd_trq.*fd_spd; %[W] Power Request at the inlet of the Final Drive
            gb_eff = interp2_lim(this.gb.EffTrqOutMap(veh_gear+1).brkp1,this.gb.EffTrqOutMap(veh_gear+1).brkp2,this.gb.EffTrqOutMap(veh_gear+1).tab_data,fd_spd,abs(fd_trq)); %[-] gb efficiency calculated depending on the output torque 
            gb_in_pwr = fd_pwr.*gb_eff.^sign(-fd_pwr);
        end

        function [Reward,isdone] = powertrainStep(this,action,fd_spd,gb_in_spd,gb_in_pwr,ice_state,batt_soc)
            % Method that apply system dynamics and simulates the environment with the 
            % given action for one step.
            % Variables Definition
            time = this.sim_input.time;
            veh_spd = this.sim_input.veh_spd;
            veh_acc = this.sim_input.veh_acc;
            % Engine and EM speed computation
            ice_spd = max(this.ice.spd_idle,gb_in_spd).*ice_state;
            em_spd = gb_in_spd*this.em.gear_ratio;
            % Physical Limitation    
            em_pwr_max = interp1_lim(this.em.pwr_max.brkp1,this.em.pwr_max.tab_data,em_spd);
            em_pwr_min = interp1_lim(this.em.pwr_min.brkp1,this.em.pwr_min.tab_data,em_spd);
            ice_pwr_max = interp1_lim(this.ice.pwr_max.brkp1,this.ice.pwr_max.tab_data,ice_spd);
            ice_pwr_min = interp1_lim(this.ice.pwr_min.brkp1,this.ice.pwr_min.tab_data,ice_spd);
            ice_trq_max = interp1_lim(this.ice.trq_max.brkp1,this.ice.trq_max.tab_data,ice_spd);
            % POWERTRAIN
            % Engine power definition
            if strcmp(this.ActionInfo.Name,"ice_trq_max_perc")
                ice_trq = action*ice_trq_max;
                ice_pwr = ice_trq*ice_spd;
            elseif strcmp(this.ActionInfo.Name,"ice_pwr_max_perc")
                ice_pwr = action*ice_pwr_max;
            elseif strcmp(this.ActionInfo.Name,"ice_trq")
                ice_trq = action;
                ice_pwr = ice_trq*ice_spd;
            elseif strcmp(this.ActionInfo.Name,"ice_pwr")
                ice_pwr = action;
            else
                error("Unavailable agent action selected")
            end
            unfeasible = 0;
            if ice_pwr>ice_pwr_max
                ice_pwr = ice_pwr_max;
                unfeasible = 1;
            end
            % Gearbox and Electric Machine power demand definition
            if gb_in_pwr>0
                pwr_dmd = gb_in_pwr;
                em_pwr = pwr_dmd - ice_pwr;
                mech_brk_pwr = 0;
            elseif gb_in_pwr<=0 && ice_state==1
                ice_pwr = ice_pwr_min;
                pwr_dmd = gb_in_pwr-ice_pwr;
                em_pwr = interp2_lim(this.em.RegMap.brkp1,this.em.RegMap.brkp2,this.em.RegMap.tab_data,veh_acc,veh_spd).*pwr_dmd; %introduco la mappa sulla rigenerazione in frenata
                mech_brk_pwr = pwr_dmd-em_pwr;
                unfeasible = 1;
            else
                ice_pwr = 0;
                pwr_dmd = gb_in_pwr;
                em_pwr = interp2_lim(this.em.RegMap.brkp1,this.em.RegMap.brkp2,this.em.RegMap.tab_data,veh_acc,veh_spd).*pwr_dmd; %introduco la mappa sulla rigenerazione in frenata
                mech_brk_pwr = pwr_dmd-em_pwr;
            end
            
            if em_pwr>em_pwr_max || em_pwr<em_pwr_min
                em_pwr = em_pwr_max*(em_pwr>em_pwr_max)+em_pwr_min*(em_pwr<em_pwr_min);
                ice_pwr = pwr_dmd-em_pwr;
                unfeasible = 1;
            end
            % Electrical Variables
            em_eff = interp2_lim(this.em.EffPwrMap.brkp1,this.em.EffPwrMap.brkp2,this.em.EffPwrMap.tab_data,em_spd,em_pwr);
            em_pwr_ele = em_pwr.*em_eff.^sign(-em_pwr); %[W] Electrical Power Requested
            batt_pwr = em_pwr_ele+this.ess.accelec; %[W] Unlimited Battery Power Request
            battc_pwr = batt_pwr/this.ess.num_cell; %[W] Battery Cell Power Request
            % Fuel Consumption
            ice_fuel_rate = interp2_lim(this.ice.FuelPwrMap.brkp1,this.ice.FuelPwrMap.brkp2,this.ice.FuelPwrMap.tab_data,ice_spd,ice_pwr); %[kg/s]
            ice_fuel_rate(ice_state==0) = 0;
            ice_fuel_rate(ice_state==1&ice_spd>this.ice.cut_off_lim&ice_pwr<0) = 0;
            ice_co2_rate = ice_fuel_rate*this.ice.co2_mol_mass/this.ice.fuel_mol_mass; %CO2 Emission rate
            % BATTERY
            battc_curr_min = max(0,min(1,(1-10*(batt_soc-this.ess.soc_max))))*this.ess.num_module_parallel*this.ess.cell_curr_max_chg;
            battc_curr_max = max(0,min(1,(1+10*(batt_soc-this.ess.soc_min ))))*this.ess.num_module_parallel*this.ess.cell_curr_max_dis;
            % Resistance
            battc_Rdis = interp1_lim(this.ess.cell_rint_dis.brkp1,this.ess.cell_rint_dis.tab_data,batt_soc);
            battc_Rchg = interp1_lim(this.ess.cell_rint_chg.brkp1,this.ess.cell_rint_chg.tab_data,batt_soc);
            battc_R = battc_Rdis;
            battc_R(battc_pwr<0) = battc_Rchg(battc_pwr<0);
            % Open circuit voltage
            battc_voc = interp1_lim(this.ess.cell_voc.brkp1,this.ess.cell_voc.tab_data,batt_soc);
            % Current
            battc_curr = (battc_voc-(battc_voc.^2-4*battc_R.*battc_pwr).^0.5)./(2*battc_R);
            battc_curr(imag(battc_curr)~=0) = 10^4; %imaginary solutions => P_batt_c exceed max => imposed high current => identified as infeasibly but SOC calculation is possible
            if battc_curr>battc_curr_max || battc_curr<battc_curr_min
                battc_curr = battc_curr_max*(battc_curr>battc_curr_max)+battc_curr_min*(battc_curr<battc_curr_min);
                unfeasible = 1;
            end
            % SoC
            battc_cap = this.ess.cell_cap*3600; %[As]
            batt_soc = batt_soc-battc_curr*this.Ts/battc_cap;
            batt_en_rate = battc_curr*this.Ts/battc_cap*this.ess.cell_cap*this.ess.volt_nom*this.ess.num_cell/1000; %[kWh]
            % batt_co2_rate = -batt_en_rate*this.ess.co2_prod/1000; %[kg]
            % Episode termination
            % is_done = time==this.driving_cycles.(this.sim_input.cycle).time_simu(end)-1 || batt_soc<this.ess.socLow || batt_soc>this.ess.socHigh || unfeasible==1;
            isdone = time==this.driving_cycles(this.sim_input.cycle_num).time_simu(end)-1 || batt_soc<this.ess.soc_low || batt_soc>this.ess.soc_high;
            this.isdone = isdone;
            % Reward collection
            Reward = this.reward.getReward(ice_fuel_rate,batt_soc,batt_pwr,time,this.driving_cycles(this.sim_input.cycle_num).veh_dist(end),this.sim_input.veh_dist_perc,this.sim_input.cycle,this.Ts,isdone,unfeasible);
            % OUTPUT DEFINITION
            this.sim_output.time_simu(time/this.Ts+1) = time;
            % Engine
            this.sim_output.ice_state(time/this.Ts+1) = ice_state;
            this.sim_output.ice_pwr(time/this.Ts+1) = ice_pwr;
            this.sim_output.ice_spd(time/this.Ts+1) = ice_spd;
            this.sim_output.ice_fuel_rate(time/this.Ts+1) = ice_fuel_rate;
            % Electric Machine
            this.sim_output.em_pwr(time/this.Ts+1) = em_pwr;
            this.sim_output.em_pwr_ele(time/this.Ts+1) = em_pwr_ele;
            this.sim_output.em_spd(time/this.Ts+1) = em_spd;
            % Battery
            this.sim_output.batt_pwr(time/this.Ts+1) = batt_pwr;
            this.sim_output.batt_soc(time/this.Ts+2) = batt_soc;
            %Gearbox
            this.sim_output.fd_spd(time/this.Ts+1) = fd_spd;
            this.sim_output.pwr_dmd(time/this.Ts+1) = pwr_dmd;
            this.sim_output.gb_in_pwr(time/this.Ts+1) = gb_in_pwr;
            %Reward
            this.sim_output.reward(time/this.Ts+1) = Reward;
            this.sim_output.unfeasibility(time/this.Ts+1) = unfeasible;
            this.sim_output.action(time/this.Ts+1) = action;
        end

        function [ice_out,ess_out,em_out,gb_out,veh_out,ems_out] = outputPostprocess(this,simulink_out)
            % Method to postprocess the output depending on the environment
            % (if Simulink or Matlab) and have a compliant output
            % expression from both the environments
            if nargin == 2
                % Postprocess Simulink out
                ess_out = simulink_out.SimulationInfo(1).batt_out;
                em_out = simulink_out.SimulationInfo(1).em_out;
                ems_out = simulink_out.SimulationInfo(1).ems_out;
                gb_out = simulink_out.SimulationInfo(1).gb_out;
                ice_out = simulink_out.SimulationInfo(1).ice_out;
                veh_out = simulink_out.SimulationInfo(1).veh_out;
                % Reward and action
                this.sim_output.reward = simulink_out.Reward;
                this.sim_output.action = simulink_out.Action.ice_trq_max_perc;
                this.sim_output.is_done = simulink_out.IsDone;
            else
                % Output struct creation
                ice_out.ice_state = timeseries(this.sim_output.ice_state',this.sim_output.time_simu);
                ice_out.ice_pwr = timeseries(this.sim_output.ice_pwr',this.sim_output.time_simu);
                ice_out.ice_spd = timeseries(this.sim_output.ice_spd',this.sim_output.time_simu);
                ice_out.ice_fuel_rate = timeseries(this.sim_output.ice_fuel_rate',this.sim_output.time_simu);
                ice_out.ice_fuel_cum = timeseries(cumtrapz(this.sim_output.time_simu,this.sim_output.ice_fuel_rate)',this.sim_output.time_simu);
                ice_co2_rate = this.sim_output.ice_fuel_rate*this.ice.co2_mol_mass/this.ice.fuel_mol_mass;
                ice_co2_cum = cumtrapz(this.sim_output.time_simu,ice_co2_rate);
                ice_out.ice_co2_rate = timeseries(ice_co2_rate',this.sim_output.time_simu);
                ice_out.ice_co2_cum = timeseries(ice_co2_cum',this.sim_output.time_simu);
                ess_out.batt_soc = timeseries(this.sim_output.batt_soc(1:end-1)',this.sim_output.time_simu);
                ess_out.batt_pwr = timeseries(this.sim_output.batt_pwr',this.sim_output.time_simu);
                em_out.em_pwr = timeseries(this.sim_output.em_pwr',this.sim_output.time_simu);
                em_out.batt_pwr = timeseries(this.sim_output.em_pwr_ele',this.sim_output.time_simu);
                em_out.em_spd = timeseries(this.sim_output.em_spd',this.sim_output.time_simu);
                gb_out.gb_in_spd = timeseries(this.sim_output.gb_in_spd(1:end-1)',this.sim_output.time_simu);
                gb_out.gb_in_pwr = timeseries(this.sim_output.gb_in_pwr(1:end-1)',this.sim_output.time_simu);
                veh_out.fd_spd = timeseries(this.sim_output.fd_spd(1:end-1)',this.sim_output.time_simu);  
                veh_out.fd_pwr = timeseries(this.sim_output.fd_pwr(1:end-1)',this.sim_output.time_simu);               
                ems_out.ice_pwr_dmd = timeseries(this.sim_output.ice_pwr',this.sim_output.time_simu);
                ems_out.em_pwr_dmd = timeseries(this.sim_output.em_pwr',this.sim_output.time_simu);
                ems_out.pwr_dmd = timeseries(this.sim_output.pwr_dmd',this.sim_output.time_simu);
                this.sim_output.reward = timeseries(this.sim_output.reward',this.sim_output.time_simu);
                this.sim_output.action = timeseries(this.sim_output.action',this.sim_output.time_simu);
                this.sim_output.unfeasibility = timeseries(this.sim_output.unfeasibility',this.sim_output.time_simu);
                % Old output remove
                this.sim_output = rmfield(this.sim_output,"time_simu");
                this.sim_output = rmfield(this.sim_output,"ice_state");
                this.sim_output = rmfield(this.sim_output,"ice_spd");
                this.sim_output = rmfield(this.sim_output,"ice_fuel_rate");
                this.sim_output = rmfield(this.sim_output,"batt_soc");
                this.sim_output = rmfield(this.sim_output,"batt_pwr");
                this.sim_output = rmfield(this.sim_output,"em_pwr_ele");
                this.sim_output = rmfield(this.sim_output,"em_spd");
                this.sim_output = rmfield(this.sim_output,"gb_in_spd");
                this.sim_output = rmfield(this.sim_output,"gb_in_pwr");
                this.sim_output = rmfield(this.sim_output,"fd_spd");
                this.sim_output = rmfield(this.sim_output,"fd_pwr");
                this.sim_output = rmfield(this.sim_output,"ice_pwr");
                this.sim_output = rmfield(this.sim_output,"em_pwr");
                this.sim_output = rmfield(this.sim_output,"pwr_dmd");
            end
            % Torque calculation
            ice_trq = ems_out.ice_pwr_dmd.Data./max(80,(ice_out.ice_spd.Data));    %[Nm]
            ice_out.ice_trq = timeseries(ice_trq,ems_out.ice_pwr_dmd.Time);
            ice_out.ice_trq.Data(isnan(ice_out.ice_trq.Data)) = 0;
            ice_out.ice_bmep = (ice_out.ice_trq*4*pi)/(this.ice.disp)*10;          %[bar]
            em_trq = ems_out.em_pwr_dmd.Data./max(80,(em_out.em_spd.Data));       %[Nm]
            em_out.em_trq = timeseries(em_trq,ems_out.em_pwr_dmd.Time);
            em_out.em_trq.Data(isnan(em_out.em_trq.Data)) = 0;
            % Operating mode calculation
            ems_out.operating_mode = strings(size(ems_out.ice_pwr_dmd.Data));
            ems_out.operating_mode((ice_out.ice_state.Data==1)&(em_out.em_trq.Data>2)&(em_out.em_spd.Data>100))="E-Boost"; 
            ems_out.operating_mode((ice_out.ice_state.Data==1)&(em_out.em_trq.Data<-2)&(em_out.em_spd.Data>100))="Load Point Moving";
            ems_out.operating_mode((ice_out.ice_state.Data==1)&(em_out.em_trq.Data<=2 & em_out.em_trq.Data>=-2)&(em_out.em_spd.Data>100))="ICE Mode"; 
            ems_out.operating_mode((ice_out.ice_state.Data==0)&(em_out.em_trq.Data>2)&(em_out.em_spd.Data>100))="EV Mode";
            ems_out.operating_mode((ice_out.ice_state.Data==0)&(em_out.em_trq.Data<-2)&(em_out.em_spd.Data>100))="Regenerative Braking";
            ems_out.operating_mode((ice_out.ice_state.Data==0)&(em_out.em_trq.Data<=2 & em_out.em_trq.Data>=-2)&(em_out.em_spd.Data>100))="Mechanical Breaking";
            ems_out.operating_mode(em_out.em_spd.Data<=100)="Idle";
            ems_out.operating_mode = categorical(ems_out.operating_mode);
            ems_out.operating_mode = timeseries(ems_out.operating_mode,ems_out.ice_pwr_dmd.Time);
            % Fuel and CO2 economy computation
            ice_out.fuel_economy = ice_out.ice_fuel_cum.Data(end)*1000/(this.driving_cycles.veh_dist(end)/1000);
            ice_out.co2_economy = ice_out.ice_co2_cum.Data(end)*1000/(this.driving_cycles.veh_dist(end)/1000);
            % Positive final drive energy computation
            fd_pwr_pos = veh_out.fd_pwr.Data;
            fd_pwr_pos(fd_pwr_pos<0) = 0;
            veh_out.fd_en = cumtrapz(veh_out.fd_pwr.Time,fd_pwr_pos);
            this.sim_output.ice_out = ice_out;
            this.sim_output.ess_out = ess_out;
            this.sim_output.em_out = em_out;
            this.sim_output.gb_out = gb_out;
            this.sim_output.veh_out = veh_out;
            this.sim_output.ems_out = ems_out;
        end

        function plotTestOutput(this,figure_folder_name)
            % Method to plot the output of a driving cycle test simulation
            run plot_settings.m
            ess_out = this.sim_output.ess_out;
            em_out = this.sim_output.em_out;
            ems_out = this.sim_output.ems_out;
            gb_out = this.sim_output.gb_out;
            ice_out = this.sim_output.ice_out;
            veh_out = this.sim_output.veh_out;
            cycle_name = this.driving_cycles.name;     
            % State of Charge
            h = figure;
            figure_name = 'SOC';
            plot(ess_out.batt_soc.Time,ess_out.batt_soc.Data);
            if strcmp(this.reward.reward_type,"SOC_BC")
                hold on
                plot(this.reward.soc_bc.brkp1,this.reward.soc_bc.tab_data,"--k");
                plot(this.reward.soc_bc.brkp1,this.ess.soc_trg-(this.reward.soc_bc.tab_data-this.ess.soc_trg),"--k","HandleVisibility","off");
                hold off
            end
            grid on
            xlabel("Time [s]")
            ylabel("SoC [-]")
            if strcmp(this.reward.reward_type,"SOC_BC")
                legend('SAC','SOC BC')
            end
            title("State of Charge")
            xlim([0 ess_out.batt_soc.Time(end)])
            % ylim([this.ess.soc_trg-0.1 this.ess.soc_trg+0.1])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Engine power
            h = figure;
            figure_name = 'ICE Power';
            plot(ems_out.ice_pwr_dmd.Time,ems_out.ice_pwr_dmd.Data/1000);
            grid on
            xlabel("Time [s]")
            ylabel("$P_{ICE}$ [kW]")
            title("Engine Power")
            xlim([0 ems_out.ice_pwr_dmd.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Engine speed
            h = figure;
            figure_name = 'ICE speed';
            plot(ice_out.ice_spd.Time,ice_out.ice_spd.Data*30/pi);
            grid on
            xlabel("Time [s]")
            ylabel("$n_{ICE}$ [rpm]")
            title("Engine speed")
            xlim([0 ice_out.ice_spd.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Fuel rate
            h = figure;
            figure_name = 'Fuel rate';
            plot(ice_out.ice_fuel_rate.Time,ice_out.ice_fuel_rate.Data);
            xlabel("Time [s]")
            ylabel("$\dot{m_f}$ [kg/s]")
            grid on
            title("Fuel rate")
            xlim([0 ice_out.ice_fuel_rate.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Fuel consumption
            h = figure;
            figure_name = 'Fuel Consumption';
            plot(ice_out.ice_fuel_cum.Time,ice_out.ice_fuel_cum.Data*1000);
            xlabel("Time [s]")
            ylabel("$m_f$ [g]")
            grid on
            title("Fuel consumption")
            xlim([0 ice_out.ice_fuel_cum.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % CO2 emissions vs SOC plot
            h = figure;
            figure_name = 'CO2 emission vs SOC plot';
            plot(ess_out.batt_soc.Data(end),ice_out.co2_economy,'.','markersize', markersize*5);
            hold on;
            plot([this.ess.soc_trg this.ess.soc_trg],[ice_out.co2_economy*0.8 ice_out.co2_economy*1.2], '--k','linewidth', linewidth,'Color',[0 0 0]+0.05*15);
            xlabel('$SOC_{end}$ [-]'); 
            ylabel('$CO_2$ emissions [g/km]');
            %set(gca,'xlim',[this.ess.soc_trg-0.1 this.ess.soc_trg+0.1],'xtick',(this.ess.soc_trg-0.1):0.05:(this.ess.soc_trg+0.1));
            %set(gca,'ylim',[ice_out.co2_economy*0.8 ice_out.co2_economy*1.2]);
            grid on;
            title("SOC vs Fuel Consumption Comparison");
            legend('SAC','$SOC_{trg}$','Location','northwest')
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % ICE bsfc Map
            h = figure;
            figure_name='BSFC map and operating points';
            this.ice.BsfcMap.resize(linspace(min(this.ice.BsfcMap.brkp1),max(this.ice.BsfcMap.brkp1),1000),linspace(min(this.ice.BsfcMap.brkp2),max(this.ice.BsfcMap.brkp2),1000));
            this.ice.plotBsfcMap('Bmep',Level_BSFC,linewidth,print_fontsize,print_font)
            hold on
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="ICE Mode")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="ICE Mode"),"go","MarkerFaceColor","g","markersize",markersize);
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="E-Boost")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="E-Boost"),"ro","MarkerFaceColor","r","markersize",markersize);
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="Load Point Moving")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="Load Point Moving"),"bo","MarkerFaceColor","b","markersize",markersize);       
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="Idle")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="Idle"),"co","MarkerFaceColor","c","markersize",markersize);   
            hold off
            lgd = legend;
            lgd.String(4:end) = {'ICE Mode' 'E-boost' 'Load Point Moving' 'Idle'};
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Time distribution BSFC map
            [BsfcTD,SpdIdxTD,TimeDistribution] = this.ice.timeDistribution("Bsfc","Bmep",ice_out.ice_bmep.Time,ice_out.ice_bmep.Data,ice_out.ice_spd.Data,ice_out.ice_state.Data);
            h = figure;
            figure_name = 'Bsfc map and time distribution';
            this.ice.plotBsfcMap("Bmep",Level_BSFC,linewidth,print_fontsize,print_font)
            hold on
            plot(10000,0,'o','linewidth',linewidth,'MarkerFaceColor','#0072BD','MarkerEdgeColor','k','markersize',15);
            for i = 1:numel(TimeDistribution)
                if TimeDistribution(i)>0
                    plot(SpdIdxTD(i)*30/pi,BsfcTD(i),'o','linewidth',linewidth,'MarkerFaceColor','#0072BD','MarkerEdgeColor','k','markersize',max(0.0001,TimeDistribution(i)*10),"HandleVisibility","off");
                end
            end
            clear first_point
            hold off
            lgd = legend;
            lgd.String(end) = {'RL'};
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % ICE fuel consumption Map
            h = figure;
            figure_name = 'Fuel map and operating points';
            this.ice.FuelMap.resize(linspace(min(this.ice.FuelMap.brkp1),max(this.ice.FuelMap.brkp1),1000),linspace(min(this.ice.FuelMap.brkp2),max(this.ice.FuelMap.brkp2),1000));
            this.ice.plotFuelMap("Bmep",Level_FC,linewidth,print_fontsize,print_font)
            hold on
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="ICE Mode")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="ICE Mode"),"go","MarkerFaceColor","g","markersize",markersize);
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="E-Boost")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="E-Boost"),"ro","MarkerFaceColor","r","markersize",markersize);
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="Load Point Moving")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="Load Point Moving"),"bo","MarkerFaceColor","b","markersize",markersize);       
            plot(ice_out.ice_spd.Data(ems_out.operating_mode.Data=="Idle")*30/pi,ice_out.ice_bmep.Data(ems_out.operating_mode.Data=="Idle"),"co","MarkerFaceColor","c","markersize",markersize);   
            hold off
            lgd = legend;
            lgd.String(4:end) = {'ICE Mode' 'E-boost' 'Load Point Moving' 'Idle'};
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Time distribution fuel consumption map
            [FcTD,SpdIdxTD,TimeDistribution] = this.ice.timeDistribution("Fuel","Bmep",ice_out.ice_bmep.Time,ice_out.ice_bmep.Data,ice_out.ice_spd.Data,ice_out.ice_state.Data);
            h = figure;
            figure_name = 'Fuel map and time distribution';
            this.ice.plotFuelMap("Bmep",Level_FC,linewidth,print_fontsize,print_font)
            hold on
            plot(10000,0,'o','linewidth',linewidth,'MarkerFaceColor','#0072BD','MarkerEdgeColor','k','markersize',15);
            for i = 1:numel(TimeDistribution)
                if TimeDistribution(i)>0
                    plot(SpdIdxTD(i)*30/pi,FcTD(i),'o','linewidth',linewidth,'MarkerFaceColor','#0072BD','MarkerEdgeColor','k','markersize',max(0.0001,TimeDistribution(i)*10),"HandleVisibility","off");
                end
            end
            clear first_point
            hold off
            lgd = legend;
            lgd.String(end) = {'RL'};
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % EM power
            h = figure;
            figure_name = 'EM Power';
            plot(ems_out.em_pwr_dmd.Time,ems_out.em_pwr_dmd.Data/1000);
            grid on
            xlabel("Time [s]")
            ylabel("$P_{EM}$ [kW]")
            title("Electric Machine Power")
            xlim([0 ems_out.em_pwr_dmd.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % EM speed
            h = figure;
            figure_name = 'EM speed';
            plot(em_out.em_spd.Time,em_out.em_spd.Data*30/pi);
            grid on
            xlabel("Time [s]")
            ylabel("$n_{EM}$ [rpm]")
            title("Electric Machine speed")
            xlim([0 em_out.em_spd.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Gearbox power
            h = figure;
            figure_name = 'Gearbox Power';
            plot(gb_out.gb_in_pwr.Time,gb_out.gb_in_pwr.Data/1000);
            grid on
            xlabel("Time [s]")
            ylabel("$P_{gb}$ [kW]")
            title("Gearbox power")
            xlim([0 gb_out.gb_in_pwr.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Gearbox energy
            h = figure;
            figure_name = 'Gearbox Energy';
            plot(gb_out.gb_in_pwr.Time,cumtrapz(gb_out.gb_in_pwr.Time,gb_out.gb_in_pwr.Data)*2.78e-7);
            grid on
            xlabel("Time [s]")
            ylabel("$E_{gb}$ [kWh]")
            title("Gearbox energy")
            xlim([0 gb_out.gb_in_pwr.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Reward
            h = figure;
            figure_name = 'Reward trend';
            subplot(2,1,1)
            plot(this.sim_output.reward.Time,this.sim_output.reward.Data);
            grid on
            xlabel("")
            ylabel("Reward [-]")
            title("Reward")
            subplot(2,1,2)
            plot(this.sim_output.reward.Time,cumsum(this.sim_output.reward.Data));
            xlabel("Time [s]")
            ylabel("Cumulative reward [-]")
            xlim([0 this.sim_output.reward.Time(end)])
            grid on
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
            % Action
            h = figure;
            figure_name = 'Action';
            plot(this.sim_output.action,"Linewidth",linewidth);
            grid on
            xlabel("Time [s]","Interpreter","latex")
            ylabel("Action [-]","Interpreter","latex")
            title("Agent action","Interpreter","latex")
            xlim([0 this.sim_output.action.Time(end)])
            print_figure(h,print_size,fullfile(figure_folder_name,[figure_name '.' print_format]),print_format,print_font,print_fontsize);
            savefig(h, fullfile(figure_folder_name,[figure_name '.fig']))
        end

        function [x_axis_swp,action_swp,batch_obs_tot,batch_act,states_comb,states_comb_norm] = statesBatchCreation(this,x_axis,sweep_var,n_el,n_act,n_level)
            % Method to create a sweep in the states to analyze the trained
            % policy
            obs_info = this.getObservationInfo;
            obs_info.Name(contains(obs_info.Name,"gb_in_spd")) = "em_spd";
            obs_max = this.normalization.states_max;
            % Settings
            action_swp = linspace(0,1,n_act);
            % States sweep level definition
            eval("x_axis_swp = linspace(0,"+round(obs_max(strcmp(obs_info.Name,x_axis)))+",n_el);");
            if strcmp(x_axis,"veh_spd")
                x_axis_swp = x_axis_swp;
            end
            for i = 1:length(obs_info.Name)
                eval(obs_info.Name(i)+"_level = [" + obs_max(i)*0.2 + " " + obs_max(i)*0.5 + " " + obs_max(i)*0.8 + "];")
            end
            batt_soc_level = [0.3 0.5 0.7];
            % x axis states normalization
            norm_vec = zeros(length(obs_info.Name),n_el);
            norm_vec(strcmp(obs_info.Name,x_axis),:) = x_axis_swp;
            obs_norm = this.normalization.normalize(norm_vec);
            eval("x_axis_swp_norm = obs_norm(strcmp(obs_info.Name," + '"' +x_axis + '"' + "),:);");
            % states sweep normalization
            norm_vec = zeros(length(obs_info.Name),n_level);
            for i = 1:length(obs_info.Name)
                eval("norm_vec(strcmp(obs_info.Name,obs_info.Name(i)),:) = "+ obs_info.Name(i)+"_level;");
            end
            obs_norm = this.normalization.normalize(norm_vec);
            for i = 1:length(obs_info.Name)
                eval(obs_info.Name(i) + "_level_norm = obs_norm(strcmp(obs_info.Name," + '"' + obs_info.Name(i) + '"' + "),:);");
            end
            % States combination table-speed sweep
            eval("states_comb = combvec("+sweep_var(1)+"_level,"+sweep_var(2)+"_level,"+sweep_var(3)+"_level);")
            states_comb = table(states_comb(1,:)',states_comb(2,:)',states_comb(3,:)',repmat(em_spd_level',[length(states_comb)/n_level 1]),'VariableNames',[sweep_var,"em_spd"]);
            eval("states_comb_norm = combvec("+sweep_var(1)+"_level_norm,"+sweep_var(2)+"_level_norm,"+sweep_var(3)+"_level_norm);")
            states_comb_norm = table(states_comb_norm(1,:)',states_comb_norm(2,:)',states_comb_norm(3,:)',repmat(em_spd_level_norm',[length(states_comb_norm)/length(veh_spd_level_norm) 1]),'VariableNames',[sweep_var,"em_spd"]);
            % Batch creation
            batch_obs_tot = zeros(size(obs_norm,1),n_el,length(states_comb_norm.Variables));
            eval("batch_obs_tot(strcmp(obs_info.Name,"+'"'+x_axis+'"'+"),:,:) = repmat(x_axis_swp_norm,length(states_comb_norm.Variables),1)';");
            for k = 1:length(states_comb_norm.Properties.VariableNames)
                eval("batch_obs_tot(strcmp(obs_info.Name,states_comb_norm.Properties.VariableNames{k}),:,:) = repmat(states_comb_norm."+states_comb_norm.Properties.VariableNames{k}+",1,n_el)';");
            end
            batch_act = repmat(action_swp,[n_el,1]);
            batch_act = batch_act(:)'; 
        end
    end

    methods (Static)
        function in = SimulinkResetFcn(in)  
            % Function input
            % driving_cycles_vec = evalin("base","driving_cycles");
            VehicleEnv = evalin("base","VehicleEnv");
            reward = VehicleEnv.reward;
            ess = VehicleEnv.ess;
            driving_cycles = VehicleEnv.driving_cycles;
            % Driving cycle selection
            if length(driving_cycles)>1
                cycle_num = randi(length(driving_cycles));
                cycle_num(cycle_num==0)=1;
            else
                cycle_num = 1;
            end
            cycle = driving_cycles(cycle_num);
            % Vehicle speed
            cycle_insimulink.veh_spd = timeseries(cycle.veh_spd,cycle.time_simu,"Name","veh_spd");
            % Vehicle gear
            cycle_insimulink.veh_gear = timeseries(cycle.veh_gear,cycle.time_simu,"Name","veh_gear");
            % Vehicle distance
            cycle_insimulink.veh_dist = timeseries(cycle.veh_dist,cycle.time_simu,"Name","veh_dist");
            % Grade
            cycle_insimulink.grade = timeseries(cycle.grade,cycle.time_simu,"Name","grade");
            
            assignin("base","cycle",cycle);
            assignin("base","cycle_insimulink", cycle_insimulink);
            assignin("base","cycleTime",cycle.time_simu(end))
        
            % SOC_BC reward constraints creation
            if strcmp(reward.reward_type,"SOC_BC")
                dsoc_init = 0.075;
                dsoc_end = 0.01;
                time_vec = driving_cycles{cycle_num}.time_simu;
                time_soc_decay_init = 3/4*time_vec(end);
                time_soc_decay_end = time_vec(end)-1;
                soc_decay = (dsoc_init-dsoc_end)/(time_soc_decay_end-time_soc_decay_init);
                soc_bc_fun = @(time) -(soc_decay)*(time-time_soc_decay_init)+(ess.soc_trg+dsoc_init);
                soc_bc = (ess.soc_trg+dsoc_init)*ones(size(time_vec));
                soc_bc(time_vec>=time_soc_decay_init) = soc_bc_fun(time_vec(time_vec>=time_soc_decay_init));
                reward.soc_bc.setLookup1DParameters(time_vec,soc_bc,[],[]);
                reward_bc = timeseries(reward.soc_bc.tab_data,reward.soc_bc.brkp1,"Name","soc_bc");
                assignin("base","reward_bc", reward_bc);
            end
            % Assign vehicle components
            assignin("base","reward", reward);
            assignin("base","veh", VehicleEnv.veh);
            assignin("base","wh", VehicleEnv.wh);
            assignin("base","fd", VehicleEnv.fd);
            assignin("base","gb", VehicleEnv.gb);
            assignin("base","em", VehicleEnv.em);
            assignin("base","ess", VehicleEnv.ess);
            assignin("base","ice", VehicleEnv.ice);

            % States selection
            blk = "VehicleBackKin/Reinforcement Learning controller/Bus states";
            obs_info = VehicleEnv.getObservationInfo;
            bus_in = get_param(blk,"InputSignalNames");
            selector_in = zeros(size(obs_info.Name));
            for i = 1:length(obs_info.Name)
                [~,idx] = find(strcmp(bus_in,"<"+obs_info.Name(i)+">"));
                if ~isempty(idx)
                    selector_in(i) = idx;
                end
            end
            assignin("base","selector_in", selector_in);
            assignin("base","transfer_fun_coeff", 7*VehicleEnv.Ts/10);
        end
    end
end
