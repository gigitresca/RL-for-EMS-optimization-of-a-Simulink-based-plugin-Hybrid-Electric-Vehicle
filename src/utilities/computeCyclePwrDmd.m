function computeCyclePwrDmd(par)
    energy_ratio=zeros(1,numel(par.driving_cycle));
    soc_end_est=zeros(1,numel(par.driving_cycle));
    soc_end_round=zeros(1,numel(par.driving_cycle));
    unfeas = zeros(1,numel(par.driving_cycle));
    time_unfeas = cell(1,numel(par.driving_cycle));
    for i = 1:numel(par.driving_cycle)
        time_simu = par.driving_cycle(i).time_simu;
        wh_spd = par.driving_cycle(i).veh_spd/par.wh.radius;
        fd_spd = wh_spd*par.fd.ratio;
        gb_in_spd = fd_spd.*par.gb.gear_ratio(par.driving_cycle(i).veh_gear+1);
        ice_spd = max(par.ice.spd_idle,gb_in_spd);
        em_spd = gb_in_spd*par.em.gear_ratio;
        veh_eq_mass = par.veh.mass+(2*par.wh.inertiaR+2*par.wh.inertiaR)*1/(par.wh.radius^2)+par.ice.inertia*par.gb.gear_ratio...
                (par.driving_cycle(i).veh_gear+1).^2*par.fd.ratio^2/(par.wh.radius^2)+... %Conventional Vehicle Inertia
                par.em.inertia*par.gb.gear_ratio(par.driving_cycle(i).veh_gear+1).^2*par.fd.ratio^2*par.em.gear_ratio^2/(par.wh.radius^2); %Additional BAS Inertia
        F_acc = veh_eq_mass.*par.driving_cycle(i).veh_acc; %[N]
        F_res = (par.veh.F0+par.veh.mass*9.81*sin(par.driving_cycle(i).grade*pi/180)+par.veh.F1*(par.driving_cycle(i).veh_spd*3.6)+par.veh.F2*(par.driving_cycle(i).veh_spd*3.6).^2); %[N]
        F_tot = F_res+F_acc; %[N]
        fd_trq = F_tot*(par.wh.radius/par.fd.ratio); %[Nm] Torque Request at the inlet of the Final Drive
        fd_pwr = fd_trq.*fd_spd; %[W] Power Request at the inlet of the Final Drive
        gb_eff=zeros(length(gb_in_spd),1);
        for j=1:length(gb_in_spd)
            gb_eff(j,1) = interp2_lim(par.gb.EffTrqOutMap(par.driving_cycle(i).veh_gear(j)+1).brkp1,par.gb.EffTrqOutMap(par.driving_cycle(i).veh_gear(j)+1).brkp2,par.gb.EffTrqOutMap(par.driving_cycle(i).veh_gear(j)+1).tab_data,gb_in_spd(j),abs(fd_trq(j))); %[-] gb efficiency calculated depending on the output torque 
        end
        gb_in_pwr = fd_pwr.*gb_eff.^sign(-fd_pwr); %[W]
        em_pwr_max = interp1_lim(par.em.pwr_max.brkp1,par.em.pwr_max.tab_data,em_spd);
        ice_pwr_max = interp1_lim(par.ice.pwr_max.brkp1,par.ice.pwr_max.tab_data,ice_spd);
        pwr_unfeas = gb_in_pwr>(em_pwr_max+ice_pwr_max);
        em_pwr=gb_in_pwr;
        em_pwr(gb_in_pwr<0) = 0; %Neglect the regenerative braking energy recover beacuse there are some assumption in the model that conduct to a smaller energy recquired in the battery
        em_eff = interp2_lim(par.em.EffPwrMap.brkp1,par.em.EffPwrMap.brkp2,par.em.EffPwrMap.tab_data,em_spd,em_pwr);
        em_pwr_ele = em_pwr.*em_eff.^sign(-em_pwr); %[W] Electrical Power Requested
        batt_pwr = em_pwr_ele+par.ess.accelec; %[W] Unlimited Battery Power Request
        batt_en = cumtrapz(par.driving_cycle(i).time_simu,batt_pwr); %[J]
        batt_en_nom = par.ess.volt_nom*par.ess.cell_cap.*par.ess.num_cell_series*3.6e+3; %[J]
        energy_ratio(i)=batt_en(end)/batt_en_nom;
        if ~isempty(time_simu(pwr_unfeas==1))
            unfeas(i) = 1;
            time_unfeas{i,1} = time_simu(pwr_unfeas==1);
            h(1) = figure(1);
            figure_name{1} = 'Gearbox power with limits';
            plot(time_simu,gb_in_pwr/1000,"LineWidth",2)
            hold on
            plot(time_simu,(ice_pwr_max+em_pwr_max)/1000,"LineWidth",2);
            plot(time_simu(pwr_unfeas==1),gb_in_pwr(pwr_unfeas==1)/1000,'o',"HandleVisibility","off","MarkerEdgeColor","k","MarkerFaceColor","#A2142F","Markersize",8)
            hold off
            grid on
            xlabel("Time [s]")
            ylabel("Power [kW]")
            legend("Power demand","Powertrain limit")
            % title(par.driving_cycle(i).name+" power unfeasibility")
        elseif energy_ratio(i)<1
            soc_end_est(i)=par.ess.soc_init-energy_ratio(i);
            if soc_end_est(i)>par.ess.soc_trg
                soc_end_round(i)=ceil(soc_end_est(i)*10)/10-0.05*((ceil(soc_end_est(i)*10)/10-soc_end_est(i))>0.05);
            else
                soc_end_round(i)=par.ess.soc_trg;
            end
        else
            soc_end_round(i)=par.ess.soc_trg;
        end
    end
    disp("Unfeasible cycle: ")
    if nnz(unfeas==1)>0
        unfeas_cycle=par.driving_cycle(soc_end_round~=par.ess.soc_trg).name;
        disp(unfeas_cycle+": powertrain power limit at "+num2str(time_unfeas{unfeas==1})+" s")
        return
    elseif nnz(soc_end_round~=par.ess.soc_trg)>0
        unfeas_cycle=par.driving_cycle(soc_end_round~=par.ess.soc_trg).name;
        disp(unfeas_cycle+": soc target unreachable")
    else
        disp('NONE')
    end
end

