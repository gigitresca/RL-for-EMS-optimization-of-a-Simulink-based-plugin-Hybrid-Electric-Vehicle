function [ice,em,gb,ess,fd,wh,veh] = generateVehicleObjects(xlsx_file)
    % Function to generate the objects of the main components of the modelled
    % vehicle
    %% PLOT PARAMETERES
     % Options
    printFormat = 'png';
    printFont = 'Arial';
    printFontSize = 20;
    printSize = [12 8];
    % Plot Parameters
    Fontsize = 20;
    LineWidth = 2;
    FontName=printFont;
    % Level Fuel Consumption contour plot [kg/h]
    Level_FC{1} = 0:1:5;
    Level_FC{2} = 7;
    Level_FC{3} = 10:5:30;
    % Level BSFC contour plot [g/kWh]
    Level_BSFC{1} = 200:5:220;
    Level_BSFC{2} = 220:10:240;
    Level_BSFC{3} = 240:20:300;
    Level_BSFC{4} = 300:50:400;
    % Level EM efficiency contour plot [-]
    Level_Eff{1} = 0.7:0.05:0.9;
    Level_Eff{2} = 0.9:0.01:0.95;
    % Level Regeneration contour plot [-]
    Level_Reg{1} = 0:0.1:0.8;
    Level_Reg{2} = 0.8:0.05:0.95;
    
    %% ENGINE OBJECT CREATION
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Engine","VariableNamingRule","preserve");
    ice = Engine(tab_data);
    % Limits
    tab_data = readcell(xlsx_file,"Sheet","Engine-Limits");
    var_info = cell2table(tab_data(2,:),"VariableNames",tab_data(1,:),"RowNames","Unit");
    tab_data = cell2table(tab_data(3:end,:),"VariableNames",tab_data(1,:));
    ice.setLimits(tab_data,var_info);
    % Fuel torque map
    tab_data = readcell(xlsx_file,"Sheet","Engine-FuelTrqMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    ice.setMap(tab_data,var_info);
    % Fuel power map
    tab_data = readcell(xlsx_file,"Sheet","Engine-FuelPwrMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    ice.setMap(tab_data,var_info);
    % Bsfc torque map
    tab_data = readcell(xlsx_file,"Sheet","Engine-BsfcTrqMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    ice.setMap(tab_data,var_info);
    % Bsfc power map
    tab_data = readcell(xlsx_file,"Sheet","Engine-BsfcPwrMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    ice.setMap(tab_data,var_info);
    % Efficiency torque map
    tab_data = readcell(xlsx_file,"Sheet","Engine-EffTrqMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    ice.setMap(tab_data,var_info);
    % Optimal Operating Line (OOL) computation
    ice.minBsfcComputation("Torque");
    ice.minBsfcComputation("Power");
    % Maps plot
    figure
    ice.plotBsfcMap("Bmep",Level_BSFC,LineWidth,Fontsize,FontName)
    figure
    ice.plotFuelMap("Bmep",Level_FC,LineWidth,Fontsize,FontName)
    
    %% ELECTRIC MOTOR OBJECT CREATION
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Electric motor","VariableNamingRule","preserve");
    em = ElectricMotor(tab_data);
    % Limits
    tab_data = readcell(xlsx_file,"Sheet","Electric motor-Limits");
    var_info = cell2table(tab_data(2,:),"VariableNames",tab_data(1,:),"RowNames","Unit");
    tab_data = cell2table(tab_data(3:end,:),"VariableNames",tab_data(1,:));
    em.setLimits(tab_data,var_info);
    % Efficiency power map
    tab_data = readcell(xlsx_file,"Sheet","Electric motor-EffPwrMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    em.setMap(tab_data,var_info);
    % Regeneration map
    tab_data = readcell(xlsx_file,"Sheet","Electric motor-RegMap");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    em.setMap(tab_data,var_info);
    em.maxEffComputation;
    % Maps plot
    figure
    em.plotEffMap(Level_Eff,LineWidth,Fontsize,FontName)
    figure
    em.plotRegMap(Level_Reg,LineWidth,Fontsize,FontName)
    
    %% GEARBOX OBJECT CREATION
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Gearbox","VariableNamingRule","preserve");
    gb = Gearbox(tab_data);
    % Efficiency maps intialization and 1st efficiency map
    tab_data = readcell(xlsx_file,"Sheet","Gearbox-EffMap1");
    tab_data{3,2} = NaN;
    var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
    var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
    tab_data = cell2mat(tab_data(3:end,2:end));
    gb.initializeEfficiencyMap(tab_data,var_info);
    gb.setEfficiencyMap(tab_data,var_info,1);
    % Efficiency maps from 2nd to last gear
    for i = 3:length(gb.gear_idx)
        sheet_name = "Gearbox-EffMap"+string(gb.gear_idx(i));
        tab_data = readcell(xlsx_file,"Sheet",sheet_name);
        tab_data{3,2} = NaN;
        var_info = [tab_data(3:4,1),tab_data(1:2,2),tab_data(1:2,1)];
        var_info = cell2table(var_info(2,:),"VariableNames",var_info(1,:),"RowNames","Unit");
        tab_data = cell2mat(tab_data(3:end,2:end));
        gb.setEfficiencyMap(tab_data,var_info,gb.gear_idx(i));
    end
    clear sheet_name
    
    %% BATTERY OBJECT CREATION
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Battery","VariableNamingRule","preserve");
    ess = Battery(tab_data);
    % Battery usage strategy
    ess.setBatteryStrategy(0.5,0.5);
    % Cell data
    tab_data = readcell(xlsx_file,"Sheet","Battery-VoC-Rint");
    var_info = cell2table(tab_data(2,:),"VariableNames",tab_data(1,:),"RowNames","Unit");
    tab_data = cell2table(tab_data(3:end,:),"VariableNames",tab_data(1,:));
    ess.setCellData(tab_data,var_info);
    figure
    ess.plotCellData(LineWidth,Fontsize,FontName)
    
    %% FINAL DRIVE OBJECT
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Final drive","VariableNamingRule","preserve");
    fd = FinalDrive(tab_data);
    
    %% WHEELS OBJECT
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Wheels","VariableNamingRule","preserve");
    wh = Wheels(tab_data);
    
    %% VEHICLE OBJECT
    % Main parameters
    tab_data = readtable(xlsx_file,"Sheet","Vehicle","VariableNamingRule","preserve");
    veh = Vehicle(tab_data);
    % Plot road resistance
    figure
    veh.plotRoadReistance(LineWidth,Fontsize,FontName)

    %% SAVING COMPONENTS OBJECTS
    save("data\VehicleComponentsObjects","ice","em","gb","ess","fd","wh","veh")
end