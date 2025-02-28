function [env, VehicleEnv] = getEnv(env_type, obs_info, act_info)
    % Vehicle components loading
    if ~exist("DrivingCycleObjects.mat","file")
        generateDrivingCycleObjects("data/training_cycles.xlsx");
        driving_cycles = load("DrivingCycleObjects.mat");
    else
        driving_cycles = load("DrivingCycleObjects.mat");
    end
    
    if ~exist("VehicleComponentsObjects.mat","file")
        [ice,em,gb,ess,fd,wh,veh] = generateVehicleObjects("data/vehicle_data.xlsx");
    else
        load VehicleComponentsObjects
    end
    % Load config file
    config = jsondecode(fileread("src/config/env.json")); % Load JSON
    if ~ismember(config.reward_type, ["Sparse","Dense","Dense_multicycle","ECMS"])
        error('Reward type not valid. Please select between "Sparse","Dense","Dense_multicycle","ECMS"')
    end
    VehicleEnv = VehicleBackKinEnv(obs_info,act_info,veh,wh,fd,gb,em,ess,ice,driving_cycles,config.reward_type);
    % Remove test driving cycles
    if config.multi_cycle
        VehicleEnv.removeTestDrivingCycles;
    else
        if ~ismember(config.cycle_name, fieldnames(driving_cycles))
            error('Driving cycle not available: please select a valid driving cycle')
        end
        VehicleEnv.setSingleCycleTraining(config.cycle_name);

    end
    % Drving cycle resample
    for i = 1:length(VehicleEnv.driving_cycles)
        if VehicleEnv.driving_cycles(i).Ts ~= 1
            VehicleEnv.driving_cycles(i).resample(1);
        end
    end
    % Select the environment
    switch env_type
        case 'Matlab'
            env = VehicleEnv;
        case 'Simulink'
            env = rlSimulinkEnv("VehicleBackKin","VehicleBackKin/Reinforcement Learning controller/RL Agent",obs_info,act_info);
            env.ResetFcn = @(in)VehicleEnv.SimulinkResetFcn(in);
        otherwise
            error('Unsupported env type: %s! Please add the environment file in "src/envs"', env_type);
    end
end