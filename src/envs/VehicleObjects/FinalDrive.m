classdef FinalDrive < matlab.mixin.Copyable
    % object that define the parameters of the four wheels of the vehicle

    properties
        ratio
        eff_max
        inertia
    end

    methods
        function this = FinalDrive(tab_data)
            % Constructor method: Intializes the scalar value of the
            % final drive
            if nargin>0
                this.ratio = tab_data.("Ratio");
                this.eff_max = tab_data.("Efficiency max");
                this.inertia = tab_data.("Inertia");
            end
        end

        function setFinalDriveParameters(this,tab_data)
            % Method to set the final drive parameters
            this.ratio = tab_data.("Ratio");
            this.eff_max = tab_data.("Efficiency max");
            this.inertia = tab_data.("Inertia");
        end
    end
end