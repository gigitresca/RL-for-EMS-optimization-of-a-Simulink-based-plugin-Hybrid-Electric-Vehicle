classdef Wheels < matlab.mixin.Copyable
    % object that define the parameters of the four wheels of the vehicle

    properties
        radius
        inertiaR
        inertiaF
    end

    methods
        function this = Wheels(tab_data)
            % Constructor method: Intializes the scalar value of the
            % wheels
            if nargin>0
                this.radius = tab_data.("Radius");
                this.inertiaR = tab_data.("Inertia R");
                this.inertiaF = tab_data.("Inertia F");
            end
        end

        function setWheelsParameters(this,tab_data)
            % Method to set the wheels parameters
            this.radius = tab_data.("Radius");
            this.inertiaR = tab_data.("Inertia R");
            this.inertiaF = tab_data.("Inertia F");
        end
    end
end