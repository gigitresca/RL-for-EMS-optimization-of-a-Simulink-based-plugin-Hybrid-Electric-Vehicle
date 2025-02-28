classdef Lookup1D < matlab.mixin.Copyable
    properties (SetAccess = private)
        brkp1
        tab_data
        var_names
        var_units
    end

    methods
        function this = Lookup1D(brkp1,tab_data,var_names,var_units)
            % Constructor method: initialize the 1D lookup table variables
            if nargin>0
                this.brkp1 = brkp1;
                this.tab_data = tab_data;
                if size(brkp1)~=size(tab_data)
                    error("Axis and table have incompatible sizes")
                end
                this.var_names = var_names;
                this.var_units = var_units;
            end
        end

        function setLookup1DParameters(this,brkp1,tab_data,var_names,var_units)
            % Method to set the 1D lookup table variables
            this.brkp1 = brkp1;
            this.tab_data = tab_data;
            if size(brkp1)~=size(tab_data)
                error("Axis and table have incompatible sizes")
            end
            this.var_names = var_names;
            this.var_units = var_units;
        end

        function val = compute(this,x1)
            val = interp1_lim(this.brkp1,this.tab_data,x1);
        end

        function resize(this,brkp1_new)
            this.tab_data = interp1_lim(this.brkp1,this.tab_data,brkp1_new);
            this.brkp1 = brkp1_new;
        end
    end
end