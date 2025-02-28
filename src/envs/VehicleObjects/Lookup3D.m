classdef Lookup3D < matlab.mixin.Copyable
    properties (SetAccess = private)
        brkp1
        brkp2
        brkp3
        tab_data
        var_names
        var_units
    end

    methods
        function this = Lookup3D(brkp1,brkp2,brkp3,tab_data,var_names,var_units)
            % Constructor method: initialize the 2D lookup talbe variables
            if nargin>0
                this.brkp1 = brkp1;
                this.brkp2 = brkp2;
                this.brkp3 = brkp3;
                this.tab_data = tab_data;
                tab_dim = [length(brkp2) length(brkp1)];
                if ~isequal(tab_dim,size(tab_data))
                    error("Axes and table have incompatible sizes")
                end
                this.var_names = var_names;
                this.var_units = var_units;
            end
        end

        function setLookup2DParameters(this,brkp1,brkp2,brkp3,tab_data,var_names,var_units)
            this.brkp1 = brkp1;
            this.brkp2 = brkp2;
            this.brkp3 = brkp3;
            this.tab_data = tab_data;
            tab_dim = [length(brkp2) length(brkp1)];
            if ~isequal(tab_dim,size(tab_data))
                error("Axes and table have incompatible sizes")
            end
            this.var_names = var_names;
            this.var_units = var_units;
        end

        function val = compute(this,x1,x2,x3)
            val = interp3_lim(this.brkp1,this.brkp2,this.brkp3,this.tab_data,x1,x2,x3);
        end

        function resize(this,brkp1_new,brkp2_new,brkp3_new)
            if ~isrow(brkp1_new)
                brkp1_new = brkp1_new';
            end
            if ~iscolumn(brkp2_new)
                brkp2_new = brkp2_new';
            end
            this.tab_data = interp3_lim(this.brkp1,this.brkp2,this.brkp3,this.tab_data,brkp1_new,brkp2_new,brkp3_new);
            this.brkp1 = brkp1_new;
            this.brkp2 = brkp2_new;
            this.brkp3 = brkp3_new;
        end
    end
end