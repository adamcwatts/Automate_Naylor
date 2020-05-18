classdef Naylor_Replicate_Basic < dynamicprops
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Time
        HeatFlux
        Temp 
        RH
        file_name = "file name"
        option_smooth {logical}
        steady_state_index {int32}
    end
    
    methods
        function obj = Naylor_Replicate_Basic(data, file_name, ss_index)
            %Construct an instance of this class
         obj.Time = data(:,1);
         obj.HeatFlux = data(:,2);
         obj.Temp = data(:,3);
         obj.RH = data(:,4);
         obj.file_name = file_name;
         obj.steady_state_index = ss_index;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

