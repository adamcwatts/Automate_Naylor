classdef Naylor_Replicate
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Time
        HeatFlux
        Temp 
        RH
        file_name = "file name"
        starting_time_minutes = 6  % look for step changes begining at 6 min in
        starting_index
        steady_state_time_minutes = 10  % min
        steady_state_index {int32}
        mean_delta_time_seconds
        option_smooth_import_data = false 
%         RH_step_idx_Local {int32}
%         RH_step_time_Local

        RH_step_and_steady_state_idx {int32}
        RH_step_and_steady_state_time
        
        RH_step_idx_Global {int32}
        RH_step_time_Global
        
    end 
   
    
    methods
         function obj = Naylor_Replicate(file_path_location)
             
         % Construct an instance of this class
         % Any code not using output argument (obj)
         
         temp_data = readmatrix(file_path_location,"Range",'A17:X1000');
         temp_data = temp_data(:,[1, 8, 15, 16]); % select 1st column, time; 2nd column, flux; 3rd column, Chamber temp; 4th column, amb RH
         data = rmmissing(temp_data,2);
         obj.Time = data(:,1);
         obj.HeatFlux = data(:,2);
         obj.Temp = data(:,3);
         obj.RH = data(:,4);
         obj.mean_delta_time_seconds = mean(diff(obj.Time));
         
         obj.starting_index =  ceil((obj.starting_time_minutes * 60) /  obj.mean_delta_time_seconds);
         obj.steady_state_index = ceil((obj.steady_state_time_minutes * 60 ) / obj.mean_delta_time_seconds);
         name = split(file_path_location,"/"); obj.file_name = name{2,1};
         end 
        
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

