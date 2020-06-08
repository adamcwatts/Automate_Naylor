classdef ISO_Replicate < handle
    %ISO_REPLICATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Time
        Temp
        Temp_Baselined
        Air_Temp
        RH
        delta_Temp_raw
        delta_Temp_smooth
        file_name
        Phi
        Phi_BaseLined
        Phi_step_idx
        Phi_Integrated_Area
        Phi_Integrated_Reference
        Energy_Joules
 
    
    end
    
    methods(Static)
        
        function Phi_step_idx = baseliner(Data, starting_idx, step_change_threshold)
            
            if isempty(step_change_threshold)
                step_change_threshold = 5e-4;  % THRESHOLD DIFF
            end 

            [~, max_idx] = max(Data);
            reversed_data = flip(Data(1:max_idx));

            step_idx = 0;
            for j=1:length(reversed_data)

                if abs(reversed_data(j)) < step_change_threshold
                   step_idx = (j - 1) ;
                   break
                end
            end

            Phi_step_idx = max_idx - step_idx + (starting_idx);  % when starting idx is 1 and greater. (corrects idx for double peaks)
        end 
        
        
    end
    
    methods
        function obj = ISO_Replicate(file_path)
            % Assumes file_path is excel file with 1 row allocated for the
            % header titles
            all_data_table = readtable(file_path);
            obj.Time = all_data_table{:,1};
            obj.Temp = all_data_table{:,2};
            obj.Air_Temp = all_data_table{:,3};
            obj.RH = all_data_table{:,4};
            obj.file_name = file_path;
            
        end
        
        function  smooth_temp_diff(obj, smoothing_window)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if isempty(smoothing_window)
                smoothing_window = 20;
            end 
            
            obj.delta_Temp_raw = diff(obj.Temp);
            obj.delta_Temp_smooth = smoothdata(obj.delta_Temp_raw,"gaussian",smoothing_window);
            
        end
        
        function integrate_flux(obj)
        Phi_Base = obj.Phi_BaseLined;
        time = obj.Time;
        step_index = obj.Phi_step_idx;
            
            
        area_to_negate = 1/2 * Phi_Base(end) * (time(end)-time(step_index)) ;
        total_area = trapz(time(step_index + 1:end), Phi_Base(step_index:end)) ;  % j /m^2
        adjusted_area = total_area - area_to_negate; % j / m^2

 
        slope = (Phi_Base(end) - Phi_Base(step_index)) /(time(end) - time(step_index + 1));
        baseline_function = @(t) Phi_Base(end) + slope.*(t - t(end));
        Reference_results = baseline_function(time);
        
        obj.Phi_Integrated_Area = adjusted_area; % j / m^2
        obj.Phi_Integrated_Reference = Reference_results; % j / m^2 baseline vector used for integration
       end
        
        
        
        
    end
end

