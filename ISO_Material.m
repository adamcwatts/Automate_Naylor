classdef ISO_Material < handle  % objects by reference, handle class 
    %ISO_MATERIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Material_Type
        Square_Length_mm
        Area_m  % meter squared
        thickness_mm
        thickness_m
        Weight_Density_g_m_sq
        rho_fabric_g
        rho_fabric_kg
        cp_fabric
        mass_g
        mass_kg
        R_cf
        file_count
        Replicate
        Plot_title
        Compiled_Integration_Results
        colors = {[0, 0.4470, 0.7410],... 
        [0.8500, 0.3250, 0.0980],...
        [0.9290, 0.6940, 0.1250],...
        [0.4660, 0.6740, 0.1880],...
        [0.4940, 0.1840, 0.5560],...
        [0.6350, 0.0780, 0.1840],...
        [0.3010, 0.7450, 0.9330],...
        };
        
    end
    
    methods(Static)
        function export_results()
            % call whos function in the base workspace 
            s = evalin('base', 'whos'); 
            
            % find the objects with the type 'ISO_Material' from the workspace:
            matches= strcmp({s.class}, 'ISO_Material');
            my_class_variables = {s(matches).name};
            n_classes = length(my_class_variables);
           
            n_max_replicates = 0;

            % Find maximum file count for each class instance
            for var = my_class_variables
                ob =  evalin('base', var{1});

                if ob.file_count > n_max_replicates
                    n_max_replicates = ob.file_count;
                    
                end 
            end 

            master_data =  nan(n_max_replicates, n_classes);
            
            for k = 1:size(my_class_variables, 2)
                ob =  evalin('base', my_class_variables{k});
                master_data(1:ob.file_count, k) = ob.Compiled_Integration_Results{:,1};  
            end 
            
            T = array2table(master_data,'VariableNames', my_class_variables);
            fn ="Integrated_Flux_" + "ISO_METHOD" + ".csv";
            writetable(T,fn);
            disp("Exported: " + fn)
            
        end 
        
    end 
    
    methods
        function obj = ISO_Material(file_path_locations_of_reps)
            FPLR = file_path_locations_of_reps;
            obj.file_count = length(FPLR); % File Count
            
           % store replicate objects into structure array 
            for i=1:obj.file_count
                obj.Replicate(i).Data = ISO_Replicate(FPLR{i});
                obj.Replicate(i).Data.smooth_temp_diff([])
            end 
            
            
        end
        
        function Material_Data(obj, Type ,Square_Length_mm, thickness_mm, Weight_Density, cp_fabric, R_cf)
            obj.Material_Type = Type;
            obj.Area_m = Square_Length_mm^2 * (1/1000)^2;  % mm^2 to m^2
            obj.thickness_mm = 0.44; % mm
            obj.thickness_m = thickness_mm / 1000; % m
            obj.Weight_Density_g_m_sq = Weight_Density; % g/m^2
            obj.rho_fabric_g = Weight_Density / obj.thickness_m;  % g/m^3
            obj.rho_fabric_kg = obj.rho_fabric_g / 1000;  % kg/m^3
            obj.cp_fabric = cp_fabric; % J / kg C
            obj.mass_g = obj.Weight_Density_g_m_sq  * obj.Area_m;  % g
            obj.mass_kg = obj.mass_g / 1000;   % kg   
            obj.R_cf = R_cf; 
        end 
        
        function Plot_Smoothing_Delta_T(obj)
        
            for i=1:obj.file_count
                figure(i); clf;  
                hold on
                plot(obj.Replicate(i).Data.Time(2:end,1), obj.Replicate(i).Data.delta_Temp_raw(:,:), "Color", '#787878')
                plot(obj.Replicate(i).Data.Time(2:end,1), obj.Replicate(i).Data.delta_Temp_smooth(:,:),"LineWidth",2, "Color", obj.colors{i})
                xlabel('time (s)','Interpreter',"latex")
                ylabel('$\frac{dT}{dt} \ \ [\frac{^{\circ}C}{s}]$','Interpreter', "latex")
                legend({'Raw Delta', "Gaussian Smoothed"})
                title({'Differential Smoothing' ,obj.Replicate(i).Data.file_name}, 'Interpreter', 'none')
            end 
            
            
        end 
        
        function Plot_Raw_Data(obj)
           
            figure(1); clf;
            subplot(3,1,1);
   
            hold on
            for i =1:obj.file_count
                plot(obj.Replicate(i).Data.Time, ...
                obj.Replicate(i).Data.Temp, ...
                "LineWidth",1.75,"DisplayName",....
                obj.Replicate(i).Data.file_name,...
                "Color", obj.colors{i})
        
            end 
            xlabel('time (s)','Interpreter',"latex")
            ylabel(['Temperature [' char(176) 'C]'])
            legend('Interpreter', 'none')
            hold off
            
            subplot(3,1,2);
            hold on
            for i =1:obj.file_count
                plot(obj.Replicate(i).Data.Time(2:end), ...
                obj.Replicate(i).Data.delta_Temp_smooth, ...
                "LineWidth",1.75,"DisplayName",....
                obj.Replicate(i).Data.file_name,...
                "Color", obj.colors{i})
        
            end 
            xlabel('time (s)','Interpreter',"latex")
            ylabel('$\frac{dT}{dt} \ \ [\frac{^{\circ}C}{s}]$','Interpreter', "latex")
            legend('Interpreter', 'none')
            hold off
            
            
            subplot(3,1,3);
            hold on
            for i =1:obj.file_count
                plot(obj.Replicate(i).Data.Time, ...
                obj.Replicate(i).Data.RH, ...
                "LineWidth",1.75,"DisplayName",....
                obj.Replicate(i).Data.file_name,...
                "Color", obj.colors{i})
        
            end 
            xlabel('time (s)','Interpreter',"latex")
            ylabel('RH [%]')
            legend('Interpreter', 'none', 'Location', "SouthEast")
            hold off
            
            
            
        end 
        
        function Calculate_Flux(obj)
            
            for i=1:obj.file_count
                phi = obj.Replicate(i).Data.delta_Temp_smooth * (obj.rho_fabric_kg * obj.thickness_m * obj.cp_fabric) + ...
                    (obj.Replicate(i).Data.Temp(2:end) - obj.Replicate(i).Data.Air_Temp(2:end)) / obj.R_cf;
                
                obj.Replicate(i).Data.Phi = phi;
            end 
        end 
        
        function Plot_Flux(obj)
            Min_RH = zeros(obj.file_count, 1);
            Max_RH = zeros(obj.file_count, 1);
            
            for i=1:obj.file_count
                Min_RH(i,1) = min(obj.Replicate(i).Data.RH);
                Max_RH(i,1) = max(obj.Replicate(i).Data.RH);
            end 
            
            Min_RH_Mean = mean(Min_RH);
            Max_RH_Mean = mean(Max_RH);
            
            str1 = sprintf('%d',round(Min_RH_Mean));
            str2 = sprintf('%d',round(Max_RH_Mean));
         
            figure(1);clf;
            hold on
            for i=1:obj.file_count
               plot(obj.Replicate(i).Data.Time(2:end),obj.Replicate(i).Data.Phi, "LineWidth", 2, "Color", obj.colors{i})
            end 
            
            obj.Plot_title = {'ISO16533', str1 + "% to " + str2 + "% RH Step Change", obj.Material_Type};
            xlabel('Time (s)', "Interpreter", "latex")
            ylabel({"Estimate Heat Flux ", '$\frac{W}{m^2}$'}, 'Interpreter',"latex")
            
%             title({'ISO16533',"15% to 85% RH Step Change" ,obj.Material_Type})
            title(obj.Plot_title)
           
            set(gca,'fontsize', 18)
            
        end 
        
        function Create_Basline(obj, starting_idx, step_change_threshold)
            
            for i=1:obj.file_count
                Phi_BaseLined = obj.Replicate(i).Data.Phi;
                smooth_delt = obj.Replicate(i).Data.delta_Temp_smooth;
                flux_step_idx = ISO_Replicate.baseliner(smooth_delt, starting_idx, step_change_threshold);
                ave_steady_state_flux = mean(Phi_BaseLined(1:flux_step_idx-1));
                baseline_flux = Phi_BaseLined - ave_steady_state_flux;
                obj.Replicate(i).Data.Phi_BaseLined = baseline_flux;
                obj.Replicate(i).Data.Phi_step_idx = flux_step_idx;
                
                obj.Replicate(i).Data.integrate_flux
                
               
            end 
        
        end 
        
        function Plot_Baseline(obj)
            
            figure(1);clf;
            hold on
            for i=1:obj.file_count
               idx = obj.Replicate(i).Data.Phi_step_idx;
               plot(obj.Replicate(i).Data.Time(2:end),obj.Replicate(i).Data.Phi_BaseLined, "LineWidth", 1.75, "Color", obj.colors{i})
               plot(obj.Replicate(i).Data.Time(idx + 1:end), obj.Replicate(i).Data.Phi_Integrated_Reference(idx + 1:end), "Color", obj.colors{i}, "LineStyle" ,':',"LineWidth", 1.75 ,"DisplayName", 'Baseline')
            end 
            xlabel('Time (s)', "Interpreter", "latex")
            ylabel({"Baselined Estimate Heat Flux ", '$\frac{W}{m^2}$'}, 'Interpreter',"latex")
            
%             title({'ISO16533',"15% to 85% RH Step Change" ,obj.Material_Type})
%             title({'ISO16533', str1 + "% to " + str2 + "% RH Step Change", obj.Material_Type})
            title(obj.Plot_title)
            set(gca,'fontsize', 18)
            hold off
            
            
        end 
        
        function Compile_Results(obj)
            
            compiled_E = zeros(obj.file_count, 1);
            for i=1:obj.file_count
                compiled_E(i, 1) = obj.Replicate(i).Data.Phi_Integrated_Area;
            end 
            
            E_Table = array2table(compiled_E,'VariableNames', {'Integrated_Enegery'});
            obj.Compiled_Integration_Results = E_Table;  
            %  UNITS of kJ/m^2
            
        end 
        
        

        
    end
end

