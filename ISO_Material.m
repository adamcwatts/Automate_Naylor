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
        regain_EQ
        Energy_Theoretical
        R_cf
        file_count
        Replicate
        Plot_title
        Compiled_Integration_Joules_per_Meter_Squared
        Compiled_Energy_Joules
        Compiled_Average_Temp
        Compiled_Baseline_Temp
        Compiled_Baseline_Delta_Temp
        Compiled_STD_Temp
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

            master_data_Normalized =  nan(n_max_replicates, n_classes);
            master_data_Energy = nan(n_max_replicates, n_classes);
            master_data_theory = nan(1, n_classes);
            
            for k = 1:size(my_class_variables, 2)
                ob =  evalin('base', my_class_variables{k});
                master_data_Normalized(1:ob.file_count, k) = ob.Compiled_Integration_Joules_per_Meter_Squared{:,1};  
                master_data_Energy(1:ob.file_count, k) = ob.Compiled_Energy_Joules{:,1};  
                master_data_theory(1,k) = ob.Energy_Theoretical;
            end 
            
            T_Norm = array2table(master_data_Normalized,'VariableNames', my_class_variables);
            T_Energy = array2table(master_data_Energy,'VariableNames', my_class_variables);
            T_Theory = array2table(master_data_theory,'VariableNames', my_class_variables);
            
            fn_norm = "Integrated_Flux_Energy_Normalized_" + "ISO_METHOD" + ".csv";
            fn_energy = "Integrated_Flux_Energy_Total_" + "ISO_METHOD" + ".csv";
            fn_theory = "Integrated_Flux_Energy_Theory_" + "ISO_METHOD" + ".csv";
            
            writetable(T_Norm,fn_norm);
            disp("Exported: " + fn_norm)
            
            writetable(T_Energy,fn_energy);
            disp("Exported: " + fn_energy)
            
            writetable(T_Theory,fn_theory);
            disp("Exported: " + fn_theory)
            
            
        end 
        
        function Plot_Ave_Material_Temps()
            % call whos function in the base workspace 
            s = evalin('base', 'whos'); 
            
            % find the objects with the type 'ISO_Material' from the workspace:
            matches= strcmp({s.class}, 'ISO_Material');
            my_class_variables = {s(matches).name};
            n_classes = length(my_class_variables);
            
            ob =  evalin('base', my_class_variables{1});
            Guide_Temp = mean(ob.Compiled_Average_Temp(1:500));
            
            figure(99);clf;
            hold on           
            for i=1:n_classes
                ob =  evalin('base', my_class_variables{i});
                
                mean_basline = mean(ob.Compiled_Average_Temp(1:500));
                
                if  mean_basline >= Guide_Temp 
                    
                    Compiled_Average_Temp = ob.Compiled_Average_Temp - (mean_basline - Guide_Temp);
                    
                elseif mean_basline < Guide_Temp 
                    Compiled_Average_Temp = ob.Compiled_Average_Temp + (Guide_Temp - mean_basline );
                end 
                
                ob.Compiled_Baseline_Temp = Compiled_Average_Temp;
                ob.Compiled_Baseline_Delta_Temp = Compiled_Average_Temp - Guide_Temp;
                plot(ob.Replicate(1).Data.Time,  Compiled_Average_Temp - Guide_Temp, ...
                    "LineWidth", 2.25, "DisplayName", ob.Material_Type)
                
            end
            xlabel('time [s]','Interpreter',"latex") 
            ylabel("$\Delta \  \mathrm{Temperature} \ [^{\circ}\mathrm{C}]$", "Interpreter", "latex")
            legend()
            
            
            figure(100);clf;
            hold on           
            for i=1:n_classes
               ob =  evalin('base', my_class_variables{i});
               plot(ob.Replicate(1).Data.Time,  ob.Compiled_Baseline_Temp, ...
                    "LineWidth", 2.25, "DisplayName", ob.Material_Type)
            end 
            xlabel('time [s]','Interpreter',"latex") 
            ylabel("$\mathrm{Temperature} \ [^{\circ}\mathrm{C}]$", "Interpreter", "latex")
            legend()
%             legend(["Viscose (Thin)", "Wool (Thick", "Wool (Thin)", "Cotton (Thin)", "Polyester (Thick)", "Polyester (Thin)"])
            
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
        
        function Material_Data(obj, Type ,Square_Length_mm, thickness_mm, Weight_Density, cp_fabric, R_ct)
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
            obj.R_cf = R_ct; 
        end 
        
        function Plot_Smoothing_Delta_T(obj)
        
            for i=1:obj.file_count
                figure(10+i); clf;  
                plot(obj.Replicate(i).Data.Time(2:end,1), obj.Replicate(i).Data.delta_Temp_raw(:,:), "Color", '#787878')
                hold on
                plot(obj.Replicate(i).Data.Time(2:end,1), obj.Replicate(i).Data.delta_Temp_smooth(:,:),'-d',"LineWidth",2, "Color", obj.colors{i},...
                    'MarkerIndices',1:200:length(obj.Replicate(i).Data.delta_Temp_smooth(:,:)), 'MarkerSize',7)
                
                xlabel('time [s]','Interpreter',"latex")
                ylabel('$\frac{dT}{dt} \ \ [\frac{^{\circ}\mathrm{C}}{s}]$','Interpreter', "latex")
                legend({'Raw Delta', "Gaussian Smoothed"})
                title({'Differential Smoothing' ,obj.Replicate(i).Data.file_name}, 'Interpreter', 'none')
                hold off
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
            xlabel({'I'},'Interpreter',"latex")
            ylabel('Temperature [$^{\circ}\mathrm{C}$]', 'Interpreter', "latex")
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
            xlabel('II','Interpreter',"latex")
            ylabel('$\frac{dT}{dt} \ \ [\frac{^{\circ}\mathrm{C}}{s}]$','Interpreter', "latex")
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
            xlabel({'III' ,'time [s]'},'Interpreter',"latex")
            ylabel('RH [\%]', 'Interpreter', "latex")
            legend('Interpreter', 'none', 'Location', "SouthEast")
            hold off
            
            
            
        end 
        
        function Impose_Temp(obj)
            
            first_rep_step_idx  = obj.Replicate(1).Data.Phi_step_idx;
            
            Guide_Temp = mean(obj.Replicate(1).Data.Temp(1:first_rep_step_idx));
            obj.Replicate(1).Data.Temp_Baselined = obj.Replicate(1).Data.Temp;
            
            for i=2:obj.file_count
               T_step_idx =  obj.Replicate(i).Data.Phi_step_idx;
               
               mean_basline = mean(obj.Replicate(i).Data.Temp(1:T_step_idx));
               
               if mean_basline > Guide_Temp
                obj.Replicate(i).Data.Temp_Baselined = obj.Replicate(i).Data.Temp - (mean_basline - Guide_Temp) ;
               
               elseif mean_basline < Guide_Temp
                 obj.Replicate(i).Data.Temp_Baselined = obj.Replicate(i).Data.Temp + (Guide_Temp - mean_basline);
               end  
               
            end 
            
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
            xlabel('time [s]', "Interpreter", "latex")
            ylabel({"Estimate Heat Flux ", '$[\frac{W}{m^2}]$'}, 'Interpreter',"latex")
            
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
                obj.Replicate(i).Data.integrate_flux(); 
                obj.Replicate(i).Data.Energy_Joules = obj.Replicate(i).Data.Phi_Integrated_Area * obj.Area_m;  % j
                                
               
            end 
        
        end 
        
        function Plot_Baseline(obj, xlimit, ylimit)
            
            figure(1);clf;
            hold on
            for i=1:obj.file_count
               idx = obj.Replicate(i).Data.Phi_step_idx;
               plot(obj.Replicate(i).Data.Time(2:end),obj.Replicate(i).Data.Phi_BaseLined, "LineWidth", 1.75, "Color", obj.colors{i})
               plot(obj.Replicate(i).Data.Time(idx + 1:end), obj.Replicate(i).Data.Phi_Integrated_Reference(idx + 1:end), "Color", obj.colors{i}, "LineStyle" ,':',"LineWidth", 1.75 ,"DisplayName", 'Baseline')
            end 
            xlabel('time [s]', "Interpreter", "latex")
            ylabel({"Baselined Estimate Heat Flux ", '$[\frac{W}{m^2}]$'}, 'Interpreter',"latex")
            
            if ~isempty(xlimit)
                xlim([xlimit(1), xlimit(2)])
            end 
            
            if ~isempty(ylimit)
                ylim([ylimit(1), ylimit(2)])
            end 
%             title({'ISO16533',"15% to 85% RH Step Change" ,obj.Material_Type})
%             title({'ISO16533', str1 + "% to " + str2 + "% RH Step Change", obj.Material_Type})
            title(obj.Plot_title)
            set(gca,'fontsize', 18)
            hold off
            
            
        end 
        
        function Calculate_Theoretical_Energy(obj, regain_EQ, RH_initial, RH_final)
            obj.regain_EQ = regain_EQ;
            % ALL arguments assumes %
            regain_EQ = regain_EQ /100;
            RH_initial = RH_initial / 100;
            RH_final = RH_final / 100;
            
            H_vap = 2418 ; % j/g
            H_sorp = @(rh) 195 .* (1.0 - rh) .* ((1.0 ./ (0.2 + rh)) + (1.0 ./ (1.05 - rh)));
            
            R = @(RH) 0.578 * RH * (1/(0.321 + RH) + 1/(1.262-RH));
            delta_regain = R(RH_final) - R(RH_initial);
            mass_vapor = delta_regain * regain_EQ * obj.Area_m * obj.Weight_Density_g_m_sq;  % grams
            
            H = H_vap + integral(H_sorp, RH_initial, RH_final); % joules/gram
            
            obj.Energy_Theoretical = mass_vapor * H;  % joules
        end 
        
        
        function Compile_Results(obj)
            % Asumes all replicates have the same length of data
            Compiled_Base_Temp = zeros(length(obj.Replicate(1).Data.Temp), obj.file_count);
            
            compiled_E_Norm = zeros(obj.file_count, 1);
            complied_E_Total = zeros(obj.file_count, 1);
            
            for i=1:obj.file_count
                compiled_E_Norm(i, 1) = obj.Replicate(i).Data.Phi_Integrated_Area;
                complied_E_Total(i, 1) = obj.Replicate(i).Data.Energy_Joules;
                Compiled_Base_Temp(:,i) = obj.Replicate(i).Data.Temp_Baselined;
            end 
            
            E_Norm_Table = array2table(compiled_E_Norm,'VariableNames', {'Integrated_Enegery_Normalized_Area'});
            obj.Compiled_Integration_Joules_per_Meter_Squared = E_Norm_Table; 
            % UNITS J/m^2
            
            E_Total_Table = array2table(complied_E_Total,'VariableNames', {'Integrated_Enegery_Joules'});
            obj.Compiled_Energy_Joules = E_Total_Table;
            %  UNITS of J
            
            obj.Compiled_STD_Temp = std(Compiled_Base_Temp,[], 2);
            obj.Compiled_Average_Temp = mean(Compiled_Base_Temp, 2);
            


        end 
        
        

        
    end
end

