classdef Naylor_Material < handle  % objects by reference, handle class 
    %{
    Written by Adam C. Watts
    5/26/2020 : Added export results() static method for easy export and compilation of integrated baseline results 
    5/19/2020 : Completed Class
    %}
    
    properties
        Replicate
        material_info
        Area_m
        Weight_Density_g_m_sq
        regain_EQ
        Energy_Theoretical
        file_count   % file count
        Compiled_Model_Parameters
        Compiled_Predictions_Integrated
        Compiled_Predictions_Energy_Joules
        baseline_flux  % Using 1st replicate
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
            
            % find the objects with the type 'Naylor_Material' from the workspace:
            matches= strcmp({s.class}, 'Naylor_Material');
            my_class_variables = {s(matches).name};
            n_classes = length(my_class_variables);
            
   
            n_max_replicates = 0;
            all_method_names = {};
            
            % Find how many unique function methods was used and maximum
            % file count for each class instance
            for var = my_class_variables
                ob =  evalin('base', var{1});
                method_names = ob.Compiled_Predictions_Integrated.Properties.VariableNames;
                size_of_names = size(method_names, 2);
  
                if size_of_names > size(all_method_names, 2)
                      all_method_names = cat(2,all_method_names, method_names);
                      all_method_names = unique(all_method_names);
                end 
                
                if ob.file_count > n_max_replicates
                    n_max_replicates = ob.file_count;
                    
                end 
            end 
            
            n_methods = length(all_method_names);
            master_data_Normalized =  nan(n_max_replicates, n_classes, n_methods);
            master_data_E = nan(n_max_replicates, n_classes, n_methods);
            master_data_theory = nan(n_methods, n_classes);

            
            for i = 1:size(all_method_names,2)
                for k = 1:size(my_class_variables, 2)
                    ob =  evalin('base', my_class_variables{k});
                    mask = cellfun(@(s) contains(all_method_names{i}, s), ob.Compiled_Predictions_Integrated.Properties.VariableNames);
                    
                    if any(mask)
                        master_data_Normalized(1:ob.file_count, k, i) = ob.Compiled_Predictions_Integrated.(all_method_names{i});
                        master_data_E(1:ob.file_count, k, i) = ob.Compiled_Predictions_Energy_Joules.(all_method_names{i});
                        master_data_theory(i, k) = ob.Energy_Theoretical;   
                    end
                  
                end 
                
                % Create table with for each function method using all replicates and class instances 
                T_Norm = array2table(master_data_Normalized(1:ob.file_count, :, i),'VariableNames', my_class_variables);
                fn_Norm = "Integrated_Flux_Area_Normalized_" + all_method_names{i} + "_NAYLOR.csv";
                writetable(T_Norm, fn_Norm);
                disp("Exported: " + fn_Norm)
                
                T_E =  array2table(master_data_E(1:ob.file_count, :, i),'VariableNames', my_class_variables);
                fn_E = "Integrated_Flux_Energy_" + all_method_names{i} + "_NAYLOR.csv";
                writetable(T_E, fn_E);
                disp("Exported: " + fn_E)
                
                T_theory =  array2table(master_data_theory(i, :),'VariableNames', my_class_variables);
                fn_theory = "Integrated_Flux_Energy_Theory_" + all_method_names{i} + "_NAYLOR.csv";
                writetable(T_theory, fn_theory);
                disp("Exported: " + fn_theory)
                
                
                
            end 
            
             
        end 
        
        function [RH_step_and_steady_state_idx, RH_step_idx] = step_changer(Data, starting_idx, ss_index, smooth_opt, span)
            %{

        Arg2: starting_idx: occasionaly largest derivative is not actual time step and needs custom idx 
        Arg3: ss_index: how many idex steps to take prior to step change which acts as steady state time

        %}

        step_change = 0.02;  % THRESHOLD DIFF

        if isempty(span)
            span = 5;
        end 

        smoothed_diff_data = abs(diff(smooth(Data, span, 'moving')));
        % smoothed_diff_data = abs(smoothdata(diff(Data(starting_idx:end, 1))));  %4th column R
        diff_data = abs(diff(Data(starting_idx:end,1)));
        % figure(1); clf;
        % plot(Data(1 + starting_idx:end,1), smoothed_diff_data);

        [~, max_idx_smooth] = max(smoothed_diff_data);
        [~, max_idx] = max(diff_data);

        reversed_data_smooth = flip(smoothed_diff_data(1:max_idx_smooth));
        reversed_data = flip(diff_data(1:max_idx));

        if smooth_opt
            reversed_data = reversed_data_smooth;
            max_idx = max_idx_smooth;
        end

        step_idx = 0;
        for j=1:length(reversed_data)

            if abs(reversed_data(j)) < step_change
               step_idx = (j - 1) ;
               break
            end
        end
        
%         local_step_idx = max_idx - step_idx;  % local to starting idx being truncated  (already takes into accout of starting index)
%         local_step_idx = [];
        
        RH_step_idx = max_idx - step_idx + (starting_idx);  % when starting idx is 1 and greater. (corrects idx for double peaks)
        RH_step_and_steady_state_idx = RH_step_idx - ss_index;  % number of steps to go back before step change
        
        
%         global_step_idx = max_idx - step_idx;  % relative to the entire array
%         
%         local_step_idx = max_idx - step_idx + (starting_idx - 1);  % when starting idx is 1 and greater. (corrects idx for double peaks)
%         local_step_idx = local_step_idx - ss_index;  % number of steps to go back before step change
        end 
        
    
        function up_flux_vals =  inverter(HeatFlux, baseline_flux, ss_step_index)
            if isempty(baseline_flux)
                baseline_flux = mean(HeatFlux(1:ss_step_index));
            end 
    
        up_flux_vals = (-1*HeatFlux(ss_step_index:end)) + baseline_flux;
        
        
        end 
        
    end 
    
    methods
        % Intialize material with naylor replicate class 
        function obj = Naylor_Material(file_path_locations_of_reps, length_cm)
            FPLR = file_path_locations_of_reps;
            obj.file_count = length(FPLR); % File Count
            obj.Area_m = length_cm^2 * (1/100)^2;
           % store replicate objects into structure array 
            for i=1:obj.file_count
                obj.Replicate(i).Data = Naylor_Replicate(FPLR{i});
            end 
            
        end 
        
      
        function Plot_Raw_Data(obj)
            % Time vs Heat Flux Plots
            figure(1); clf;
            hold on
            for i =1:obj.file_count
            plot(obj.Replicate(i).Data.Time, ...
            obj.Replicate(i).Data.HeatFlux, ...
            "LineWidth",1.75,"DisplayName",....
            obj.Replicate(i).Data.file_name)
        
            end 
            xlabel('time [s]','Interpreter',"latex")
            ylabel({'Heat Flux due to Sorption' '$[\frac{W}{m^2}]$'}, 'Interpreter',"latex")
            legend()
            title("Raw Heat Flux vs Time Data")
            hold off

            % Time vs RH Plots
            figure(2); clf;
            hold on
            for i =1:obj.file_count
            plot(obj.Replicate(i).Data.Time,...
                obj.Replicate(i).Data.RH,...
                "LineWidth",1.75, "DisplayName",...
                obj.Replicate(i).Data.file_name)
            end 
            xlabel('time [s]','Interpreter',"latex")
            ylabel('Relative Humidity [\%]', 'Interpreter',"latex")
            legend('Location','SouthEast')
            title("Raw RH vs Time Data")
            hold off

        end 
        
        function Update_Starting_Index(obj,replicate_number, updated_time)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            n = replicate_number;
            [~, closest_index] = min(abs(obj.Replicate(n).Data.Time - updated_time));
           obj.Replicate(n).Data.starting_index = closest_index;
                      
        end
        
        function RH_Indices_And_Plot(obj)
        %UNTITLED4 Summary of this function goes here
        %   Detailed explanation goes here
        RH_step_and_steady_state_idx = zeros(1, obj.file_count);
        RH_step_idx = zeros(1, obj.file_count);

        for i =1:obj.file_count 
             [RH_step_and_steady_state_idx(1,i), RH_step_idx(1,i)] = Naylor_Material.step_changer(obj.Replicate(i).Data.RH,  obj.Replicate(i).Data.starting_index , obj.Replicate(i).Data.steady_state_index, obj.Replicate(i).Data.option_smooth_import_data, []);
             obj.Replicate(i).Data.RH_step_idx = RH_step_idx(1,i);
             obj.Replicate(i).Data.RH_step_and_steady_state_idx = RH_step_and_steady_state_idx(1,i);
             obj.Replicate(i).Data.RH_step_time = obj.Replicate(i).Data.Time(RH_step_idx(1,i), 1) ;
             obj.Replicate(i).Data.RH_step_and_steady_state_time = obj.Replicate(i).Data.Time(RH_step_and_steady_state_idx(1,i), 1) ;
        end 

        figure(3); clf;
        hold on
        for i =1:obj.file_count
            plot(obj.Replicate(i).Data.Time(RH_step_and_steady_state_idx(1, i):end, 1), obj.Replicate(i).Data.RH(RH_step_and_steady_state_idx(1, i):end,1),"LineWidth",1.75, "DisplayName", obj.Replicate(i).Data.file_name)
 
        end 
        xlabel('time [s]','Interpreter',"latex")
        ylabel('Relative Humidity [\%]', 'Interpreter',"latex")
        legend('Location','SouthEast')
        
        ss_time = obj.Replicate(1).Data.starting_time_minutes  ;
        str = sprintf('%d', ss_time);
        title({"Truncated RH Curves", "Steady-State Time of " + str + " Minutes"}) 
        hold off 

        end 
        
        function Create_Baseline(obj)
        % IMPOSED FLUX DATA CREATION
        rh_idx = obj.Replicate(1).Data.RH_step_and_steady_state_idx;
        ss_idx = obj.Replicate(1).Data.steady_state_index;

        % Use first replicate as the official baseline heat flux
        base_line_flux = obj.Replicate(1).Data.HeatFlux(rh_idx:end, 1); % truncate beginning 
        base_line_flux = mean(base_line_flux(1:ss_idx));  % average baseline to ss_idx


            for i =1:obj.file_count
            rh_idx = obj.Replicate(i).Data.RH_step_and_steady_state_idx; 
            ss_index = obj.Replicate(i).Data.steady_state_index;
            time = obj.Replicate(i).Data.Time(rh_idx:end, 1);
            
            time_0 = time(1);
            baseline_time = time - time_0;
            flux = obj.Replicate(i).Data.HeatFlux(rh_idx:end, 1);   % truncate beginning 

                if i > 1
                    ss_flux = mean(flux(1:ss_index)); % average baseline to prior_idx

                    if ss_flux > base_line_flux
                        adjust_base_line = ss_flux - base_line_flux;
                        flux = flux - adjust_base_line;

                    elseif ss_flux <= base_line_flux
                        adjust_base_line = base_line_flux - ss_flux;
                        flux = flux + adjust_base_line;
                    end
                end

            baseline_temp = obj.Replicate(i).Data.Temp(rh_idx:end, 1);
            baseline_RH = obj.Replicate(i).Data.RH(rh_idx:end, 1);
            pass_data = [baseline_time, flux, baseline_temp, baseline_RH];

            obj.Replicate(i).Imposed_Data =  Naylor_Replicate_Basic(pass_data, obj.Replicate(i).Data.file_name, obj.Replicate(i).Data.steady_state_index);
            
            end 
        obj.baseline_flux = base_line_flux;
        end 
        
        function Plot_Basline(obj)
    
        figure(10); clf;
        hold on
        for i =1:obj.file_count
        plot(obj.Replicate(i).Imposed_Data.Time, obj.Replicate(i).Imposed_Data.HeatFlux,"LineWidth",1.75,"Color", obj.colors{i}, "DisplayName", obj.Replicate(i).Imposed_Data.file_name)
        end 
        xlabel('time [s]','Interpreter',"latex")
        ylabel({'Heat Flux due to Sorption' '$[\frac{W}{m^2}]$'}, 'Interpreter',"latex")
        legend()
        plot(obj.Replicate(i).Imposed_Data.Time, obj.baseline_flux*ones(1,length(obj.Replicate(i).Imposed_Data.Time)), "Color", 'k',"LineStyle" ,':',"LineWidth", 1.75 ,"DisplayName", 'Baseline')
        title("Superimposed Heat Flux vs Time")
        hold off


        figure(11); clf;
        hold on
        for i =1:obj.file_count
        plot(obj.Replicate(i).Imposed_Data.Time, obj.Replicate(i).Imposed_Data.RH, "LineWidth",1.75,"Color",obj.colors{i}, "DisplayName", obj.Replicate(i).Imposed_Data.file_name)
        end 
        xlabel('time [s]','Interpreter',"latex")
        ylabel('Relative Humidity [\%]', 'Interpreter',"latex")
        legend('Location','SouthEast')
        title("Superimposed RH vs Time")
        hold off
      
        end

        function Fit(obj,n_guesses, fit_type)
            
             if isempty(fit_type)
                fit_type = obj.model_types{1};
            end 
            
            for i=1:obj.file_count
                Time = obj.Replicate(i).Imposed_Data.Time;
                obj.Replicate(i).(fit_type) = Convolution_Model(n_guesses, fit_type, Time);  % Instantiate Model Predictor
                obj.Replicate(i).(fit_type).Stochastic_Solver(obj.Replicate(i).Imposed_Data) % Runs Stochastic Solver
            end 
            
        end 
        
        function models_used = Check_Model(obj)
          models_used = {'None'};
          all_field_names = fieldnames(obj.Replicate)';
          
          if length(all_field_names) <= 2
              error("Error: No Predictions were made. Please use the Predict Class Method")  
          end 
          
          all_models = Convolution_Model.model_list();
          
          for k=3:length(all_field_names)
            
              mask = cellfun(@(s) ~isempty(strfind(all_field_names{k}, s)), all_models); % check whether the field item is in the list of all models
              
              if any(mask)
                name_cell =  all_models(mask);
                models_used{k-2} = name_cell{1}; % returns the model
              end 
          end 
          
          
        end 
        
        function Plot_Convolution(obj)

          models_used = Check_Model(obj);
            
         for i=1:obj.file_count
            Time = obj.Replicate(i).Imposed_Data.Time;
            HeatFlux = obj.Replicate(i).Imposed_Data.HeatFlux;
            RH = obj.Replicate(i).Imposed_Data.RH;
%             BEST_HF_MODEL = obj.Replicate(i).Imposed_Data.Model_Prediction;
            file_name = obj.Replicate(i).Imposed_Data.file_name;
       
            
            figure(i+10);clf
            subplot(2,1,1);
            hold on
            plot(Time, HeatFlux,'r',"DisplayName", "Raw Data: " + file_name,"LineWidth",1.5,"Color", obj.colors{1})
            
            for ii = 1:length(models_used)
              plot(Time(1:end-1), obj.Replicate(i).(models_used{ii}).Model_Fit, "DisplayName", "$\hat{\Psi}(t,$\boldmath$\theta$): " + obj.Replicate(i).(models_used{ii}).Model_Choice , "Color" ,obj.colors{ii+1}  ,"LineWidth", 1.8, "LineStyle", ":")
%               plot(Time(1:end-1), obj.Replicate(i).(models_used{ii}).Model_Fit, "DisplayName",  obj.Replicate(i).(models_used{ii}).Model_Choice, "Color" ,obj.colors{ii+1}  ,"LineWidth", 1.8, "LineStyle", ":")

            end 
            xlabel('I',"Interpreter","latex")
            ylabel('Heat Flux $[\frac{W}{m^2}]$',"Interpreter","latex")
            lgd1 = legend("Interpreter","latex");
            lgd1.Location = 'southeast';
%             formatSpec = "Best Fitted Model for " + obj.Replicate(i).Imposed_Data.file_name + "\nUsing %d Stochastic Initial Guesses";
%             A1 = obj.Replicate(i).Imposed_Data.Initial_Guesses;
%             str = sprintf(formatSpec,A1);
            title({file_name, "Convolution Model Results"})
            hold off

            subplot(2,1,2);
            plot(Time,RH,"Color",obj.colors{end},"LineWidth",2,"DisplayName","Relative Humidity [\%]")
            xlabel({"II", 'time [s]'},"Interpreter","latex")
            ylabel('RH [\%]',"Interpreter","latex")
            lgd2 = legend("Interpreter","latex");
            lgd2.Location = 'southeast';
         end 
        end
        
        function Predict(obj)
            
            models_used = Check_Model(obj);
            n_models = length(models_used);

            for i=1:obj.file_count
                model_predict_flux = zeros(length(obj.Replicate(i).Imposed_Data.Time) - 1 , n_models);
                model_pred_inverted_flux = zeros(length(obj.Replicate(i).Imposed_Data.Time) - 1, n_models);
                baseline_results = zeros(length(obj.Replicate(i).Imposed_Data.Time) - 1, n_models);
                adjusted_area = zeros(1, n_models);
                
                % New Instance of class
                pass_data = [obj.Replicate(i).Imposed_Data.Time, obj.Replicate(i).Imposed_Data.HeatFlux, obj.Replicate(i).Imposed_Data.Temp, obj.Replicate(i).Imposed_Data.RH];
                obj.Replicate(i).Predicted_Data = Naylor_Replicate_Basic(pass_data, obj.Replicate(i).Data.file_name, obj.Replicate(i).Data.steady_state_index);
                
                try 
                obj.Replicate(i).Predicted_Data.addprop('Integration_Baseline');
                obj.Replicate(i).Predicted_Data.addprop('Integration_values');
                obj.Replicate(i).Predicted_Data.addprop('HeatFlux_Model_Predictions');
                obj.Replicate(i).Predicted_Data.addprop('HeatFlux_Model_Predictions_Inverted');
                obj.Replicate(i).Predicted_Data.addprop('RH_STEP');
                catch 
                     
                end   
                    
                n_original_rep = length(obj.Replicate(i).Imposed_Data.Time);
                
                ss_index = obj.Replicate(i).Imposed_Data.steady_state_index + 1; % model now starts at t=0 not t=10 seconds, needs extra index
                
                RH = obj.Replicate(i).Imposed_Data.RH;
                Temp = obj.Replicate(i).Imposed_Data.Temp;
                new_RH_ss = mean(RH(1:ss_index,1)) * ones(ss_index, 1);
                new_RH_step = max(RH) * ones(n_original_rep - ss_index, 1);
                new_RH = [new_RH_ss; new_RH_step]; % concatenate vertically 
                
                obj.Replicate(i).Predicted_Data.RH_STEP = new_RH;
                
                RHdiff = diff(new_RH);
                original_HeatFlux = obj.Replicate(i).Imposed_Data.HeatFlux;  
                HF0 = mean(original_HeatFlux(1:ss_index));
                
                for ii=1:length(models_used)
                    Time = obj.Replicate(i).Imposed_Data.Time(1:end-1);
                    model_ob = obj.Replicate(i).(models_used{ii});
                    params = model_ob.Best_Convolution_Model_Params;
                    
        % overwrite naylor function with the time domain for the prediction
                    model_ob.Naylor_Function = @(params)model_ob.Naylor_Base_Function(params, Time); 
                    predicted_flux =  model_ob.naylor_curve(params, Temp, ss_index, HF0, RHdiff ); 
                    
                    
                    model_predict_flux(:, ii) = predicted_flux;
                   
                    inverted_pred_flux = Naylor_Material.inverter(predicted_flux, [], 1);  % dont pass baseline, use index of 1
                    model_pred_inverted_flux(:, ii) = inverted_pred_flux;
                    
                    new_step_index = ss_index;
                    
%                     [~, new_step_index] = Naylor_Material.step_changer(inverted_pred_flux, 1, ss_index, false, []);
                    [baseline_results(:, ii), adjusted_area(1, ii)] = Convolution_Model.integrate_flux(inverted_pred_flux, Time, new_step_index);  
                    
                    
                end 
                obj.Replicate(i).Predicted_Data.HeatFlux_Model_Predictions = array2table(model_predict_flux,'VariableNames', models_used);
                obj.Replicate(i).Predicted_Data.HeatFlux_Model_Predictions_Inverted = array2table(model_pred_inverted_flux, 'VariableNames', models_used);
                obj.Replicate(i).Predicted_Data.Integration_Baseline = array2table(baseline_results, 'VariableNames', models_used);
                obj.Replicate(i).Predicted_Data.Integration_values = array2table(adjusted_area, 'VariableNames', models_used);
                                 
            end 
            
        end 
        
        function Plot_Prediction(obj)
            models_used = Check_Model(obj);
            n_models = length(models_used);
            corrected_model_names = cell(1, n_models);
            
            for i=1:n_models
                corrected_model_names{i} = strrep(models_used{i}, "_", " ");
            end 
            
            for i=1:obj.file_count 
                Time = obj.Replicate(i).Predicted_Data.Time;
                RH_STEP = obj.Replicate(i).Predicted_Data.RH_STEP;
    %             BEST_HF_MODEL = obj.Replicate(i).Imposed_Data.Model_Prediction;
                file_name = obj.Replicate(i).Imposed_Data.file_name;
                Predicted_HeatFlux = table2array(obj.Replicate(i).Predicted_Data.HeatFlux_Model_Predictions);
                Predicted_HeatFlux_Inverted = table2array(obj.Replicate(i).Predicted_Data.HeatFlux_Model_Predictions_Inverted);
                Y_baselines = table2array(obj.Replicate(i).Predicted_Data.Integration_Baseline);
                ss_index = obj.Replicate(i).Imposed_Data.steady_state_index + 1;

                figure(i+10);clf
                subplot(2,1,1);
                hold on

                for ii = 1:n_models
%                     plot(Time(1:end-1), Predicted_HeatFlux(:,ii), "DisplayName", "$\hat{\Psi}(t,$\boldmath$\theta$): Predicted" + corrected_model_names{ii}, "Color" ,obj.colors{ii+1}  ,"LineWidth", 1.8, "LineStyle", ":")
                    plot(Time(1:end-1), Predicted_HeatFlux(:,ii), "DisplayName", "$\Psi(t,$\boldmath$\theta$): Predicted with " + corrected_model_names{ii}, "Color" ,obj.colors{ii+1}  ,"LineWidth", 1.8, "LineStyle", ":")

                end 
                xlabel('I',"Interpreter","latex")
                ylabel('Heat Flux $[\frac{W}{m^2}]$',"Interpreter","latex")
%                 ylabel("$\hat{\Psi}(t,$\boldmath$\theta$): Predicted","Interpreter","latex")
                lgd1 = legend("Interpreter", 'latex');
                lgd1.Location = 'southeast';
                title({file_name, "Predicted Step-Change Convolution Model Results"})
                hold off

                subplot(2,1,2);
                plot(Time,RH_STEP,"Color",obj.colors{end},"LineWidth",2,"DisplayName","Relative Humidity [\%]")
                xlabel({"II", 'time [s]'},"Interpreter","latex")
                ylabel('RH [\%]',"Interpreter","latex")
                lgd2 = legend("Interpreter","latex");
                lgd2.Location = 'southeast';
                
                
                figure(i+20);clf
                subplot(2,1,1);
                hold on

                for ii = 1:n_models
                    plot(Time(1:end-1), Predicted_HeatFlux_Inverted(:,ii), "DisplayName", "$\Psi(t,$\boldmath$\theta$): Predicted with " + corrected_model_names{ii}, "Color" ,obj.colors{ii+1}  ,"LineWidth", 1.8, "LineStyle", "-")
                    plot(Time(ss_index:end-1), Y_baselines(ss_index: end, ii),'LineStyle', ':', "Color", 'k',"LineWidth", 2,"DisplayName","Integrating Baseline")
                end 
               
                xlabel('I',"Interpreter","latex")
                ylabel({'Heat Flux',  ' $(Exotherm \rightarrow) \  [\frac{W}{m^2}]$'},"Interpreter","latex")
% % % %                 ylabel("$\hat{\Psi}(t,$\boldmath$\theta$): Predicted","Interpreter","latex")
                lgd1 = legend("Interpreter", 'latex');
                lgd1.Location = 'southeast';
                title({file_name, "Predicted Step-Change Convolution Model Results"})
                hold off

                subplot(2,1,2);
                plot(Time,RH_STEP,"Color",obj.colors{end},"LineWidth",2,"DisplayName","Relative Humidity [\%]")
                xlabel({"II", 'time [s]'},"Interpreter","latex")
                ylabel('RH [\%]',"Interpreter","latex")
                lgd2 = legend("Interpreter","latex");
                lgd2.Location = 'southeast';
                
                
            end 
        end
        
        function Calculate_Theoretical_Energy(obj, regain_EQ, RH_initial, RH_final, Weight_Density_g_m_sq)
            obj.Weight_Density_g_m_sq = Weight_Density_g_m_sq;
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
           models_used = Check_Model(obj);
           n_models = length(models_used); 
           
           model_param_counts = zeros(1, n_models);
    
           for i=1:length(models_used)
               model_param_counts(1, i) = obj.Replicate(1).(models_used{i}).n_params;
           end 
          column_count = sum(model_param_counts);  % total columns for parameters
           
%           param_array = zeros(obj.file_count, column_count); % instantiate empty array
          
          param_cell = cell(obj.file_count, n_models);
          predict_cell = zeros(obj.file_count, n_models);
          
          for i=1:obj.file_count
              predict_cell(i, :) = table2array(obj.Replicate(i).Predicted_Data.Integration_values);
              for ii=1:n_models
                  param_cell{i, ii} = obj.Replicate(i).(models_used{ii}).Best_Convolution_Model_Params;
                  
              end 
          end 
        
           obj.Compiled_Model_Parameters = cell2table(param_cell,'VariableNames', models_used);
           obj.Compiled_Predictions_Integrated = array2table(predict_cell,'VariableNames', models_used);
           obj.Compiled_Predictions_Energy_Joules = array2table( obj.Area_m * predict_cell,'VariableNames', models_used);
        end 
        
    end
end

