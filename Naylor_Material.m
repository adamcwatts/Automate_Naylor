classdef Naylor_Material < handle  % objects by reference, handle class 
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Replicate
        file_count   % file count
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
        function [local_step_idx, global_step_idx] = step_changer(Data, starting_idx, ss_index, smooth_opt, span)
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
               step_idx = j;
               break
            end
        end
        
%         local_step_idx = max_idx - step_idx;  % local to starting idx being truncated  (already takes into accout of starting index)
%         local_step_idx = [];
%         
%         global_step_idx = max_idx - step_idx + (starting_idx);  % when starting idx is 1 and greater. (corrects idx for double peaks)
%         global_impose_idx = global_step_idx - ss_index;  % number of steps to go back before step change
        
        
        global_step_idx = max_idx - step_idx;  % relative to the entire array
        
        local_step_idx = max_idx - step_idx + (starting_idx - 1);  % when starting idx is 1 and greater. (corrects idx for double peaks)
        local_step_idx = local_step_idx - ss_index;  % number of steps to go back before step change
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
        function obj = Naylor_Material(file_path_locations_of_reps)
            FPLR = file_path_locations_of_reps;
            obj.file_count = length(FPLR); % File Count
            
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
            xlabel('time (s)','Interpreter',"latex")
            ylabel({'Heat Flux due to Sorption' '$\frac{W}{m^2}$'}, 'Interpreter',"latex")
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
            xlabel('time (s)','Interpreter',"latex")
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
        RH_step_idx_Local = zeros(1, obj.file_count);
        RH_step_idx_Global = zeros(1, obj.file_count);

        for i =1:obj.file_count 
             [RH_step_idx_Local(1,i), RH_step_idx_Global(1,i)] = Naylor_Material.step_changer(obj.Replicate(i).Data.RH,  obj.Replicate(i).Data.starting_index , obj.Replicate(i).Data.steady_state_index, obj.Replicate(i).Data.option_smooth_import_data, []);
             obj.Replicate(i).Data.RH_step_idx_Local = RH_step_idx_Local(1,i);
             obj.Replicate(i).Data.RH_step_idx_Global = RH_step_idx_Global(1,i);
             obj.Replicate(i).Data.RH_step_time_Local = obj.Replicate(i).Data.Time(RH_step_idx_Local(1,i), 1) ;
             obj.Replicate(i).Data.RH_step_time_Global = obj.Replicate(i).Data.Time(RH_step_idx_Global(1,i), 1) ;
        end 

        figure(3); clf;
        hold on
        for i =1:obj.file_count
            plot(obj.Replicate(i).Data.Time(RH_step_idx_Local(1, i):end, 1), obj.Replicate(i).Data.RH(RH_step_idx_Local(1, i):end,1),"LineWidth",1.75, "DisplayName", obj.Replicate(i).Data.file_name)
 
        end 
        xlabel('time (s)','Interpreter',"latex")
        ylabel('Relative Humidity [\%]', 'Interpreter',"latex")
        legend('Location','SouthEast')
        title("Truncated Steady-State RH Curves")
        hold off

        end 
        
        function Create_Baseline(obj)
        % IMPOSED FLUX DATA CREATION
        rh_idx = obj.Replicate(1).Data.RH_step_idx_Local;
        ss_idx = obj.Replicate(1).Data.steady_state_index;

        % Use first replicate as the official baseline heat flux
        base_line_flux = obj.Replicate(1).Data.HeatFlux(rh_idx:end, 1); % truncate beginning 
        base_line_flux = mean(base_line_flux(1:ss_idx));  % average baseline to ss_idx


            for i =1:obj.file_count
            rh_idx = obj.Replicate(i).Data.RH_step_idx_Local; 
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
        xlabel('time (s)','Interpreter',"latex")
        ylabel({'Heat Flux due to Sorption' '$\frac{W}{m^2}$'}, 'Interpreter',"latex")
        legend()
        plot(obj.Replicate(i).Imposed_Data.Time, obj.baseline_flux*ones(1,length(obj.Replicate(i).Imposed_Data.Time)), "Color", 'k',"LineStyle" ,':',"LineWidth", 1.75 ,"DisplayName", 'Baseline')
        title("Superimposed Heat Flux vs Time")
        hold off


        figure(11); clf;
        hold on
        for i =1:obj.file_count
        plot(obj.Replicate(i).Imposed_Data.Time, obj.Replicate(i).Imposed_Data.RH, "LineWidth",1.75,"Color",obj.colors{i}, "DisplayName", obj.Replicate(i).Imposed_Data.file_name)
        end 
        xlabel('time (s)','Interpreter',"latex")
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
                obj.Replicate(i).(fit_type).Stochastic_Solver(obj.Replicate(i).Imposed_Data)
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
                plot(Time(1:end-1), obj.Replicate(i).(models_used{ii}).Model_Fit, "DisplayName",obj.Replicate(i).(models_used{ii}).Model_Choice, "Color" ,obj.colors{ii+1}  ,"LineWidth", 1.8, "LineStyle", ":")
            end 
            xlabel('Time (s)',"Interpreter","latex")
            ylabel('Heat Flux $\frac{W}{m^2}$',"Interpreter","latex")
            lgd1 = legend();
            lgd1.Location = 'southeast';
%             formatSpec = "Best Fitted Model for " + obj.Replicate(i).Imposed_Data.file_name + "\nUsing %d Stochastic Initial Guesses";
%             A1 = obj.Replicate(i).Imposed_Data.Initial_Guesses;
%             str = sprintf(formatSpec,A1);
            title({file_name, "Convolution Model Results"})
            hold off

            subplot(2,1,2);
            plot(Time,RH,"Color",obj.colors{end},"LineWidth",2,"DisplayName","Relative Humidity [%]")
            xlabel('Time (s)',"Interpreter","latex")
            ylabel('RH [\%]',"Interpreter","latex")
            lgd2 = legend();
            lgd2.Location = 'southeast';
         end 
        end
        
        function Predict(obj)
            
            models_used = Check_Model(obj);

            for i=1:obj.file_count
                model_predict_flux = zeros(length(obj.Replicate(i).Imposed_Data.Time) - 1 ,length(models_used));
                model_pred_inverted_flux = zeros(length(obj.Replicate(i).Imposed_Data.Time) - 1,length(models_used));
                
                % New Instance of class
                pass_data = [obj.Replicate(i).Imposed_Data.Time, obj.Replicate(i).Imposed_Data.HeatFlux, obj.Replicate(i).Imposed_Data.Temp, obj.Replicate(i).Imposed_Data.RH];
                obj.Replicate(i).Predicted_Data = Naylor_Replicate_Basic(pass_data, obj.Replicate(i).Data.file_name, obj.Replicate(i).Data.steady_state_index);
                
                try 
                obj.Replicate(i).Predicted_Data.addprop('HeatFlux_Model_Predictions');
                obj.Replicate(i).Predicted_Data.addprop('HeatFlux_Model_Predictions_Inverted');
                obj.Replicate(i).Predicted_Data.addprop('RH_STEP');
                catch 
                     
                end   
                    
                n_original_rep = length(obj.Replicate(i).Imposed_Data.Time);
                
                ss_index = obj.Replicate(i).Imposed_Data.steady_state_index;

                Temp = obj.Replicate(i).Imposed_Data.Temp;
                new_RH_ss = obj.Replicate(i).Imposed_Data.RH(1,1) * ones(ss_index, 1);
                new_RH_step = max(obj.Replicate(i).Imposed_Data.RH) * ones(n_original_rep - ss_index, 1);
                new_RH = [new_RH_ss; new_RH_step]; % concatenate vertically 
                
                obj.Replicate(i).Predicted_Data.RH_STEP = new_RH;
                
                RHdiff = diff(new_RH);
                original_HeatFlux = obj.Replicate(i).Imposed_Data.HeatFlux;  
                HF0 = mean(original_HeatFlux(1:ss_index));
                
                for ii=1:length(models_used)
                    model_ob = obj.Replicate(i).(models_used{ii});
                    params = model_ob.Best_Convolution_Model_Params;
                    predicted_flux =  model_ob.naylor_curve(params, Temp, ss_index, HF0, RHdiff ); 
                    model_predict_flux(:, ii) = predicted_flux;
                   
                    inverted_pred_flux = Naylor_Material.inverter(predicted_flux, [], 1);  % dont pass baseline, use index of 1
                    model_pred_inverted_flux(:, ii) = inverted_pred_flux;
                    
                end 
                obj.Replicate(i).Predicted_Data.HeatFlux_Model_Predictions = array2table(model_predict_flux,'VariableNames', models_used);
                obj.Replicate(i).Predicted_Data.HeatFlux_Model_Predictions_Inverted = array2table(model_pred_inverted_flux, 'VariableNames', models_used);
                                 
            end 
            
        end 
        
        function Plot_Prediction(obj)
            models_used = Check_Model(obj);
            
            for i=1:obj.file_count
                
                
                
            end 
        end 
        

        
        

    end
end

