classdef Convolution_Model < handle 
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Model_Types = {'three_exponentials', 'rational_22', 'rational_33', 'rational_43', 'rational_55', 'poly4'}
        Model_Choice_Default = 'three_exponentials'
        Model_Choice
        RMSE
        Best_Convolution_Model_Params
        Model_Fit
        Initial_Guesses = 30
        n_params
        Naylor_Function
    end
    
    methods(Static)
        function m = model_list()
            m = {'three_exponentials', 'rational_22', 'rational_33', 'rational_43', 'rational_55', 'poly4'};
        end 
        
        function TE = triple_exp(params, Time)
        ff = params(1);     % if ff=1, only 1 exotherm timeconstant used in model
        tau1 = params(2);
        tau2 = params(3);  %  2nd time constant if needed
        K = params(4);
        tauReg = params(5);   

        fs= 1-ff;    % fraction of 2nd time constant if needed

        TE = 1*(ff*(1-exp(-Time/tau1)) + fs*(1-exp(-Time/tau2))).*K.*exp(-Time/tauReg);

        end
    end 
    
    methods
        function obj = Convolution_Model(n_guesses, fit_type, Time)
            % Construct an instance of this class
            if isempty(fit_type)
                obj.Model_Choice_Default = 'three_exponentials';
                fit_type = obj.Model_Choice_Default;
            end 
            obj.Initial_Guesses = n_guesses;
            
           
            switch fit_type
            case 'rational_22'
                obj.Naylor_Function =  @(p)(p(1).*Time.^2 + p(2).*Time + p(3)) ./ (Time.^2 + p(4).*Time.^1 + p(5));
                obj.n_params =  5;
                obj.Model_Choice = 'Rational 2,2 Polynomial';
            case 'rational_33'
                obj.Naylor_Function =  @(p)(p(1).*Time.^3 + p(2).*Time.^2 + p(3).*Time + p(4)) ./ (Time.^3 + p(5).*Time.^2 + p(6).*Time.^1 + p(7));
                obj.n_params =  7;
                obj.Model_Choice = 'Rational 3,3 Polynomial';
            case 'rational_43'
                obj.Naylor_Function = @(p)(p(1).*Time.^3 + p(2).*Time.^2 + p(3).*Time + p(4)) ./ (Time.^4 + p(5).*Time.^3 + p(6).*Time.^2 + p(7).*Time + p(8));
                obj.n_params =  8;
                obj.Model_Choice = 'Rational 4,3 Polynomial';
            case 'rational_55'
                obj.Naylor_Function = @(p) (p(1).*Time.^5 + p(2).*Time.^4 + p(3).*Time.^3 + p(4).*Time.^2 + p(5).*Time + p(6)) ./...
                       (Time.^5 + p(7).*Time.^4 + p(8).*Time.^3 + p(9).*Time.^2 + p(10).*Time + p(11));
               obj.n_params =  11;
               obj.Model_Choice = 'Rational 5,5 Polynomial';
            case 'poly4'
                obj.Naylor_Function = @(p) p(1).*Time^.4 + p(2).*Time.^3 + p(3).*Time.^2 + p(4).*Time + p(5);
                obj.n_params =  5;
                obj.Model_Choice = '4th Degree Polynomial';
            case 'three_exponentials'
                obj.Naylor_Function = @(p)Convolution_Model.triple_exp(p,Time);  % simply function down to just 1 variable for the model parameters
                obj.n_params =  5;
                obj.Model_Choice = 'Three Exponentials';
            end 
            
            
        end
        
        function Stochastic_Solver(obj, Material)

            HeatFlux = Material.HeatFlux;
            Temp = Material.Temp;
            RH = Material.RH;
            N = obj.Initial_Guesses;
            ss_idx = Material.steady_state_index;
            RHdiff = diff(RH);
            HF0 = mean(HeatFlux(1:ss_idx));
            lower_bounds = 1e-16;
            upper_bounds = 1e3;

%             if isempty(obj.Initial_Guesses)
%                 N = 30; % Stochastic Sample Number
%             end

            X0 = (upper_bounds - lower_bounds) .* rand(N,obj.n_params) + lower_bounds ;
            opt_array = zeros(N, obj.n_params);

            fmin_options = optimoptions('fmincon','Algorithm','interior-point', 'FiniteDifferenceType', 'central' ,'Display',"off");

%             test_result = optimize_naylor(obj, X0(1, :), Temp, ss_idx, HF0, RHdiff, HeatFlux)
            naylor_trunc = @(p)optimize_naylor(obj, p, Temp, ss_idx, HF0, RHdiff, HeatFlux);  % returns the non-linear least squares as only a function of p, the parameters
    
            lower_lim = 1e-6*ones(1, obj.n_params);
            parfor j=1:N
                opt_array(j,:) = fmincon(@(p)naylor_trunc(p), X0(j,:), [],[],[],[], lower_lim,[],[], fmin_options);
%                         opt_array(j,:) = fmincon(@(x)Naylor_Material.optimize_naylor(x,Time,Temp,ss_idx, HF0, RHdiff, HeatFlux), X0(j,:), [],[],[],[],[0, 1e-6, 1e-6, 0, 1e-6],[],[],fmin_options);

            end

            opt_array_unique = unique(opt_array,'rows');

            error = zeros(N,1);
            parfor j=1:N
                error(j,1) = naylor_trunc(opt_array_unique(j,:));
            end

            [~, min_idx] = min(error,[],1);
            x_best = opt_array_unique(min_idx,:);

            global_best_fit = naylor_curve(obj, x_best, Temp, ss_idx, HF0, RHdiff);

            obj.Initial_Guesses = N;
            obj.Model_Fit = global_best_fit;
            obj.Best_Convolution_Model_Params = x_best;
            obj.RMSE = sum((HeatFlux(1:end-1) - global_best_fit).^2);

        end

        function HFmodel = naylor_curve(obj, params, Temp, ss_idx, HF0, RHdiff)

            Naylor_Result = obj.Naylor_Function(params);
            tester = conv(RHdiff,Naylor_Result,'full');  % convolve 2 vectors together: output length is length(RHdiff)+length(Naylor)-1
            nn=length(tester);
    %         T_end = Time(end);
    %         convTime = 0:ceil(T_end/nn):T_end;

    %         simple_conv_time = 2 * convTime(1: ceil(nn/2));  % sampling every 5 sec -> every 10 seconds  (cuts of end domain where not needed)
            simple_conv_resp = tester(1:ceil(nn/2));

            initTemp_2 = mean(Temp(1:ss_idx));
            TempEffect_2 = HF0*initTemp_2./Temp;
            % disp([length(TempEffect_2), length(simple_conv_resp)])
            HFmodel = TempEffect_2(1: end -1) - simple_conv_resp;  % Naylor Curve

        end
    
        function opt = optimize_naylor(obj, params, Temp, prior_idx, HF0, RHdiff, Flux)
            HFmodel = naylor_curve(obj, params, Temp, prior_idx, HF0, RHdiff);
            opt = sum((HFmodel - Flux(1:end-1)).^2);  % Naylor Curve
        end 
        
        
        function total_sum = integrate_flux(obj)
        end 
    end 
        
end

