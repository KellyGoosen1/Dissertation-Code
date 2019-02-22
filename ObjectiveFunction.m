classdef ObjectiveFunction
% ObjectiveFunction: Objective function used in ABM calibration
% experiments, constructed using the method of simulated moments, as
% presented by Winker et al. (2007). 
%
% The class has attributes:
%
% Name:     Description:                                              
% W         Weight matrix for method of simulated moments              
% R         Measured time series
% meanR     Mean estimate based on R
% stdR      Standard deviation estimate based on R
% kurtR     Kurtosis estimate based on R
% hurstR    Generalized Hurst exponent estimate based on R
%
% In our implementation of the method of simulated moments, we
% consider the mean, standard deviation, kurtosis and generalized Hurst 
% exponent of the measured and simulated series, as well as the 
% Kolmogorov-Smirnov test comparing their empirical CDFs.
%
% When creating a new objective function object (for a new set of
% Thomson-Reuters tick data), the class constructor method should be
% called with the required attributes. This should only be done once,
% however, and this object should be saved to a .mat file for future
% use when calibrating as to avoid the lengthy processes of loading 
% data. It is recommended that ObjectiveFunction Objects be saved with
% the following file name format:
%
% <RIC>_<start date>_<end date>_Bars_<T>.mat
%
% where T is the same as the argument provided to the constructor. 
% A valid example of this format given as follows:
%
% AGLJ_01Nov2013_05Nov2013_Bars_2300.mat
%
% For method details:
% See also ObjectiveFunction/ObjectiveFunction, weightMatrix,
% minuteBars, quoteRepeat.
%
% References:
% 1. Winker P, Gilli M, Jeleskovic V (2007) An Objective Function for 
% Simulation Based Inference on Exchange Rate Data. SSRN:964131
%
% 2. Di Matteo T, Aste T, Dacorogna MM (2003) Scaling behaviors in 
% differently developed markets. Physica A 324:183-188
%
% Donovan Platt
% School of Computer Science and Applied Mathematics
% University of the Witwatersrand, Johannesburg, South Africa
 
    %% Class Attributes
    properties
        W       %Weight matrix for method of simulated moments
        R       %Measured time series
        meanR   %Mean estimate based on R
        stdR    %Standard deviation estimate based on R
        kurtR   %Kurtosis estimate based on R
        hurstR  %Generalized Hurst exponent estimate based on R
        
    end
    
    %% Class Methods
    methods
        %% Class Constructor Method
        function f = ObjectiveFunction(fileName, T, b, nBoot)
        % ObjectiveFunction: Class constructor method.
        %
        % f = ObjectiveFunction(fileName, T, b, nBoot) takes in a file
        % name, fileName, identifying a set of Thomson-Reuters tick data
        % to which an agent-based model is to be calibrated, a measured
        % time series length, T, a window length for a moving block 
        % bootstrap, b, and the desired number of bootstrapped samples to 
        % be generated in the bootstrapping process, nBoot. After receiving 
        % input, the method proceeds to: load the specified data set, 
        % call the minuteBars function to convert the loaded data set into 
        % a series of one minute price bars, perform any required conversions 
        % of these price bars, compute required sample statistics for the 
        % measured series and finally call the weightMatrix method to 
        % construct the weight matrix required by the method of simulated 
        % moments. The weight matrix, measured time series and associated 
        % sample statistics are stored in the appropriate attributes and 
        % the constructed ObjectiveFunction object is returned.
        %
        % See also minuteBars, weightMatrix.
           
           % Load Tick Data
           load(fileName);
           
           % Convert Data to One Minute Price Bars (Final Quote Mid Price 
           % for Each Minute of Data)
           prices = f.minuteBars(data, T);
           
           % Compute Log Prices Based on Measured Prices
           f.R = log(prices);
           
           % Compute Required Sample Statistics of Measured Series
           f.meanR = mean(f.R);
           f.stdR = std(f.R);
           f.kurtR = kurtosis(f.R);
           f.hurstR = f.hurstExponent(f.R);
           
           % Calculate Weight Matrix
           f.W = f.weightMatrix(b, nBoot);
           
        end
        
        %% Calculate Weight Matrix
        function W = weightMatrix(f, b, nBoot)
        % weightMatrix: Use a moving block bootstrap to generate the 
        % weight matrix required by the method of simulated moments to 
        % construct an objective function for calibration.
        %
        % W = weightMatrix(f, b, nBoot) takes in an ObjectiveFunction 
        % object, f, containing a measured time series (prices or returns),
        % R, obtained from empirical data (at any frequency) and proceeds 
        % to create nBoot bootstrapped samples from R using a moving 
        % block bootstrap with a block length of b. Thereafter, the 
        % distributions of the estimates of the moments and test 
        % statistic values considered in the objective function are 
        % determined by computing the moments and test statistic values 
        % for the various bootstrapped samples. Finally, the weight 
        % matrix is set to be the inverse of the variance-covariance 
        % matrix of the estimated moments and test statistics, as 
        % determined by the bootstrapped distributions, such that 
        % higher weights are assigned to moments and test statistics 
        % estimated with greater uncertainty and vice-versa. The weight 
        % matrix, W, is then returned.
        
            % Step 1: Apply a Moving Block Bootstrap to the Measured Series 

            % Calculate Length of Measured Time Series
            n = length(f.R);

            % Create Data Structure to Store Bootstrapped Samples
            bSamples = zeros(n, nBoot); 

            % Generate Indices for First Elements of All Blocks used in
            % Bootstrapping
            blockIndices = 1 : n - b + 1;

            % Generate Bootstrapped Samples in Parallel
            parfor i = 1 : nBoot

                % Sample n/b Blocks of Length b with Replacement
                randBlocks = randsample(blockIndices, n/b, true)';

                % In the line above, we represent blocks by the index of 
                % their first element.

                % Create a Vector of Indices Representing the Elements of the
                % Measured Series to be Extracted to Create a New Sample
                sampleIndices = repmat(randBlocks, 1, b)';
                sampleIndices = sampleIndices(:);
                additionVector = repmat(0 : b - 1, 1, n/b)';
                sampleIndices = sampleIndices + additionVector;

                % In the above, we construct a vector of the form [i, i + 1, 
                % i + 2,..., i + b, j, j + 1,..., j + b,...], where i and j 
                % represent the indices of the first elements of two sampled 
                % blocks.

                % Extract the Sampled Blocks from the Measured Series
                bSamples(:, i) = f.R(sampleIndices);

            end

            % Step 2: Calculate Distributions for Each Moment and Test 
            % Statistic

            % Matrix to Store Moment and Test Statistic Distributions
            dist = zeros(nBoot, 5);

            % Mean and Standard Deviation (General Data Spread)
            dist(:, 1) = mean(bSamples)';
            dist(:, 2) = std(bSamples)';
            
            % Kurtosis (Shape)
            dist(:, 3) = kurtosis(bSamples)';

            % KS Test (Emprical CDF Comparison)
            parfor i = 1 : nBoot
                [~, ~, dist(i, 4)] = kstest2(f.R, bSamples(:, i));
            end
            
            % Hurst Exponent (Scaling)
            parfor i = 1 : nBoot
                dist(i, 5) = f.hurstExponent(bSamples(:, i));
            end
            
            % Step 3: Calculate the Weight Matrix for the Method of 
            % Simulated Moments

            % We set the weight matrix to be the inverse of the approximation 
            % of the variance-covariance matrix of the empirically measured 
            % movements, as determined by the bootstrap distribution. This 
            % assigns larger weights to moments and test statistics which 
            % tend to be estimated with greater uncertainty and vice-versa.
                        
            % Generate Weight Matrix
            W = inv(cov(dist));

        end
        
        %% Evaluate Objective Function
        function functionValue = evaluate(f, simR)
        % evaluate: Evaluate the objective function associated with the
        % object, given a measured time series simulated by the ABM to be 
        % calibrated.
        %
        % functionValue = evaluate(f, simR) takes in an ObjectiveFunction 
        % object, f, and a measured time series, simR, and proceeds to 
        % evaluate the objective function associated with f using the method 
        % of simulated moments. The calculated objective function value is 
        % then returned.
            
            % Create Data Structure to Store G
            G = zeros(size(simR, 2), 5);
            
            % Repeat for Number of Replications
            for i = 1 : size(simR, 2)
                
                % Sample Statistics for Simulated Series
                meanSim = mean(simR(:, i));
                stdSim = std(simR(:, i));
                kurtSim = kurtosis(simR(:, i));
                hurstSim = f.hurstExponent(simR(:, i));

                % Kolmogorov-Smirnov Test Statistic for Simulated and Measured
                % Series
                [~, ~, ksStat] = kstest2(f.R, simR(:, i));

                % Calculate G for Method of Simulated Moments
                G(i, :) = ([meanSim; stdSim; kurtSim; ksStat; hurstSim] - [f.meanR; ...
                    f.stdR; f.kurtR; 0; f.hurstR])';

            end
           
            % Reformat G
            if size(G, 1) > 1
                G = mean(G)';
            else
                G = G';
            end
            
            % Use Calculated Weight Matrix and Method of Simulated Moments 
            % to Obtain an Objective Function Value
            functionValue = G' * f.W * G;
                
        end
        
        %% Evaluate Generalized Hurst Exponent
        function H = hurstExponent(f, R)
        % hurstExponent: Generalized Hurst exponent for a given time
        % series.
        %
        % H = hurstExponent(f, R) takes in an ObjectiveFunction object, f,
        % and a measured time series R, and returns the generalized Hurst
        % exponent for R with q = 1 and tau_max = 19, H.
        
            % Parameters for Estimation Process (Defaults from Literature)
            q = 1;
            max_tau = 19;
            
            % Initialize Variables
            k = 0;         
            H  = [];
            L = length(R);
            
            % Iterate tau_max from 5 to max_tau
            for tau_max = 5 : max_tau
                
                % Increment k
                k = k + 1;
                
                % Values of tau_max Considered in Current Run
                x = 1 : tau_max;
                
                % k_q_t for each tau_max
                k_q_t = zeros(tau_max, 1);
                
                % Iterate tau from 1 to tau_max
                for tau = 1 : tau_max
                    
                    % Calculate Numerator and Denominator for k_q_t
                    numerator = R((tau + 1) : tau : L) - R(((tau + 1) : tau : L) - tau);
                    denominator = R(((tau + 1) : tau : (L + tau)) - tau)';
                    
                    % Determine Drift
                    N = length(numerator) + 1;
                    X = 1 : N;
                    Y = denominator;
                    mx = sum(X) / N;
                    SSxx = sum(X .^ 2) - N * mx ^ 2;
                    my = sum(Y) / N;
                    SSxy = sum(X .* Y) - N * mx * my;
                    cc(1) = SSxy / SSxx;
                    cc(2) = my - cc(1) * mx;
                    
                    % Subtract Drift
                    numerator  = numerator - cc(1);
                    denominator  = denominator - cc(1) .* (1 : N) - cc(2);

                    % Calculate k_q_t
                    k_q_t(tau) = mean(abs(numerator) .^q) / mean(abs(denominator) .^q);
                    
                end
                
                % Calculate Hurst Exponent for Current Iteration
                mx = mean(log10(x));
                SSxx = sum(log10(x) .^ 2) - tau_max * mx^2;
                my = mean(log10(k_q_t));
                SSxy = sum(log10(x) .* log10(k_q_t')) - tau_max * mx * my;
                H(k, 1) = SSxy/SSxx;
                
            end
            
            % Determine Mean Hurst Exponent
            H = mean(H) ./ q;
            
        end
        
        %% Convert Tick Data to One Minute Price Bars
        function processedData = minuteBars(f, rawData, T)
        % minuteBars: Convert Thomson-Reuters tick data to one minute 
        % price bars corresponding to the final quote mid price in each 
        % minute.
        %
        % processedData = minuteBars(f, rawData, T) takes in a set of 
        % Thomson-Reuters tick data, rawData, and a desired number of 
        % one-minute price bars, T, and determines a price time series 
        % corresponding to the final quote mid price for each minute of 
        % trading. The mid price for each quote is calculated as the 
        % average of the L1BidPrice and L1AskPrice Columns. Following 
        % the creation of the time series of one minute price bars for 
        % the entirety of the data, the first T price values are returned, 
        % producing a series of T prices. This method also performs all 
        % necessary data cleaning procedures. In the input arguments
        % above, f is simply an ObjectiveFunction object that allows
        % this method to be called, but its attribute values are not
        % considered in the calculations performed in the method.
        
            % Step 1: Data Cleaning
            
            % Removal of erroneous data items, in this case, NaN rows.
            rawData(cellfun(@(x) isnan(x), rawData(:, 6)), :) = [];

            % Step 2: Repeating Quote Data
            
            % Since the data in the L1BidPrice and L1AskPrice columns only
            % represents relative changes, with a zero indicating no change, 
            % we repeat quote information such that each quote contains the 
            % current best bid and ask.
            rawData = f.quoteRepeat(rawData);

            % Step 3: Extracting Required Data
            
            % Considering that we are only interested in the final quote 
            % mid price for each minute, we are required to extract only 
            % quotes from the data set. We also only extract columns 3, 9 
            % and 11 as these are the only columns we require, since we 
            % only need a timestamp and the value for the best bid and ask 
            % for each quote.
            quoteData = rawData(cellfun(@(x) strcmp(x, 'Quote'), rawData(:, 5) ...
                ), [3 9 11]);

            % From the obtained quote data, it is now required that only 
            % quotes during the trading day are considered, in particular, 
            % quotes between 09:00 and 17:00.
            quoteData(cellfun(@(x) (x(1) == '0' && (str2double(x(2)) < 9))    ...
                || (x(1) == '1' && (str2double(x(2)) >= 7)) || (x(1) == '2'), ... 
                quoteData(:, 1)), :) = [];
                        
            % Step 4: Calculate Mid Price
            % We calculate the mid-price of each quote in order to obtain a 
            % series of prices.
            quoteData(:, 2) = num2cell((cell2mat(quoteData(:, 2)) + ... 
                cell2mat(quoteData(:, 3))) / 2);

            % Step 5: Convert Timestamps to Bar Numbers

            % Initialize Bar Counter
            bar = 1;

            % Set First Bar Number in the Series
            quoteData(1, 3) = {bar};

            % Above, we set the value of the third column in the first row 
            % to be the initial bar number, initialized in the previous step. 
            % Because quote mid prices have already been calculated and 
            % placed in column 2, column 3 no longer serves any purpose and 
            % can be used to store the required bar numbers. We could not 
            % use column 1 as we still require the timestamp for row i-1 
            % when determining the bar number for row i.

            % Iterate Until the Bar Numbers of All Quotes Have been Determined
            for i = 2 : length(quoteData)

                % Timestamps for Current and Preceding Quotes
                currentTime = quoteData{i, 1};
                previousTime = quoteData{i - 1, 1};

                % Check if Current and Preceding Quotes Fall Within the Same
                % Minute
                if ~strcmp(currentTime(4 : 5), previousTime(4 : 5))
                    
                    % Determine Minutes Between Quotes
                    minuteDiff = str2double(currentTime(4 : 5)) - ... 
                        str2double(previousTime(4 : 5));
                    
                    % Check if the Current Quote in in the Following Hour
                    if minuteDiff < 0
                        minuteDiff = minuteDiff + 60;
                    end
                    
                    % Increment Bar Number
                    bar = bar + minuteDiff;
            
                end

                % Allocate Session Number to Current Quote
                quoteData(i, 3) = {bar};
                
            end

            % Delete Time Column, No Longer Required
            quoteData(:, 1) = [];

            % Determine Final Mid-Price for Each Minute
            
            % Obtain Indices for Final Quote Mid Prices of All Trading Minutes
            [~, indices] = unique(cell2mat(quoteData(:, 2)), 'last');
            
            % Extract Final Time Stamps
            timeData = quoteData(indices, 2);
            
            % Extract Final Prices
            quoteData = quoteData(indices, 1);     
            
            % Data Structure to Store Final Prices
            processedData = zeros(timeData{end}, 1);
            
            % Set Calculated Mid Prices
            processedData(cell2mat(timeData)) = cell2mat(quoteData);
            
            % Fill in Gaps
            for i = 1 : length(processedData)
                
                % If No Data for Any Minute, Set Price Equal to Previous
                % Minute's Price
                if processedData(i) == 0
                    processedData(i) = processedData(i - 1);
                end
                
            end
            
            % Remove First and Last 10 Minutes of Each Day of Trading
            processedData = processedData([11 : 470, 491 : 950, 971 : 1430,... 
                1451 : 1910, 1931 : 2390]);
            
            % Check if Full Week Requested
            if T == 2301
                processedData = [processedData; processedData(end)];
            end
            
            % Return the Price Time Series Up to The Desired Number of Bars
            processedData = processedData(1 : T) / 100;
            
            % In the above, we divide by 100 in order to convert from cents
            % to Rands.

        end
        
        %% Repeat Quote Information in Tick Data
        function data = quoteRepeat(f, data)
        % quoteRepeat: Repeat L1BidPrice and L1AskPrice values in a set 
        % of Thomson-Reuters tick data such that each quote contains the 
        % current prevailing best bid and ask.
        % 
        % data = quoteRepeat(f, data) takes in a set of Thomson-Reuters 
        % tick data and repeats the values in the L1AskPrice and 
        % L1BidPrice columns to ensure that each quote contains the 
        % current prevailing best bid and ask. The modified data set is 
        % then returned. In the input arguments above, f is simply an 
        % ObjectiveFunction object that allows this method to be called, 
        % but its attribute values are not considered in the calculations 
        % performed in the method.
        % 
        % Since the data in the L1BidPrice and L1AskPrice columns only 
        % represents changes from the last update, with a 0 representing 
        % no change, it is necessary to repeat this data. The values in 
        % the L1BidSize and L1AskSize columns are also updated for 
        % consistency, though this is not essential.
           
            % Step 1: Repeat Quote Data
            % Check if the L1BidPrice at the index after the current index 
            % is zero, if so, set the values of the L1BidSize and L1BidPrice 
            % columns at the next row index to those at the current index. 
            % Similarly for L1AskPrice and L1AskSize

            % Iterate Until All Quotes are Checked
            for i = 1 : length(data(:, 11)) - 1

                % Ask Prices and Sizes
                if data{i + 1, 11} == 0
                    data{i + 1, 11} = data{i, 11};
                    data{i + 1, 12} = data{i, 12};
                end

                % Bid Prices and Sizes
                if data{i + 1, 9} == 0 || strcmp(data{i + 1, 5}, 'Trade') %Trades contain erroneous values in the L1BidPrice column
                    data{i + 1, 9} = data{i, 9};   %Update L1BidPrice
                    data{i + 1, 10} = data{i, 10}; %Update L1BidSize
                end
                
            end
        end
        
    end
    
end