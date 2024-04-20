classdef CssAdjustment
    %CSSADJUSTMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bench_names
        benchmarks
        constrains
        design
        deltas
        apriori_sigmas
        aposteriori_sigmas
        A                   % Design matrix for the observations
        L                   % Obs vector (delta g values)
        P                   % Weight matrix of the observations
        X
        R                   % Residuals
        So                  % A posteriori variance of unit weight
        Ak                  % Design matrix for the stochastic constraints           
        Lk                  % Constraints (absolute gravity values)
        Pk                  % Weight matrix of the constraints
        instrument_index
        adjusted_g
        adjusted_g_sigma
        fa_anomalies
        ba_anomalies
        ondulations
        residuals
        delta_residuals
        network
        lines
        instruments
        instruments_index
        outliers
        delta_info
    end
    
    methods
        function self = CssAdjustment(lines, benchmarks)
            % get a list of the benchmark names
            self.bench_names = sort({benchmarks.name});
            self.benchmarks  = benchmarks;
            self.lines       = lines;
            self.network     = CssNetwork(lines, benchmarks);
            
            % initialize an empty structure to store the residuals
            self.residuals   = struct();
            self.residuals.start_benchmark = [];
            self.residuals.end_benchmark   = [];
            self.residuals.residual        = [];
            self.residuals.instrument      = [];
            self.residuals.line_name       = [];
            self.residuals.delta_names     = [];
        end
        
        function self = GetDesign(self)
            % this function builds the global design matrix based on the
            % invidivual design matrices of each line-instrument pair
            %global_d = [];
            %data     = [];
            
            h = waitbar(0,'Building design matrix...');
            
%             for i = 1:length(self.lines)
%                 d = self.lines(i).design;
%                 self.lines(i) = self.lines(i).UpdateDesign();
%                 if ~all(cellfun(@isequal, d, self.lines(i).design))
%                     disp(['design of ' self.lines(i).line_name ' changed'])
%                 end
%             end
            % estimate the number of equations to initialize matrix to
            % increase the speed of the function
            t = {self.lines.design}'; % get all the design matrices
            s = cellfun('size',vertcat(t{:}),1); % size of each design matrix
            
            % enumerate all instruments
            t = {self.lines.instruments}';
            self.instruments = unique(vertcat(t{:}));
            self.instruments_index = cell(length(self.instruments),1);
            self.delta_residuals   = cell(length(self.instruments),1);
            self.apriori_sigmas    = [];
            
            global_d = zeros(sum(s), size(self.bench_names,2));
            data     = nan(sum(s), 1);
            
            % counter for the s vector (size of design matrices).
            cc = 1;
            
            for i = 1:length(self.lines)
                for j = 1:length(self.lines(i).instruments)
                    
                    % design matrix for this instrument
                    tdesign = self.lines(i).design{j};
                    
                    if and(~isempty(tdesign),~any(isnan(self.lines(i).deltas{j})))
                        % if design is not empty, then there are deltas to
                        % put into the adjustment
                        
                        % d = zeros(size(tdesign,1),size(self.bench_names,2));
                        
                        % all benchmarks in the line
                        ben = self.lines(i).benchmarks;
                        
                        % find the index in the global benchmark matrix
                        index_st = ismember(self.bench_names, {ben.name});
                        
                        % sort the benchmarks of the line in the same order
                        % as in the global matrix
                        [~,c] = sort({ben.name});
                        
                        % assign the values.
                        %d(:,index_st) = tdesign(:,c);
                        
                        % cumsum(s(1:i)) returns the number of rows in the
                        % global design matrix from 1 to i
                        index = sum(s(1:cc))-s(cc)+1;
                        
                        global_d(index:index+size(tdesign,1)-1,index_st) = tdesign(:,c);
                        data(index:index+size(tdesign,1)-1) = self.lines(i).deltas{j};
                        
                        % check that the size of deltas corresponds to the
                        % size of the design matrix (just as a confirmation
                        % that everything worked ok).
                        if size(self.lines(i).deltas{j},1) ~= size(tdesign,1)
                            disp(['problem with ' self.lines(i).line_name])
                        end
                        
                        % fill the information regarding the instruments.
                        % This basically creates a references of which line
                        % of the design matrix belongs to each instrument
                        inst_i = find(ismember(self.instruments, self.lines(i).instruments(j)));
                        self.instruments_index{inst_i} = [self.instruments_index{inst_i}; (index:index+size(tdesign,1)-1)'];
                        % delta_residuals contains a list of all the
                        % residuals for each instrument so that we can
                        % estimate an a priori value for sigma
                        self.delta_residuals{inst_i}   = [self.delta_residuals{inst_i}; self.lines(i).residuals{j}(~isnan(self.lines(i).residuals{j}))];
                        
                        % fill the start and end benchmarks for the
                        % residuals structure. This facilitates the plot of
                        % the residuals after the adjustment.
                        [delta_name_cell, delta_names] = self.lines(i).GetDeltaName(j, 0);
                        
                        self.residuals.start_benchmark = [self.residuals.start_benchmark; CssBenchmark.ReturnBenchmark(self.benchmarks, {delta_names{:,1}})];
                        self.residuals.end_benchmark   = [self.residuals.end_benchmark;   CssBenchmark.ReturnBenchmark(self.benchmarks, {delta_names{:,2}})];
                        self.residuals.instrument      = [self.residuals.instrument; repmat(self.lines(i).instruments(j), size(delta_names,1), 1)];
                        self.residuals.line_name       = [self.residuals.line_name; repmat({self.lines(i).line_name}, size(delta_names,1), 1)];
                        self.residuals.delta_names     = [self.residuals.delta_names; delta_name_cell'];
                        %global_d = [global_d; d];
                        %data    = [data; self.lines(i).deltas{j}];
                    end
                    cc = cc + 1;
                end
                waitbar(i/length(self.lines))
            end
            close(h);
            
            % remove any nan rows from data (and global_d). These rows
            % exist because sometimes the design matrix is not empty and it
            % get counted in s = cellfun('size',vertcat(t{:}),1), but due
            % to a problem in the line, the deltas are NaN.
            inan = ~isnan(data);
            
            self.design = global_d(inan,:);
            self.deltas = data(inan,:);

            % remove those indeces from instruments_index
            for i = 1:length(self.instruments)
                % translate the indeces to a boolean variable with all the
                % rows in the design matrix. This is necessary to have
                % consistent index mapping between self.design and global_d
                f = false(size(global_d,1));
                f(self.instruments_index{i}) = true;
                
                self.instruments_index{i} = f(inan);
                
                % calculate the apriori sigma for this instrument
                self.apriori_sigmas(i) = std(self.delta_residuals{i});
            end
        end
        
        function self = Invert(self, constrains)
            
            self.constrains = constrains;
            
            % remove the floating lines (orphan benchmarks) from the design
            % matrix
            orphans     = self.network.GetOrphans(constrains);
            no_orphans1 = ~ismember(self.bench_names, orphans);
            self.A      = self.design(:,no_orphans1);
            
            % remove any zero cols from the design matrix
            % self.A = self.design(:,~all(self.design == 0));
            self.L = self.deltas;
            
            % asign the apriori sigmas for each instrument
            self.P = zeros(size(self.A,1),1);
            
            for i = 1:length(self.instruments)
                self.P(self.instruments_index{i}) = 1./self.apriori_sigmas(i).^2;
            end
            
            self.P = diag(self.P);
            
            % We preallocate the matrices of the stochastic constraints, which will be stacked in the final adjustment  
            self.Ak = zeros(length(self.constrains),sum(no_orphans1));
            self.Lk = zeros(length(self.constrains),1);
            self.Pk = zeros(length(self.constrains),length(self.constrains));

            % add the constrains
            for i = 1:length(self.constrains)
                
                % find the constrain
                index = ismember(self.bench_names(no_orphans1), self.constrains(i).name);

                self.Ak(i, index) = 1;                                  % Assign a 1 into the position of the design matrix corresponding to the AG constraint
                self.Lk(i) = self.constrains(i).absolute_g;             % Assign the constraint value (AG value)
                self.Pk(i,i) = 1./self.constrains(i).uncertainty^2;     % Assign the constraint weight

            end
            
            [self.X, S, self.So, self.R, r, ~, ~, a_sigmas, adjust_outliers] = self.adjust_lsq();
            
            % self.X = self.A'*self.P*self.A\self.A'*self.P*self.L;

            % delta residuals
            % self.R = self.A*self.X - self.L;
            
            % sort the residuals in the same order as in the benchmark
            % array. Also remove the benchmarks that are orphans!
            bnames = {self.benchmarks.name};
            
            % change no_orphans to match the order of the benchmarks vector
            no_orphans2 = ~ismember({self.benchmarks.name}, orphans);
            
            [~,c] = ismember(bnames(no_orphans2), self.bench_names(no_orphans1));
            
            % adjusted gravity
            self.adjusted_g = nan(length(self.bench_names),1);
            self.adjusted_g(no_orphans2) = self.X(c);
            
            % uncertainty
            self.adjusted_g_sigma = nan(length(self.bench_names),1);
            self.adjusted_g_sigma(no_orphans2) = S(c);
            
            % the delta residuals should be shown as line segments painted
            % in the color scale. Each residual row corresponds to a row
            % in A where we have the start and end becnhmark
            
            self.residuals.residual = nan(size(self.design,1),1);
            self.residuals.residual(r) = self.R;

            % outliers should be ordered in the same way as residuals
            self.outliers = false(size(self.design,1),1);
            self.outliers(r) = adjust_outliers;
            
            % calculate the aposteriori sigmas of each instrument excluding
            % outliers
            for i = 1:length(self.instruments)
                self.aposteriori_sigmas(i) = nanmean(a_sigmas(and(self.outliers(r), self.instruments_index{i}(r))));
            end
            
        end
        
        function [C, S, So, V, r, dof, cst_pass, sigma, index] = adjust_lsq(self)

            limit = 2.5;

            % select the rows that are not zero
            % some rows end up as zeros due to the elimination of orphans
            r = ~all(self.A == 0,2);

            % don not take into account outliers
            Ai  = self.A(r,:);
            Li = self.L(r);
            Pi = self.P(r,r);   
            Po = Pi;

            Ak = self.Ak;
            Lk = self.Lk;
            Pk = self.Pk;

            cst_pass = false;
            iter = 1;
            factor = 1;

            while and(~cst_pass, iter <= 10) 

                % invert for the parameters, stacking Observations and Constraints matrices
                % A solution obtained by stacking Normal Equations is given by:
                % N = (A' P A)^-1   ;   c = (A' P L)            % Example of Normal Eqs for the Observations 
                % Nk = (Ak' Pk Ak)^-1   ;   ck = (Ak' Pk Lk)    % Example of Normal Eqs for the Constraints
                % x = (N + Nk)^-1 * (c + ck)                    % Example of Least Squares solution stacking Normal Eqs

                C = (Ai'*Pi*Ai + Ak'*Pk*Ak)\(Ai'*Pi*Li + Ak'*Pk*Lk);        % Constrained solution

                % use the input A to get the information about outliers too
                V = Li - Ai*C;          % Residual vector of the Observations
                Vk = Lk - Ak*C;         % Residual vector of the Constraints
                
                dof = (size(Ai,1) - size(Ai,2) + length(Lk));       % dof = unknowns-obs+constraints

                So = sqrt((V'*Pi*V + Vk'*Pk*Vk)/dof);

                x = So.^2.*dof;

                % find the a priori sigma for the observations
                factor = factor.*So;
                % find how many sigmas away the outliers are
                s = abs(V./std(V));

                fprintf('%s',['Iteration: ' num2str(iter) ' std residuals: ' sprintf('%.3f', std(V(s <= limit))) ' variance of unit weight: ' sprintf('%.3f', So)])
                                
                % careful! This function returns the opposite value of alpha as on
                % Leick, page 143
                X1 = chi2inv(1-0.05/2,dof);
                X2 = chi2inv(0.05/2,dof);

                if or(x < X2, x > X1)
                    % if it falls in here it's because it didn't pass the Chi2 test
                    cst_pass = false;

                    if So < 1
                        % weights are too pesimistic, just inform the user
                        fprintf('%s\n', char(hex2dec('25B2')));
                    else
                        % weights are too optimistic, just inform the user
                        fprintf('%s\n', char(hex2dec('25BC')));
                    end

                    % reweigh by Mike's method of equal weight until 2 sigma
                    f = ones(size(V));

                    f(s > limit) = 1./(10.^(limit - s(s > limit)));
                    % do not allow sigmas > 100 m, which is basicaly not putting
                    % the observation in. Otherwise, due to a model problem
                    % (missing jump, etc) you end up with very unstable inversions
                    f(f > 100) = 100;

                    Pi = Po.*diag(1./((factor.*f).^2));     % Note that the weight of the Stochastic Condition Eqs is not scaled by the new So (Constraints are not observations, therefore, should not be scaled)
                else
                    cst_pass = true;
                end

                iter = iter + 1;

            end

            fprintf('\nPassed Chi square test. Done adjusting.\n')
            
            %%%%%%%%%%%% statistics %%%%%%%%%
            % sigmas for the adjusted parameters!
            S = diag(So*sqrt(inv(Ai'*Pi*Ai + Ak'*Pk*Ak)));
            % a posteriori sigmas
            sigma = diag(1./sqrt(Pi));

            % mark observations with residuals > limit
            index = true(size(V));
            index(s > limit) = false;

        end
        
        function self = GetAnomalies(self)
            % obtain the value of N for each benchmark
            self.write_egm2008_input_file(self);
            
            if isunix()
                system('cd egm2008; ./interp_1min');
            else
                system('cd egm2008 & interp_1min.exe');
            end
            
            % read the output
            fid = fopen('egm2008/OUTPUT.DAT');
            self.ondulations = fscanf(fid,' %f %f %f ',[3 inf])';
            self.ondulations = self.ondulations(:,3);
            
            % calculate the anomalies
            [self.fa_anomalies, self.ba_anomalies] = self.combined_gravity_anomaly(self.adjusted_g, [self.benchmarks.lat]'*pi/180, [self.benchmarks.height]' - self.ondulations);
        end
    end
    
    methods(Static)
        
        function write_egm2008_input_file(self)
            fid = fopen('egm2008/INPUT.DAT','w');
            fprintf(fid,'%6.8f %12.8f\n', [[self.benchmarks.lat]; [self.benchmarks.lon]]);
            fclose(fid);
        end
        
        function [delta_ga, delta_gb] = combined_gravity_anomaly(g_meas_v,lat_v,h_v)

            % Hofmann-Wellenhof & Moritz -- Physical Geodesy p. 136
            % from Kevin' code

            for i = 1:numel(g_meas_v)
                g_meas = g_meas_v(i);
                lat = lat_v(i);
                h = h_v(i);

                if g_meas < 10
                    g_meas = g_meas*10^5;
                end

                % delta_gfa = free_air_reduction(g_meas,lat,h);
                % delta_gb = bouguer_gravity(h);
                % 
                % gb = g_meas - delta_gb + delta_gfa;

                % calculate normal gravity relative to wgs84 ellipsoid

                g0 = CssAdjustment.wgs84gravity(lat);

                gb = g_meas + 0.1967*h;

                % use either the first order only term or 1st/2nd Order term
                ga_correction = -0.3086*h;

                delta_ga(i) = g_meas - ga_correction - g0;
                delta_gb(i) = gb - g0; % mGal
            end
            delta_ga = delta_ga';
            delta_gb = delta_gb'; 
        end
        
        function g = wgs84gravity(lat)

            % lat must be in radians!

            if nargin ~= 1
                error('Usgage Error: wgs84gravity(lat) takes 1 inputs')
            end

            % Constants for WGS84 Ellipsoid -- Hofmann-Wellenhof & Moritz, Physical
            % Geodesy, Springer 2005
            GM = 3986004.418 * 10^8;
            J2 = 0.00108262998905;
            a = 6378137;
            omega = 7292115 * 10^-11;

            % g0 = normalgravity(lat,GM,J2,a,omega);
            % 
            % % IAG GRS 80 --> WGS84
            % g0_iag = 9.7803267714*((1+0.00193185138639*sin(lat)^2)/(sqrt(1-0.00669437999013*sin(lat)^2)));
            % g0_grs80 = 9.7803267715*(1+0.0052790414*sin(lat)^2 + 0.0000232718*sin(lat)^4 ...
            %     + 0.0000001262*sin(lat)^6 + 0.0000000007*sin(lat)^8);

            g = nan(size(lat));

            for i = 1:numel(lat)
                % Somigliana Formula
                g_equator = 9.7803253359;
                g_pole    = 9.8321849378;
                b         = 6356752.3142;
                e         = 8.1819190842622*10^-2;
                k = (b*g_pole - a*g_equator)/(a*g_equator);

                g0_somigliana = g_equator*((1+k*sin(lat(i))^2)/(sqrt(1-e^2*sin(lat(i))^2)));
                g(i) = g0_somigliana*10^5;
            end
        end
    end
end

%                         for k = 1:size(tdesign,1)
%                             % start benchmarks
%                             st = GravityLines(i).benchmarks(tdesign(k,:) ==  1);
% 
%                             % end benchmarks 
%                             en = GravityLines(i).benchmarks(tdesign(k,:) == -1);   
% 
%                             % find the locations in the global names array
%                             index_st = ismember(self.names, {st.name});
%                             index_en = ismember(self.names, {en.name});
% 
%                             % assign values to design matrix
%                             d(k,index_st) =  1;
%                             d(k,index_en) = -1;
%                         end