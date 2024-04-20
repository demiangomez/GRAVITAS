classdef CssGravityLine
    % CLASS CssGravityLine
    %   handles the gravity line operations such as delta calculation,
    %   instrument drifts, importing from Kevin's g2Gravity structure, etc.
    %   Usage: you either pass a string with the location of the line
    %   structure information or you pass a start and end city with a
    %   directory name.
    
    properties
        line_name
        line_filename
        benchmarks
        observations
        directions
        default_dir % for lines with only one direction, save which one if default.
        start_benchmark
        end_benchmark
        instruments
        drifts
        deltas
        residuals
        status   % saves the benchmark-instrument status (discarded obs). Rows = instruments; cols = benchmarks
        design   % cell array containing the delta-benchmark-instrument design matrix. In each cell one instrument as follows: Rows = obs; cols = benchmarks
        comments % self explanatory field
    end
    
    methods
        function self = CssGravityLine(varargin)
            
            % initialize the children objects
            self.benchmarks      = [];
            self.observations    = [];
            self.directions      = [];
            self.start_benchmark = [];
            self.end_benchmark   = [];
            self.drifts          = [];
            self.status          = [];
            self.design          = [];
            self.default_dir     = [];
            self.comments        = '';
            
            switch nargin
                case 1
                    % calling the class to load a structure file
                    line_file = varargin{1};
                    
                    % class constructor
                    if ~isempty(line_file)
                        % load the data from a line structure
                        S = load(line_file);

                        self = S.line_struct;
                        [~, self.line_filename, ~] = fileparts(varargin{1});
                        self.line_filename = [self.line_filename '.mat'];
                    end
                case 3
                    % calling the class to load a g2Gravity file line and make
                    % the corresponding structure of it
                    path  = varargin{1};
                    scity = varargin{2};
                    ecity = varargin{3};

                    self = load_g2grav_line(self, path, scity, ecity);
            end
        end
        
        function self = SaveGravityLine(self, folder, varargin)
            
            if nargin == 3
                overwrite = varargin{1};
            else
                overwrite = false;
            end
            % save the gravity structure to a mat file
            
            % populate the deltas and drifts
            self = GetDeltas(self);
            
            line_struct = self;
            
            % save to folder
            if isempty(self.line_filename)
                % if the line was not opened from a file, create a file
                % name with the current information
                % DDG: this prevents the line to change name if there is a
                % modification to the start and end points or to the
                % observation timestamps
                line_date = median([self.observations.timestamp]);
                line_yr   = sprintf('%04d',year(line_date));
                line_dy   = sprintf('%02d',day(line_date));
                line_mo   = sprintf('%02d',month(line_date));

                self.line_filename = [self.line_name '-' line_yr '-' line_mo '-' line_dy '.mat'];
            end
            
            if ~overwrite
                if exist(fullfile(folder, self.line_filename),'file')
                    choice = questdlg(['The line file ' fullfile(folder, self.line_filename) ' already exists. Overwrite?'],'','Yes','No','No');

                    if strcmp(choice, 'No')
                        return
                    end
                end
            end
            save(fullfile(folder, self.line_filename), 'line_struct');
        end
        
        function self = AddObservation(self, observation)
            % function to insert an observation given an instrument,
            % direction and all relevant information about the observation
            
            % first, put the observation in (as a copy of the incomming)
            % later, if benchmark doesn't exist, make a handle to this copy
            % of CssObservation 
            self.observations = [self.observations; copy(observation)];
            
            if ~CssBenchmark.exists(self.benchmarks, observation.benchmark.name)
                % benchmark doesn' exist in line, create it!
                % DDG: sep 28 1027 -> If the benchmark doesn't exist but
                % there is more then one instrument, then we are adding a
                % benchmark that is not present in any of the other
                % instruments. Observations to this benchmark should be added  on all other instr
                if length(self.instruments) > 1
                    h = warndlg('This line has more than one instrument and the selected benchmark is not part of the other instrument(s)!');
                    waitfor(h)
                    self.observations(end) = [];
                    return
                else
                    self.benchmarks = [self.benchmarks; CssBenchmark.ReturnBenchmark([self.observations.benchmark], observation.benchmark.name)];
                end
            end
            
            % if we are creating an observation in a line with both
            % directions, force the creation of an observation in the
            % reverse line to match the added observation
            if and(length(self.directions) == 2, observation.direction == CssDirections.reverse)
                
                % get the observations for the flipped direction
                obs = self.GetObservationsByInstrumentDirection(observation.instrument, Flip(observation.direction), true);
                
                if ~ismember(observation.benchmark, [obs.benchmark])
                    % make the time stamp of the obs that of the flipped
                    % direction - one hour (therefore goes to the top).
                    % and change the direction!
                    observation = observation.SwitchDirection();
                    observation.timestamp = min([obs.timestamp]) - hours(1);
                    observation.epoch = cal2jd(year(observation.timestamp),month(observation.timestamp),day(observation.timestamp) + hour(observation.timestamp)/24 + minute(observation.timestamp)/1440);
                    self = self.AddObservation(observation);
                    %h = warndlg('Attention! A forward observation with the same timestamp as the reverse observation was created for this benchmark. Please remember to modify the observation accordingly!');
                    %waitfor(h);
                end
                
            elseif and(length(self.directions) == 1, and(self.default_dir == CssDirections.reverse, observation.direction == CssDirections.forward))
                % a line that was by default a reverse line but now we are 
                % adding a forward direction. We need to copy all
                % observation to the forward direction
                response = questdlg('You are about to add an observation to a reverse-only line, thus transforming this line into a two-way line. To do this, the observations in the reverse line have to exist in the forward. Do you want to copy and flip the current reverse observations to the forward direction? Or do you want to transfer all observations to the forward direction? (In this case, the timestamps will not be modified)','About to convert line', 'Copy and Flip', 'Transfer forward','Cancel', 'Cancel');
                
                if or(strcmp(response,'Copy and Flip'), strcmp(response,'Transfer forward'))
                    transfer = strcmp(response,'Transfer forward');
                    % user selected to copy and flip!
                    obs = self.GetObservationsByInstrumentDirection(observation.instrument, self.default_dir, true);
                    
                    for i = 1:length(obs)
                        % CAREFUL! do operation over a copy of the
                        % observation to leave the originals intact.
                        if transfer
                            obs(i) = obs(i).SwitchDirection();
                        else
                            nobs = copy(obs(i));
                            % flip by making the timestamp == to the min - 1
                            % hour
                            nobs = nobs.SwitchDirection();
                            nobs.timestamp = min([obs.timestamp]) - hours(i);
                            nobs.epoch = cal2jd(year(nobs.timestamp),month(nobs.timestamp),day(nobs.timestamp) + hour(nobs.timestamp)/24 + minute(nobs.timestamp)/1440);
                            % add observation
                            self.observations = [self.observations; copy(nobs)];
                        end
                    end
                    h = warndlg('Copy and flip done. Please remember to modify the observations accordingly!');
                    waitfor(h);
                else
                    % remove from line
                    benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks, observation.benchmark.name);
                    self.benchmarks(self.benchmarks == benchmark) = [];
                    
                    self.observations(end) = [];
                    return
                end
            end
            
            self = self.sort_benchmarks();
        end
        
        function self = DeleteObservation(self, observations)
            % delete the observation passed as argument
            self.observations(ismember(self.observations,observations)) = [];
            % resort observations and benchmarks
            self = self.sort_benchmarks();
        end
        
        function self = GetDeltas(self)
           
            if length(self.directions) == 2
                % this line has two directions 
            
                for i = 1:length(self.instruments)
                    
                    ndeltas = size(self.design{i},1);
                    % compute the drift rate and deltas for each instrument
                    
                    % get the forward and reverse lines from the
                    % observation array for this instrument
                    fwd = self.GetObservationsByDirection(self.GetObservationsByInstrument(self.observations, self.instruments{i}), CssDirections.forward);
                    rev = self.GetObservationsByDirection(self.GetObservationsByInstrument(self.observations, self.instruments{i}), CssDirections.reverse);

                    % if number of fwd and rev observations is not equal
                    % then this is a line under construction. Do not
                    % calculate the deltas yet
                    if length(fwd) ~= length(rev)
                        self.deltas{i,1}    = NaN;
                        self.residuals{i,1} = NaN;
                        self.drifts(i,1)    = NaN;
                        continue
                    end
                    
                    % if branch to handle instruments with all benchmarks
                    % disabled
                    if ndeltas ~= 0
                        % identify lines with 2 directions but with multiple
                        % observations using same gravimeter
                        if ndeltas ~= length(diff([fwd.reduced_g]'))
                            % separate the multiple passes in indivual columns
                            clear dfwd tfwd tfwd drev trev;
                            k      = 1;
                            c      = 1;
                            tnames = {};
                            for j = 1:length(fwd)
                                if ~ismember(fwd(j).benchmark.name, tnames)
                                    dfwd(k,c) = fwd(j).reduced_g;
                                    tfwd(k,c) = fwd(j).timestamp;
                                    tnames(k) = {fwd(j).benchmark.name};
                                    k         = k + 1;
                                else
                                    % we are back at the beginning of the line
                                    % reset variables
                                    tnames = {};
                                    c = c + 1;
                                    dfwd(1,c) = fwd(j).reduced_g;
                                    tfwd(1,c) = fwd(j).timestamp;
                                    k = 2;
                                end
                            end

                            % do the same for the reverse
                            k      = 1;
                            c      = 1;
                            tnames = {};
                            for j = 1:length(rev)
                                if ~ismember(rev(j).benchmark.name, tnames)
                                    drev(k,c) = rev(j).reduced_g;
                                    trev(k,c) = rev(j).timestamp;
                                    tnames(k) = {rev(j).benchmark.name};
                                    k         = k + 1;
                                else
                                    % we are back at the beginning of the line
                                    % reset variables
                                    tnames = {};
                                    c = c + 1;
                                    drev(1,c) = rev(j).reduced_g;
                                    trev(1,c) = rev(j).timestamp;
                                    k = 2;
                                end
                            end

                            % build the difference vectors
                            dfwd = diff(dfwd);
                            tfwd = hours(diff(tfwd));
                            % for the reverse line
                            drev = diff(drev);
                            trev = hours(diff(trev));

                            A = [-tfwd(:)         repmat(diag(ones(ndeltas,1)) ,size(dfwd,2),1);
                                 -trev(:) repmat(flipud(-diag(ones(ndeltas,1))),size(trev,2),1)];
                        else
                            dfwd = diff([fwd.reduced_g]');
                            drev = diff([rev.reduced_g]');

                            tfwd = hours(diff([fwd.timestamp]'));
                            trev = hours(diff([rev.timestamp]'));

                            % flipud in diag is used to match the order or the ida matrix
                            % therefore, results are expressed in ida order
                            A = [-tfwd         diag(ones(size(tfwd,1),1)); 
                                 -trev flipud(-diag(ones(size(trev,1),1)))];
                        end

                        x = A'*A\A'*[dfwd(:); drev(:)];

                        % save the results in the structure
                        self.drifts(i,1) = x(1);

                        % save the deltas and residuals
                        self.deltas{i,1} = x(2:end);

                        v = [dfwd(:); drev(:)] - A*x;
                        % save the residuals of the forward direction only
                        self.residuals{i,1} = v(1:ndeltas);
                    else
                        % there are no activated benchmarks for this
                        % instrument
                        self.deltas{i,1}    = NaN;
                        self.residuals{i,1} = NaN;
                        self.drifts(i,1)    = NaN;
                    end
                end
            else
                % line with only one direction
                % example: given the following observations:
                % A -> B -> C -> A -> D -> E ...
                % calculate the drift and the deltas between the benchmarks
                % of the line, which is assumed to be A,B,C,D,E ...
                
                for i = 1:length(self.instruments)
                    
                    % get only the "active" observations
                    obs = self.GetObservationsByInstrument(self.observations, self.instruments{i});
                    
                    % find the unique points that define the delta g
                    k      = 1;
                    h      = 1;
                    tnames = {};
                    uobs   = [];
                    robs   = [];
                    for j = 1:length(obs)
                        % verify that the benchmark is NOT a member if
                        % tnames 
                        if ~ismember(obs(j).benchmark.name, tnames)
                            uobs      = [uobs; obs(j)];
                            tnames(k) = {obs(j).benchmark.name};
                            k         = k + 1;
                        else
                            % repeated observations
                            % find the location of the first occurence of
                            % this repeated observation
                            bench     = [uobs.benchmark];
                            index     = strcmp({bench.name}, obs(j).benchmark.name);
                            robs      = [robs; uobs(index)];
                            robs      = [robs;      obs(j)];
                            h         = h + 2;
                        end
                    end
                    
                    if ~isempty(robs) 
                        % get the time and g difference for these observations
                        d_uobs = diff([uobs.reduced_g]');
                        t_uobs = hours(diff([uobs.timestamp]'));

                        % take difference of reobserved benchmarks
                        d = [robs.reduced_g]';
                        t = [robs.timestamp]';
                        d_robs = d(2:2:end) - d(1:2:end);
                        t_robs = hours(t(2:2:end) - t(1:2:end));

                        % # of rows of the design matrix = deltas we should
                        % expect from this adjustment
                        A = [-t_uobs diag(ones (size(t_uobs,1),1)); 
                             -t_robs zeros(size(t_robs,1),size(t_uobs,1))];

                        x = A'*A\A'*[d_uobs; d_robs];

                        % save the results in the structure
                        self.drifts(i,1) = x(1);

                        % save the deltas and residuals
                        self.deltas{i,1} = x(2:end);

                        v = [d_uobs; d_robs] - A*x;
                        % save the residuals of the forward direction only
                        self.residuals{i,1} = v(1:length(d_uobs));
                    else
                        % missing re-observed information, can't calculate
                        % anything
                        self.deltas{i,1}    = NaN;
                        self.residuals{i,1} = NaN;
                        self.drifts(i,1)    = NaN;
                    end
                end
            end
        end
        
        function out = GetObservationsByBenchmark(self, benchmark, varargin)
            % this function returns the observations from a particular
            % instruments, specified in pinstruments. This function will by
            % default return only the active observations (those not 
            % filtered). This behavior can be overriden by setting an
            % optional argument (varargin) that allows to get ALL the
            % observations, including those flagged as inactive.

            if ~isempty(varargin)
                all = varargin{1};
            else
                all = false;
            end
            
            if ~isempty(self.observations)
                ben   = [self.observations.benchmark];
                if isa(benchmark,'CssBenchmark')
                    array = self.observations(strcmp({ben.name}, benchmark.name));
                elseif isa(benchmark,'char')
                    array = self.observations(strcmp({ben.name}, benchmark));
                else
                    error('invalid input')
                end
                
                % is the user requesting all observations?
                if ~all
                    % find the status of these observations
                    out = array(self.GetObservationStatus(array));
                else
                    out = array;
                end
            else
                out = [];
            end
        end
        
        function out = GetObservationsByInstrument(self, pobservations, pinstrument, varargin)
            % this function returns the observations from a particular
            % instruments, specified in pinstruments. This function will by
            % default return only the active observations (those not 
            % filtered). This behavior can be overriden by setting an
            % optional argument (varargin) that allows to get ALL the
            % observations, including those flagged as inactive.

            if ~isempty(varargin)
                all = varargin{1};
            else
                all = false;
            end
            
            if ~isempty(pobservations)
                array = pobservations(strcmp({pobservations.instrument}, pinstrument));

                % is the user requesting all observations?
                if ~all
                    % find the status of these observations
                    out = array(self.GetObservationStatus(array));
                else
                    out = array;
                end
            else
                out = [];
            end
            
        end
        
        function out = GetObservationsByDirection(self, pobservations, pdirection, varargin)
            % this function returns the observations from a particular
            % direction, specified in pdirection. This function will by
            % default return only the active observations (those set to 
            % 1 in the status table). This behavior can be overriden by 
            % setting an optional argument (varargin) that allows to get 
            % ALL the observations, including those flagged as inactive.
            
            if ~isempty(varargin)
                all = varargin{1};
            else
                all = false;
            end
            
            if ~isempty(pobservations)
                array = pobservations([pobservations.direction] == pdirection);

                % is the user requesting all observations?
                if ~all
                    % find the status of these observations
                    out = array(self.GetObservationStatus(array));
                else
                    out = array;
                end
            else
                out = [];
            end
        end
        
        function out = GetObservationsByInstrumentDirection(self, pinstrument, pdirection, varargin)
            % uses GetObservationsByDirection and
            % GetObservationsByDirection in a sequence to obtain the
            % requested observations
            
            if ~isempty(varargin)
                all = varargin{1};
            else
                all = false;
            end
            
            obs = self.GetObservationsByInstrument(self.observations, pinstrument, all);
            obs = self.GetObservationsByDirection(obs, pdirection, all);
            
            out = obs;
        end
        
        function out = GetObservationStatus(self, observation)
            % return the status of an observation based on the
            % instrument-benchmark status table
            out = true(size(observation));
            for i = 1:length(observation)
                ins_index = strcmp(self.instruments, observation(i).instrument);

                % find the benchmark in the benchmarks array
                ben_index = self.benchmarks == observation(i).benchmark;
                
                out(i) = self.status(ins_index, ben_index);
            end
        end
        
        function self = UpdateDesign(self)
            % this function updates the design matrix based on the
            % status matrix
            % columns = benchmarks
            % rows    = deltas
            
            self.design = {};
            
            for i = 1:length(self.instruments)
                % get all the observations, not just the active ones
                obs = self.GetObservationsByInstrument(self.observations, self.instruments{i}, true);
                % if a single direction available, then the
                % observations that determine the deltas are those NOT
                % repeated.
                % this code also works for lines with two directions
                % that have been observed multiple times (for the same
                % direction) using the same gravimeter.
                k      = 1;
                tnames = {};
                tobs   = [];
                for j = 1:length(obs)
                    % verify that the benchmark is NOT a member if
                    % tnames
                    if and(~ismember(obs(j).benchmark.name, tnames), ...
                            or(obs(j).direction == CssDirections.forward, length(self.directions) == 1)) 
                        tobs   = [tobs; obs(j)];
                        tnames(k) = {obs(j).benchmark.name};
                        k         = k + 1;
                    end
                end
                % replace the current array
                if ~isempty(tobs)
                    [~, ind] = sort([tobs.timestamp]);
                    obs = tobs(ind);

                    % the the observation count for this instrument
                    obs_count   = self.GetObservationStatus(obs);
                    delta_count = sum(obs_count) - 1;   % the deltas are the active observations - 1
                    delta_order = find(obs_count == 1); % an ordered index array showing the delta members

                    d = zeros([delta_count length(self.benchmarks)]);

                    for j = 1:delta_count
                        % assign the values (1 or -1) to each benchmark pair
                        d(j, delta_order(j  )) = -1;
                        d(j, delta_order(j+1)) =  1;
                    end
                else
                    d = [];
                end
                % save the design matrix for this instrument
                self.design{i,1} = d;
            end
        end
        
        function self = DeactivateObservationPair(self, instrument, benchmark)
            % function used to deactivate a benchmark observation-pair from
            % a particular instrument
            ins_index = strcmp(self.instruments, instrument);
            obs_index = self.benchmarks == benchmark;
            self.status(ins_index, obs_index) = false;
            
            % trigger the update of the line
            self = self.UpdateDesign();
            self = self.GetDeltas();
        end
        
        function self = ActivateObservationPair(self, instrument, benchmark)
            % function used to activate a benchmark observation-pair from
            % a particular instrument
            ins_index = strcmp(self.instruments, instrument);
            obs_index = self.benchmarks == benchmark;
            self.status(ins_index, obs_index) = true;
            
            % trigger the update of the line
            self = self.UpdateDesign();
            self = self.GetDeltas();
        end
        
        function self = UpdateObservation(self, observation, field, newValue, varargin)
            % update a field of the observation object member of this class
            % possible fields: offset, timestamp, reading[123], active.
            
            switch lower(field)
                case 'timestamp'
                    if isempty(varargin)
                        error('this operation requires the instrument calibration')
                        return
                    else
                        observation = observation.UpdateTimeStamp(varargin{1}, newValue);
                        self = self.sort_benchmarks();
                        self = self.GetDeltas();
                    end
                case 'offset'
                    if isempty(varargin)
                        error('this operation requires the instrument calibration')
                        return
                    else
                        observation = observation.UpdateOffset(varargin{1}, newValue);
                        self = self.GetDeltas();
                    end
                case 'reading'
                    if length(varargin) < 2
                        error('this operation requires the instrument calibration and the reading index as arguments')
                        return
                    else
                        observation = observation.UpdateReading(varargin{1}, varargin{2}, newValue);
                        self = self.GetDeltas();
                    end
                case 'active'
                    stat = GetObservationStatus(self, observation);
                    if ~stat
                        self = self.ActivateObservationPair(observation.instrument, observation.benchmark);
                    else
                        self = self.DeactivateObservationPair(observation.instrument, observation.benchmark);
                    end
            end
        end
        
        function self = UpdateBenchmark(self, original_benchmark, replacement_benchmark, instrument_array)
            % Updates a benchmark or replaces one benchmark for another.
            % When an existing benchmark is replaced, it is replaced from ALL
            % observations where it appears.

            if isa(original_benchmark, 'CssBenchmark')
                original_benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks, original_benchmark.name);
            elseif isa(original_benchmark, 'char')
                original_benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks, original_benchmark);
            end

            if ~CssBenchmark.exists(self.benchmarks, replacement_benchmark.name)
                % not in the list! add it
                self.benchmarks = [self.benchmarks; copy(replacement_benchmark)];
                % find it from self.benchmarks
                replacement_benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks, replacement_benchmark.name);
            else
                % this is an existing benchmark, probably a
                % modified version.
                if strcmp(original_benchmark.name, replacement_benchmark.name)
                    % only replace if 
                    % it's the same benchmaker (with different properties)
                    self.benchmarks(self.benchmarks == original_benchmark) = copy(replacement_benchmark);
                end
                % find the benchmark in self.benchmarks
                replacement_benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks, replacement_benchmark.name);
            end

            % DDG: New behavior: leave offset untouched! Offset is a
            % benchmark variable but I allow it to have a different value
            % in each line rather than the value declared in the
            % database.benchmarks list.
            replacement_benchmark.offset = original_benchmark.offset;
            
            % determine if offsets are different, in which case we
            % have to recalculate the deltas
            if replacement_benchmark.offset ~= original_benchmark.offset
                trigger_recalc = true;
            else
                trigger_recalc = false;
            end

            for i = 1:length(self.observations)
                if self.observations(i).benchmark == original_benchmark
                    % find the corresponding instrument
                    instrument = instrument_array(strcmp({instrument_array.name},self.observations(i).instrument));
                    % change the benchmark!
                    self.observations(i).ChangeBenchmark(replacement_benchmark, instrument.calibration);
                end
            end
            % remove original_benchmark from the list of benchmarks
            self.benchmarks(self.benchmarks == original_benchmark) = [];
            % reorganize the benchmarks according to the
            % observation order
            self = self.sort_benchmarks();
            if trigger_recalc
                self = self.GetDeltas();
            end
        end
        
        function self = UpdateTimeStamp(self, instrument, benchmark_name, oldTimeStamp, newTimeStamp)
            % this function updates the time stamp of an observation given
            % the instrument, benchmark name and timestamp
            %
            obs = self.GetObservationsByInstrument(self.observations, instrument.name);
            % find the observation with the "oldTimeStamp" and
            % benchmark_name
            index = [obs.timestamp] == oldTimeStamp;
            
            if and(strcmp(obs(index).benchmark.name, benchmark_name), sum(index) == 1)
                % the benchmark and timestamp agree, update information
                % no need to assign this variable to anything. It's a
                % handle pointing to self.observations.
                obs(index) = obs(index).UpdateTimeStamp(instrument.calibration, newTimeStamp);
                % recalculate deltas and residuals to propagate this change
                self = self.GetDeltas();
            else
                error(['The requested benchmark name ' benchmark_name ' and time stamp did not yield a unique observation'])
            end
        end
        
        function [out, vect] = GetDeltaName(self, instrument, delta)
            % this function returns the name of a delta (X-Y) given the
            % instrument number (as in the instruments list) and delta
            % if delta == 0, then it returns all the delta names
            
            if ischar(instrument)
                % passed the instrument name, turn it into an index
                instrument = strcmp(self.instruments, instrument);
                if ~any(instrument)
                    out     = [];
                    vect    = [];
                    return
                end
            end
            
            tdesign = self.design{instrument};
            out     = [];
            vect    = [];
            if delta == 0
                for i = 1:size(tdesign,1)
                    out = [out, {[self.benchmarks(tdesign(i,:) == 1).name '-' self.benchmarks(tdesign(i,:) == -1).name]}];
                    
                    if nargout > 1
                        vect = [vect; {self.benchmarks(tdesign(i,:) == 1).name} {self.benchmarks(tdesign(i,:) == -1).name}];
                    end
                end
            else
                if or(delta > size(tdesign,1), delta < 0)
                    error('Invalid delta value!')
                else
                    out = {[self.benchmarks(tdesign(delta,:) == 1).name '-' self.benchmarks(tdesign(delta,:) == -1).name]};
                    
                    if nargout > 1
                        vect = [{self.benchmarks(tdesign(delta,:) == 1).name} {self.benchmarks(tdesign(delta,:) == -1).name}];
                    end
                end
            end
        end
        
        function [deltas, residuals] = GetDeltasResiduals(self, instrument_name)
            % find the instrument in the deltas list
            i = find(strcmp(self.instruments,instrument_name));
            if i ~= 0
                deltas    = self.deltas{i};
                residuals = self.residuals{i};
            else
                deltas    = [];
                residuals = [];
            end
        end
        
        function drift = GetDrifts(self, instrument_name)
            % find the instrument in the deltas list
            i = find(strcmp(self.instruments,instrument_name));
            if i ~= 0
                drift = self.drifts(i);
            else
                drift = [];
            end
        end
        
        function out = ToList(self, observations)
            % convert list of observations to a list compatible with a 
            % uitable matlab control
            for i = 1:length(observations)
                
                % find if observations for this benchmark are active
                stat = self.GetObservationStatus(observations(i));
                
                for j = 1:length(observations(i).raw_data)
                    readings(j) = observations(i).raw_data(j);
                end
                
                if length(readings) < 3
                    readings(3) = 0;
                end
                
                out(i,:) = {observations(i).benchmark.name, datestr(observations(i).timestamp,'yyyy-mm-dd HH:MM'), ...
                            sprintf('%.3f', observations(i).benchmark.offset), sprintf('%.3f', readings(1)), sprintf('%.3f', readings(2)), sprintf('%.3f', readings(3)), ...
                            sprintf('%.3f', observations(i).reading), sprintf('%.3f', observations(i).reduced_g), stat};
            end
            
            if ~exist('out', 'var')
                out = [];
            end
        end
        
        function plot(self, instrument, varargin)
            
            if isempty(varargin)
                stdmax = 0.1;
            else
                stdmax = varargin{1};
            end
            
            for i = 1:length(self.instruments)
                if strcmp(self.instruments(i), instrument)
                
                    % plot deltas
                    if ~all(self.status(i,:) == false)
                        cla('reset')
                        
                        d = self.deltas{i};
                        r = self.residuals{i};
                        
                        % get the delta vector of size == benchmarks
                        [d, crossoutd] = self.check_deactivated_obs(d, self.status(i,:));
                        [r, crossoutr] = self.check_deactivated_obs(r, self.status(i,:));
                        
                        dh = stairs([d;d(end)],'-o');
                        hold on
                        
                        % mark the deactivated benchmark (deltas)
                        if any(self.status(i,:) == false)
                            plot(crossoutd(:,1), crossoutd(:,2),'xr','MarkerSize',10)
                        end
                        
                        ylabel('Delta G [mGal]')
                        yyaxis right
                        rh = stairs([r; r(end)],'--o');
                        ylabel('Residuals [mGal]')

                        % plot the limits
                        sh = plot([1 size(d,1)+1]', 3*[stdmax stdmax]','--r');
                        plot([1 size(d,1)+1]',-3*[stdmax stdmax]','--r')
                        
                        % mark the deactivated benchmark (residuals)
                        if any(self.status(i,:) == false)
                            plot(crossoutr(:,1), crossoutr(:,2),'xr','MarkerSize',10)
                        end
                        
                        title(['Deltas and Residuals - Instrument: ' self.instruments{i} ' (drift rate: ' num2str(self.drifts(i)) ')'])

                        set(gca,'XTick',1:length([self.benchmarks]))
                        set(gca,'XTickLabel',{self.benchmarks.name})
                        grid on
                        axis tight
                        xtickangle(25)
                        legend([dh,rh,sh],'Deltas','Residuals','3\sigma limit')
                    else
                        title(['Deltas and Residuals - Instrument: ' self.instruments{i} ' (NO OBSERVATIONS)'])
                        grid on
                    end
                end
            end
        end
        
        function plotComparison(self)
            
            hold on
            c = 1;
            instr = [];
            style = {'-o','--o',':o','-.o','-o','--o',':o','-.o'};
            
            for i = 1:length(self.instruments)
                % plot deltas
                if ~all(self.status(i,:) == false)
                    d = self.deltas{i};
                    
                    % check if there is any deactivated observation
                    [d,crossout] = self.check_deactivated_obs(d, self.status(i,:));
                    
                    h(c) = stairs([d;d(end)],style{i});
                    
                    % mark the deactivated benchmark
                    if any(self.status(i,:) == false)
                        plot(crossout(:,1), crossout(:,2),'xr','MarkerSize',10)
                    end
                            
                    c = c+1;
                    instr = [instr; self.instruments(i)];
                end
            end
            if ~isempty(h)
                legend(h, instr)
            end
            
            set(gca,'XTick',1:length([self.benchmarks]))
            set(gca,'XTickLabel',{self.benchmarks.name})
            xtickangle(25)
                    
            title('Instruments Delta Comparison')
            ylabel('Delta G [mGal]')
            grid on
            axis tight
        end
        
        function plotRaw(self, instrument)
            
            style = {'-o','--o',':o','-.o'};

            for i = 1:length(self.instruments)
                
                if strcmp(self.instruments(i), instrument)
                    cla('reset')
                    tobs = self.GetObservationsByInstrument(self.observations, self.instruments{i}, true);

                    h = [];
                    c = 1;
                    plot_directions = [];
                    
                    for j = 1:length(self.directions)

                        obs = self.GetObservationsByDirection(tobs, self.directions(j), true);

                        if ~isempty(obs)
                            if self.directions(j) ~= self.default_dir
                                obs = flip(obs);
                            end

                            h(c) = plot(0:length(obs)-1, [obs.reduced_g], style{j});
                            c = c +1;
                            plot_directions = [plot_directions; self.directions(j)];
                            
                            hold on
                            % plot deactivated observations
                            if any(self.status(i,:) == false)
                                plot(find(self.status(i,:) == false)-1, [obs(self.status(i,:) == false).reduced_g],'xr','MarkerSize',10)
                            end
                            
                            set(gca,'XTick',0:length(obs)-1)
                            benchs = [obs.benchmark];
                            set(gca,'XTickLabel',{benchs.name})
                    
                        end
                    end

                    if ~all(self.status(i,:) == false)
                        title(['Observations - Instrument: ' self.instruments{i} ' (drift rate: ' num2str(self.drifts(i)) ')'])
                    else
                        title(['Observations - Instrument: ' self.instruments{i} ' (NO OBSERVATIONS)'])
                    end

                    grid on
                    axis tight
                    ylabel('Delta G [mGal]')

                    if ~isempty(plot_directions)
                        legend(h, str(plot_directions))
                    end
                end
            end
        end
        
        function plotLineTime(self, instrument, varargin)
            
            style = {'-o','--o',':o','-.o'};
            
            if ~isempty(varargin)
                direction = varargin{1};
            else
                direction = [CssDirections.forward; CssDirections.reverse];
            end
            
            plot_directions = [];
            
            % plot time in the X axis
            for i = 1:length(self.instruments)
                
                if strcmp(self.instruments(i), instrument)
                    cla('reset')
                    
                    tobs = self.GetObservationsByInstrument(self.observations, self.instruments{i}, true);

                    h = [];

                    for j = 1:length(self.directions)

                        obs = self.GetObservationsByDirection(tobs, self.directions(j), true);
                        
                        if or(~ismember(self.directions(j), direction), isempty(obs))
                            % check if this direction should be ploted or
                            % not
                            continue
                        end
                        
                        plot_directions = [plot_directions; self.directions(j)];
                        
                        % plot deltas
                        if self.directions(j) ~= self.default_dir
                            obs = flip(obs);
                        end

                        % get the distances between points
                        bench = [obs.benchmark];

                        lat = [bench.lat]';
                        lon = [bench.lon]';

                        d = cumsum([0; m_idist(lon(1:end-1),lat(1:end-1),lon(2:end),lat(2:end))])/1000;

                        ts = [obs.timestamp]';
                        h(length(plot_directions)) = plot(ts, d, style{j}); 

                        hold on

                        % plot the velocities
                        if length(ts) >= 2
                            t = text(diff(ts)./2 + ts(1:end-1,1), diff(d)./2 + d(1:end-1,1), sprintfc('\\leftarrow%.0f km/h', abs(diff(d)./hours(diff(ts))) ));
                            set(t,'Rotation',90,'FontSize',7.5);
                        end
                        % plot the benchmark labels
                        %text([obs.timestamp], d, strcat('\leftarrow', {bench.name}), 'FontSize', 8);

                        % plot deactivated observations
                        if any(self.status(i,:) == false)
                            plot([obs(self.status(i,:) == false).timestamp], d(self.status(i,:) == false),'xr','MarkerSize',10)
                        end

                        % put labels with the benchmark names
                        % adjust the y axis size first
                        ylim('auto')
                        % get limits
                        y = ylim;
                        % plot one dummy label to get the extent
                        t = text(obs(1).timestamp, 0, 'MMMM\rightarrow', 'FontSize', 11, 'FontName', 'Courier', 'HorizontalAlignment','right');
                        set(t,'Rotation',90);
                        extent = get(t,'Extent');
                        % reasign the limits to fit the label
                        ylim([extent(2)+extent(2)/2.8 y(2)])
                        % delete the dummy label
                        delete(t);

                        t = text([obs.timestamp]', zeros(size([obs.timestamp]',1),1), strcat({bench.name}, '\rightarrow'), 'FontSize', 11, 'FontName', 'Courier', 'HorizontalAlignment', 'right');
                        set(t,'Rotation',90);
                        % draw a dotted line to the point
                        plot([[obs.timestamp]; [obs.timestamp]], [zeros(size([obs.timestamp]',1),1)'; d'], ':k')
                        set(gca,'TickLength',[0 0])
                    end
                    title(['Observations - Instrument: ' self.instruments{i} ' (drift rate: ' num2str(self.drifts(i)) ')'])

                    ylabel('Distance [km]')
                    xlabel('Observation Time Stamp [UTC]')
                    grid on
                    if ~isempty(plot_directions)
                        legend(h, str(plot_directions))
                    end
                end
            end
        end
        
        function self = sort_benchmarks(self)
            
            if isempty(self.observations)
                % no observations, do nothing!
                return
            end
            
            % sort OBSERVATIONS by instrument/epoch
            [~,i] = sortrows([{self.observations.instrument}; {self.observations.epoch}]', [1 2]);
            self.observations = self.observations(i);
            
            % Determine the directions and instruments available for this 
            % line. Some lines may not have both directions
            self.directions  = unique([self.observations.direction]');
            self.instruments = unique({self.observations.instrument}');
            
            % determine which direction we have
            if length(self.directions) == 2
                self.default_dir = CssDirections.forward;
            else
                if ismember(CssDirections.forward, self.directions)
                    self.default_dir = CssDirections.forward;
                else
                    self.default_dir = CssDirections.reverse;
                end
            end
            
            fObservations = [];
            % determine the end points using one of the instruments
            for i = 1:length(self.instruments)
                % use the instrument that yields the largest number of
                % observations
                obs = self.GetObservationsByInstrumentDirection(self.instruments{i}, self.default_dir, true);
                
                if length(obs) > length(fObservations)
                    fObservations = obs;
                end
            end
            
            %fObservations = self.GetObservationsByInstrument(self.observations, self.instruments{1}, true);
            % filter the desired direction (forward) to determine the
            % start and end benchmarks
            %fObservations = self.GetObservationsByDirection(fObservations, self.default_dir, true);
            % some lines have multiple observations for the same instrument
            % and direction (even when there is 2 directions) -> get the 
            % unique benchmarks
            % this also happens in lines with a single direction
            % the same benchmark name appears
            % multiple times for a single direction (reobservations)
            % thus, without this check the function creates multiple copies
            % of the same benchmark
            k      = 1;
            tnames = {};
            for j = 1:length(fObservations)
                % verify that the benchmark is NOT a member if
                % tnames and that the observation is active
                if ~ismember(fObservations(j).benchmark.name, tnames)
                    tnames(k) = {fObservations(j).benchmark.name};
                    k         = k + 1;
                end
            end

            % DDG: ATTENTION!!
            % when a line is being built (data input) the default direction
            % (forward) might might have a missing observation from a
            % benchmark because an observation was added to the reverse
            % and, therefore, the status matrix is built with a column
            % missing! This problem is solved by forcing the same
            % benchmark on the forward if a reverse is added. This problem 
            % does not occur when a line has only one direction (reverse or
            % forward) or if the benchmark is added to the forward first.
            
            % sort the benchmarks in the default direction
            % all other ordered arrays (deltas, residuals, status, etc)
            % will have cols that correspond to the order in the benchmarks
            % array
            [~,i] = ismember(tnames,{self.benchmarks.name}');
            self.benchmarks = self.benchmarks(i);

            % save the benchmark names
            self.start_benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks,tnames{1});
            self.end_benchmark   = CssBenchmark.ReturnBenchmark(self.benchmarks,tnames{end});
            
            % rename the line accordingly
            self.line_name = [self.start_benchmark.name '-' self.end_benchmark.name];
            
            % create/update the status matrix. It will contain a column for 
            % each instrument and a row for each benchmark. If a benchmark
            % observation for an instrument is bad, it gets flagged in this
            % matrix.
            
            self.status = true([length(self.instruments) length(self.benchmarks)]);
            
            % create/update the design matrix.
            self = self.UpdateDesign();
            self = self.GetDeltas();
        end
    end
    
    methods (Access = private)
        
        function self = load_g2grav_line(self, folder, scity, ecity)
            % private method to load a g2Gravity structure based on a
            % folder and a start and end city
            
            % sort the names of the cities so that we pickup both
            % directions during the scanning of the files
            r_city = sortrows([{scity};{ecity}]);
            
            % determine the possible directions (fwd, rev)
            cities = [{scity},{ecity}];
            
            % limit the search to files with this start and end city
            files1 = dir([folder '/*_' scity '_' ecity '_*.mat']);
            files2 = dir([folder '/*_' ecity '_' scity '_*.mat']);
            
            files = [files1; files2];
            % remove the directories
            files([files.isdir] == 1) = [];
            
            for i = 1:length(files)
                % split the fields using the delimiter
                fields = strsplit(files(i).name,'_');

                if length(fields) < 3
                    % to avoid fu$%@ DS_STORE files in the fu$%@ Macs
                    continue
                end
                
                point_name = self.rename_benchmark(lower(fields{1}));
                % get the sorted end points
                end_points = sortrows([fields(2);fields(3)]);
                % the line name should be the end points
                
                % verify that this file belongs to the start and end city
                % that has been requested
                
                if and(strcmp(r_city{1}, end_points{1}), strcmp(r_city{2}, end_points{2}))
                
                    load(fullfile(files(i).folder, files(i).name));
                    
                    % check the consistency of the monument name against
                    % the name of the matlab structure
                    if ~strcmp(point_name, self.rename_benchmark(lower(G.monument)))
                        disp(['WARNING! Inconsistent monument name found (' self.rename_benchmark(lower(G.monument)) ' - ' point_name ')'])
                    end
                    
                    if ~CssBenchmark.exists(self.benchmarks, point_name)
                        % if benchmark doesn't exist already, add it
                        self.benchmarks = [self.benchmarks; CssBenchmark(point_name, G.lat, G.lon, G.ht, G.voffset)];
                    end
                    
                    % obtain the handle to the current benchmark
                    current_benchmark = CssBenchmark.ReturnBenchmark(self.benchmarks, point_name);
                    
                    % determine the direction
                    % these directions are conventional w.r.t. the
                    % filenames in the mat files
                    if and(strcmp(fields{2}, cities{1}), strcmp(fields{3}, cities{2}))
                        direction = CssDirections.forward;
                    else
                        direction = CssDirections.reverse;
                    end
                    
                    % insert the observation
                    if ~isempty(current_benchmark)
                        self.observations = [self.observations; CssObservation(current_benchmark, direction, G.year, G.month, G.day, G.hour, G.minute, G.instrument, G.raw_data, G.reading, G.reduced_gravity)];
                    else
                        error(['Could not locate benchmark ' point_name ' in the benchmarks array.'])
                    end
                    
                    clear line_struct
                end
            end
            
            self = self.sort_benchmarks();
            
            self.comments = ['Line imported on ' datestr(datetime) ' from g2Gravity structures corresponding to cities ' scity ' - ' ecity];
            
        end
    end
    
    methods(Static)
        
        function [deltas, crossout] = check_deactivated_obs(deltas, status)
            
            [~,do] = find(status == false);
            crossout = [];
            for j = 1:length(do)
                % duplicate the previous delta to extend it over
                % the missing observations
                switch do(j)
                    case 1
                        deltas = [deltas(1); deltas];
                        coy = deltas(1);
                    case length(status)
                        deltas = [deltas; deltas(end)];
                        coy = deltas(end);
                    otherwise
                        if length(deltas) < do(j)-1
                            % all observations > do(j) until end are
                            % deactivated!
                            deltas = [deltas(1:end); deltas(end)];
                            coy = deltas(end);
                        else
                            deltas = [deltas(1:do(j)-1); deltas(do(j)-1); deltas(do(j):end)];
                            coy = deltas(do(j)-1);
                        end
                end
                % mark do to cross it out
                crossout = [crossout; do(j) coy];
            end
        end
        
        function out = rename_back(benchmark_name)
            translation = {'p492','p49b';
                            'gi06','p09c';
                            'k034','p16c';
                            'li00','p50c';
                            'p104','p10d';
                            'p105','p10e';
                            'pc07','p10g';
                            'pc08','p10h';
                            'pc21','p22a';
                            'pc32','p13b';
                            'pc33','p13c';
                            'pc34','p13d';
                            'pc51','p53a';
                            'pc52','p53b';
                            'pc54','p50d';
                            'pc55','p50e';
                            'pc74','p07d';
                            'pc82','p08b';
                            'pc83','p08c';
                            'pc84','p08d';
                            'pc85','p08e';
                            'pc86','p08f';
                            'pc88','p08h';
                            'pc93','p49c';
                            'pf21','p12a';
                            'pf91','p19a';
                            'pf92','p19b'};
            
            i = find(ismember(translation(:,2),benchmark_name));
            
            if i ~= 0 
                out = translation{i,1};
            else
                out = benchmark_name;
            end
        end
        
        function out = rename_benchmark(benchmark_name)
            % EXTRACTED FROM rename_existing_PC_points()
            translation = {'p492','p49b';
                            'gi06','p09c';
                            'k034','p16c';
                            'li00','p50c';
                            'p104','p10d';
                            'p105','p10e';
                            'pc07','p10g';
                            'pc08','p10h';
                            'pc21','p22a';
                            'pc32','p13b';
                            'pc33','p13c';
                            'pc34','p13d';
                            'pc51','p53a';
                            'pc52','p53b';
                            'pc54','p50d';
                            'pc55','p50e';
                            'pc74','p07d';
                            'pc82','p08b';
                            'pc83','p08c';
                            'pc84','p08d';
                            'pc85','p08e';
                            'pc86','p08f';
                            'pc88','p08h';
                            'pc93','p49c';
                            'pf21','p12a';
                            'pf91','p19a';
                            'pf92','p19b'};

            i = find(ismember(translation(:,1),benchmark_name));
            
            if i ~= 0 
                out = translation{i,2};
            else
                out = benchmark_name;
            end

        end
    end
end

