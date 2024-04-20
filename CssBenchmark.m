classdef CssBenchmark < matlab.mixin.Copyable
    % CSSBENCHMARKS 
    % Class to save the information regarding the locations of gravity
    % observations. This class just contains metadata of the benchmark
    % Can be initialized with lat lon h or x y z. The constructor detects
    % which one is being used and converts.
    % X Y Z should be given in meters
    % lat lon h should be given in deg deg meters.
    
    properties
        name
        lat
        lon
        height
        offset
        x
        y
        z
        absolute_g
        uncertainty
    end
    
    methods
        function self = CssBenchmark(name, varargin)
            self.name   = name;
            
            % save args to a temp var
            arg1 = varargin{1};
            arg2 = varargin{2};
            arg3 = varargin{3};
            
            if sqrt(arg1.^2 + arg2.^ 2 + arg3.^2) > 6300e3
                % XYZ: convert to lat lon
                lla = self.x2e([arg1 arg2 arg3]);
                
                self.lat    = lla(1)*180/pi;
                self.lon    = lla(2)*180/pi;
                self.height = lla(3);
                self.x      = arg1;
                self.y      = arg2;
                self.z      = arg3;
            else
                % lat lon h: convert to XYZ
                xyz = self.e2x([arg1*pi/180 arg2*pi/180 arg3]);
                
                self.x      = xyz(1);
                self.y      = xyz(2);
                self.z      = xyz(3);
                self.lat    = arg1;
                self.lon    = arg2;
                self.height = arg3;
            end
            
            self.offset = varargin{4}; % offset from point to gravity obs
        end
        
        function self = UpdateCoords(self, benchmark)
            % but leave offset untouched
            self.x = benchmark.x;
            self.y = benchmark.y;
            self.z = benchmark.z;
            self.lat = benchmark.lat;
            self.lon = benchmark.lon;
            self.height = benchmark.height;
        end
        
        function self = AssignAbsGravity(g, ug)
            self.absolute_g  = g;
            self.uncertainty = ug;
        end
        
        function self = RunPPP(self, path2rinex, AntH)
            % this function executes the PPP stuff (which has to be linked
            % in the ppp folder)
            % other dependencies of this function: crz2rnx (has to be in
            % the path directory)
            % links necessary in this folder:
            
            [local, rinex] = self.CopyCrinexRinex(self.name, path2rinex);

            [year,month,day,h] = self.ReadRinex(rinex, AntH);
            
            if ~isempty(year)
                [week, wkday] = gpsweek(year,month,day,0);
                
                week  = pad(num2str(week),4,'left','0');
                wkday = num2str(wkday);
                
                % get the orbits!
                path2orbits = ['ppp/orbits/' week];
                
                % determine the type of orbit
                orbits = self.find_orbit(path2orbits, week, wkday, local);
                
                if ~isempty(orbits)
                    % ready to run
                    if isnan(AntH)
                        AntH = h;
                    end
                    
                    self.make_files(local, rinex, self.name, orbits{1}, orbits{2}, orbits{3}, AntH);
                    % run the damm thing
                    system(['cd ' local '; ulimit -t 60; ../../ppp < input.inp' ]);
                    
                    % parse the output
                    [~,f,~] = fileparts(rinex);
                    self = self.parse_summary(local, f);
                else
                    error(['Could not find any orbit file for GPS week ' week ' day ' wkday '!'])
                end
            else
                error('Could not determine GPS date from RINEX file!')
            end
        end
        
        function self = parse_summary(self, local, filename)
            % parse the summary file
            summary = fileread(fullfile(local, [filename '.sum']));
            
            expression = '.*([XYZ]).\(m\)\s{1,14}[\s-]\w+.\w+\s{1,14}([\s-]\w+.\w+)\s{1,14}([\s-]\w+.\w+)';
            tokens = regexp(splitlines(summary), expression, 'tokens');
            tokens(cellfun(@isempty,tokens)) = [];
            tokens = [tokens{:}]';
            
            % get the x y z coordinates
            self.x = str2double(tokens{1}(2));
            self.y = str2double(tokens{2}(2));
            self.z = str2double(tokens{3}(2));
            
            lla = self.x2e([self.x self.y self.z]);
                
            self.lat    = lla(1)*180/pi;
            self.lon    = lla(2)*180/pi;
            self.height = lla(3);
        end
    end    
        
    methods (Static)
        
        function [folder, rinex] = CopyCrinexRinex(station, path2rinex)
            if ~exist('ppp/production', 'dir')
                mkdir('ppp/production')
            end
            
            if ~exist(['ppp/production/' station], 'dir')
                mkdir(['ppp/production/' station])
            else
                rmdir(['ppp/production/' station],'s');
                mkdir(['ppp/production/' station])
            end
            
            folder = ['ppp/production/' station];
            rinex  = [];
            
            % copy the rinex file to the folder
            copyfile(path2rinex,folder)
            
            [~, filename, extension] = fileparts(path2rinex);
            
            switch true
                case strcmp(extension, '.Z')
                    crinex = fullfile(folder,[filename extension]);
                    rinex  = fullfile(folder,[filename(1:end-1) 'o']);
                    
                    % uncompress and unhatanaka
                    [~,~] = system(['export PATH=$PATH:ppp/;ppp/crz2rnx -f -d ' crinex]);
                case strcmpi(extension(end), 'o')
                    rinex  = fullfile(folder,[filename extension]);
                    % already a rinex, nothing to be done
            end
        end
        
        function [year, month, day, height] = ReadRinex(rinex, varargin)
            year = [];
            month = [];
            day = [];
            height = [];
            % open rinex and determine the date
            fid   = fopen(rinex,'r');
            found = false(2,1);
            while and(~feof(fid), ~all(found))
                line = fgets(fid);
                first_obs  = sscanf(line,'%f %f %f %f %f %f %*3s %[TIME OF FIRST OBS]');
                ant_height = sscanf(line,'%f %f %f %[ANTENNA: DELTA H/E/N]');
                
                if length(first_obs) >= 23
                    % this might be it, make sure
                    key = char(first_obs(end-16:end))';
                    if strcmp(strip(key), 'TIME OF FIRST OBS')
                        year  = first_obs(1);
                        month = first_obs(2);
                        day   = first_obs(3);
                        
                        found(1) = true;
                    end
                end
                
                if length(ant_height) >= 23
                    key = char(ant_height(end-19:end))';
                    
                    if strcmp(strip(key), 'ANTENNA: DELTA H/E/N')
                        height  = ant_height(1);
                        
                        found(2) = true;
                    end
                end
                
                if strcmp(strip(sscanf(line,'%*s %[END OF HEADER]')), 'END OF HEADER')
                    break
                end
            end
            % check to see if requested antenna height is ~= than the
            % rinex header height. Change if different.
            if ~isempty(varargin)
                if str2double(height) ~= varargin{1}
                    new_antenna_h = sprintf('%14.4f%14.4f%14.4f%18sANTENNA: DELTA H\\/E\\/N',varargin{1},0,0,'');

                    [~, out] = system(['grep -n -E ''([- ]*\d{1,7})\.\d{4}\s{18}ANTENNA'' ' rinex ' | cut -d : -f 1']);
                    system(['sed -i '''' ''' strip(out) 's/.*/' new_antenna_h '/'' ' rinex]);
                end
            end
        end
        
        function orbits = find_orbit(path2orbits, week, wkday, local)
            % find orbits in the local orbits folder
            types = {'igs','ig2','igr','jpl','jp2','jpr'};
            
            orbits = {};
            
            for prod = {'sp3', 'clk', 'erp'}
                if strcmp(prod, 'erp')
                    wkday = '7';
                end
                
                for i = 1:length(types)
                    if exist(fullfile(path2orbits, [types{i} week wkday '.' prod{1} '.Z']),'file')
                        copyfile(fullfile(path2orbits, [types{i} week wkday '.' prod{1} '.Z']), local);
                        
                        system(['uncompress ' fullfile(local, [types{i} week wkday '.' prod{1} '.Z'])]);
                    else
                        if exist(fullfile(path2orbits, [types{i} week wkday '.' prod{1}]),'file')
                            copyfile(fullfile(path2orbits, [types{i} week wkday '.' prod{1}]), local);
                        end
                    end
                    
                    if exist(fullfile(local, [types{i} week wkday '.' prod{1}]),'file')
                        orbits = [orbits; fullfile(local, [types{i} week wkday '.' prod{1}])];
                        break
                    end
                end
            end
        end
        
        function make_files(path, rinex, name, sp3, clk, eop, ant_h)
            % only keep the file name
            [~,f,e] = fileparts(rinex);
            rinex = [f e];
            [~,f,e] = fileparts(sp3);
            sp3 = [f e];
            [~,f,e] = fileparts(clk);
            clk = [f e];
            [~,f,e] = fileparts(eop);
            eop = [f e];
            
            % make all the necessary files to run ppp
            fid = fopen(fullfile(path, 'gpsppp.def'),'w');
            
            fprintf(fid,'''LNG'' ''ENGLISH''\n');
            fprintf(fid,'''TRF'' ''../../gpsppp.trf''\n');
            fprintf(fid,'''SVB'' ''../../gpsppp.svb_gps_yrly''\n');
            fprintf(fid,'''PCV'' ''../../antennas.atx''\n');
            fprintf(fid,'''FLT'' ''../../gpsppp.flt''\n');
            fprintf(fid,'''OLC'' ''%s.olc''\n', name);
            fprintf(fid,'''MET'' ''../../gpsppp.met''\n');
            fprintf(fid,'''ERP'' ''%s''\n', eop);
            fprintf(fid,'''GSD'' ''%s''\n', 'OSU');
            fprintf(fid,'''GSD'' ''%s''\n', 'Gravity adjustment program');
            fclose(fid);
            
            fid = fopen(fullfile(path, 'commands.cmd'),'w');
            
            fprintf(fid,[''' UT DAYS OBSERVED                      (1-45)''               1\n', ...
                        ''' USER DYNAMICS         (1=STATIC,2=KINEMATIC)''               1\n', ...
                        ''' OBSERVATION TO PROCESS         (1=COD,2=C&P)''               2\n', ...
                        ''' FREQUENCY TO PROCESS        (1=L1,2=L2,3=L3)''               3\n', ...
                        ''' SATELLITE EPHEMERIS INPUT     (1=BRD ,2=SP3)''               2\n', ...
                        ''' SATELLITE PRODUCT (1=NO,2=Prc,3=RTCA,4=RTCM)''               2\n', ...
                        ''' SATELLITE CLOCK INTERPOLATION   (1=NO,2=YES)''               1\n', ...
                        ''' IONOSPHERIC GRID INPUT          (1=NO,2=YES)''               1\n', ...
                        ''' SOLVE STATION COORDINATES       (1=NO,2=YES)''               2\n', ...
                        ''' SOLVE TROP. (1=NO,2-5=RW MM/HR) (+100=grad) ''             105\n', ...
                        ''' BACKWARD SUBSTITUTION           (1=NO,2=YES)''               1\n', ...
                        ''' REFERENCE SYSTEM            (1=NAD83,2=ITRF)''               2\n', ...
                        ''' COORDINATE SYSTEM(1=ELLIPSOIDAL,2=CARTESIAN)''               2\n', ...
                        ''' A-PRIORI PSEUDORANGE SIGMA               (m)''           5.000\n', ...
                        ''' A-PRIORI CARRIER PHASE SIGMA             (m)''            .010\n', ...
                        ''' LATITUDE  (ddmmss.sss,+N) or ECEF X      (m)''          0.0000\n', ...
                        ''' LONGITUDE (ddmmss.sss,+E) or ECEF Y      (m)''          0.0000\n', ...
                        ''' HEIGHT (m)                or ECEF Z      (m)''          0.0000\n', ...
                        ''' ANTENNA HEIGHT                           (m)''          %6.4f\n', ...
                        ''' CUTOFF ELEVATION                       (deg)''           5.000\n', ...
                        ''' GDOP CUTOFF                                 ''          20.000\n'], ant_h);
             fclose(fid);
             
             fid = fopen(fullfile(path, 'input.inp'),'w');
             
             fprintf(fid,['%s\n', ...
                            'commands.cmd\n', ...
                            '0 0\n', ...
                            '0 0\n', ...
                            '%s\n', ...
                            '%s\n'], rinex, sp3, clk);
             fclose(fid);
        end
        
        function out = exists(benchmarks, name)
            for i = 1:length(benchmarks)
                if strcmp(benchmarks(i).name, name)
                    out = true;
                    return
                end
            end
            out = false;
        end
        
        function out = ReturnBenchmark(benchmarks, names)
%             for i = 1:length(benchmarks)
%                 if strcmp(benchmarks(i).name, name)
%                     out = benchmarks(i);
%                     return
%                 end
%             end
            if or(isempty(benchmarks),isempty(names))
                out = [];
            else
                % benchmarks should be returned in the same order used for the request
                % sort them in the same order
                [~, index] = ismember(names, {benchmarks.name});
                if index ~= 0
                    out = benchmarks(index);
                else
                    out = [];
                end
            end
        end
        
        function e = x2e(x,f,a)
            % x2e  Transform geocentric cartesian coords to ellipsoidal coords
            % Converts global cartesian coordinates into  ellipsoidal (geodetic)
            % coordinates for a single point. 
            %   Function e2x can be referenced in two ways
            %        e=x2e(x,f,a)  
            %   or   e=x2e(x)     in which case f,a take defaults values
            %
            %   INPUT ARGUMENTS
            %   x              global cartesian coordinate vector {x,y,z}
            %                  in same units as semi-major axis length (a).
            %   f  (optional)  flattening of ellipsoid 
            %                  default = 1/298.257223563  (WGS84)
            %   a  (optional)  semi-major axis length for ellipsoid, a.k.a.
            %                  equatorial radius - default is 6378137 meters (WGS84)
            %   OUTPUT ARGUMENT
            %   e              e(1)= ellipsoidal latitude in radians
            %                  e(2)= ellipsoidal longitude in radians
            %                  e(3)= ellipsoidal height in same units as "a"

            %  Note: Transformation from x to e is achieved iteratively
            %  Reference: Hoffman-Wellenhof et al (1994)  pages 255-258
            %  by Mike Bevis 1995
            tol=2*eps;  % convergence parameter
            itmax=10;    % max number of iterations
            if nargin~=1 & nargin~=3
              error('e2x takes one or three input arguments')
            end
            if nargin ==1
              % Set flattening (f) and semimajor axis length (a) to defaults
              f=1/298.257223563;  % Flattening in WGS-84 
              a=6378137.0;   % in meters. WGS-84 equatorial radius
            end
            esq=2*f-f^2;
            p=sqrt(x(1)^2+x(2)^2); 
            % compute first approximation for lat
            lat=atan2( x(3)/(1-esq),p );
            olat=lat;
            it=0; diff=2*tol;
            while diff > tol  & it < itmax     % iterate!
              it=it+1;
              N=a/sqrt(1 - esq*sin(lat)^2);
              h=p/cos(lat) - N;                      
              lat=atan2(x(3),p*(1 - esq*N/(N+h)) );  
              diff=abs(lat-olat);  
              olat=lat;
            end
            if it >= itmax
              error('x2e did not converge')
            end
            lon=atan2(x(2),x(1));
            e=x; % copy shape
            e(1)=lat; e(2)=lon; e(3)=h;
        end
        
        function x = e2x(e,f,a)
            % E2X converts ellipsoidal coordinates into global cartesian coordinates
            % Function e2x can be referenced in two ways
            %                x=e2x(e)                 uses a default ellipsoid
            %   or         x=e2x(e,f,a)            uses an ellipsoid specified by given (f,a)
            %
            %   INPUT ARGUMENTS
            %         e(1) = ellipsoidal latitude in radians
            %         e(2) = ellipsoidal longitude in radians
            %         e(3) = ellipsoidal height in same units as "a"
            %              f   = (optional)  flattening of ellipsoid
            %                      default value is f=1/298.257223563    (WGS-84 value)              
            %             a  = (optional)  semi-major axis length for ellipsoid, 
            %                     a.k.a. equatorial radius - default is 6378137 meters (WGS-84)
            %   OUTPUT ARGUMENT
            %              x  =   global cartesian coordinates {x,y,z}
            %                        in same units as "h" and "a".

            %  Reference: Hoffman-Wellenhof et al (1994)  pages 255-258
            %  by Mike Bevis 1995

            if nargin~=1 & nargin~=3
              error('e2x takes one or three input arguments')
            end
            if nargin ==1
              % Set flattening (f) and semimajor axis length (a) to defaults
              f=1/298.257223563;  % Flattening according to WGS-84    
              a=6378137.0;        % Equatorial radius in meters. WGS-84 value (same as GRS 1980)
            end
            esq=2*f-f^2;  
            slat=sin(e(1)); clat=cos(e(1)); slon=sin(e(2)); clon=cos(e(2));
            h=e(3);
            N=a/sqrt(1 - esq*slat^2);  % radius of curvature in prime vertical
            x=e; % copy shape of vector e
            x(1)=(N+h)*clat*clon;
            x(2)=(N+h)*clat*slon;
            x(3)=( (1-esq)*N + h)*slat;
        end
    end
end
