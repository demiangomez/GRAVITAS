classdef CssObservation < matlab.mixin.Copyable
    %CSSOBSERVATIONS class to hold one single observation. 
    % It references the benchmark class to know where the measurement 
    % was done. It also holds information like instrument, direction, raw_data, etc.
    
    properties (SetAccess = private)
        benchmark
        instrument
        direction
    end
    
    properties
        epoch
        timestamp
        year
        month
        day
        hour
        minute
        raw_data
        reading
        reduced_g
    end
    
    methods
        function self = CssObservation(benchmark, direction, year, month, day, hour, minute, instrument, raw_data, varargin)
            
            % compute the calculated fields
            self.timestamp = datetime(year,month,day,hour,minute,0);
            % cal2jd returns the same values as Kevin's AJD.
%           self.epoch     = cal2jd(year,month,day + hour/24 + minute/1440);    % ORIGINAL BEFORE FRANCO EDITED
            self.epoch     = self.JD1900(year,month,day,hour,minute);    % FRANCO EDITED 02/21/2024
            
            if length(varargin) == 2
                % all information provided
                self.benchmark  = benchmark;
                self.direction  = direction;
                self.year       = year;
                self.month      = month;
                self.day        = day;
                self.hour       = hour;
                self.minute     = minute;
                self.instrument = instrument;
                self.raw_data   = raw_data;
                self.reading    = varargin{1};
                self.reduced_g  = varargin{2} + benchmark.offset.*0.3086;
            elseif isempty(varargin)
                % no reading and reduced gravity, compute
                % in this case, instrument should be a cell array with {1}
                % being the name of the instrument and {2} the calibration
                % data
                if ~isstruct(instrument)
                    error('The instrument structure, which includes the calibration data, must be given for the instrument when invoking CssObservation without reduced_g')
                end
                
                self.benchmark  = benchmark;
                self.direction  = direction;
                self.year       = year;
                self.month      = month;
                self.day        = day;
                self.hour       = hour;
                self.minute     = minute;
                self.instrument = instrument.name;
                self.raw_data   = raw_data;
                
                self.compute_reduced_g(instrument.calibration);
                % apply the offset correction
                self.reduced_g = self.reduced_g + benchmark.offset.*0.3086;
            else
                error('Invalid number of arguments to initialize CssObservation');
            end
        end
        
        function self = UpdateTimeStamp(self, calibration, varargin)
            switch length(varargin)
                case 1
                    self.timestamp = varargin{1};
%                   self.epoch     = cal2jd(year(varargin{1}),month(varargin{1}),day(varargin{1}) + hour(varargin{1})/24 + minute(varargin{1})/1440);  % ORIGINAL BEFORE FRANCO EDITED
                    self.epoch     = self.JD1900(year(varargin{1}),month(varargin{1}),day(varargin{1}),hour(varargin{1}),minute(varargin{1}));   % FRANCO EDITED 02/21/2024
                case 5
                    self.timestamp = datetime(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},0);
%                   self.epoch     = cal2jd(varargin{1},varargin{2},varargin{3} + varargin{4}/24 + varargin{5}/1440);   % ORIGINAL BEFORE FRANCO EDITED
                    self.epoch     = self.JD1900(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});   % FRANCO EDITED 02/21/2024
                otherwise
                    error('invalid number of arguments!')
            end
            
            self.compute_reduced_g(calibration);
            self.reduced_g = self.reduced_g + self.benchmark.offset.*0.3086;
        end

%%%%% NEW FUNCTION TO RECTIFY REDUCED_G WITH CORRECT TIDE MODEL (DEMIAN, FEB 2024) %%%%%
        function self = ForceUpdate(self)
            tide_wrong  = self.gtide(self);
            self.epoch  = self.JD1900(self.year, self.month, self.day, self.hour, self.minute);
            tide        = self.gtide(self);
            self.reduced_g = self.reduced_g - tide_wrong + tide;
        end
%%%%%%%%%%%%%%%%%

        function self = UpdateReading(self, calibration, reading_index, value)
            if isnumeric(value)
                self.raw_data(reading_index) = value;
                self.compute_reduced_g(calibration);
                self.reduced_g = self.reduced_g + self.benchmark.offset.*0.3086;
            else
                error('value has to be a number')
            end
        end
        
        function self = UpdateOffset(self, calibration, value)
            if isnumeric(value)
                self.benchmark.offset = value;
                self.compute_reduced_g(calibration);
                self.reduced_g = self.reduced_g + self.benchmark.offset.*0.3086;
            else
                error('value has to be a number')
            end
        end
        
        function self = ChangeBenchmark(self, newBenchmark, calibration)
            % change the read only property benchmark
            if self.benchmark.offset ~= newBenchmark.offset
                self.benchmark = newBenchmark;
                self = self.UpdateOffset(calibration, newBenchmark.offset);
            else
                self.benchmark = newBenchmark;
            end 
        end
        
        function self = SwitchDirection(self)
            if self.direction == CssDirections.forward
                self.direction = CssDirections.reverse;
            else
                self.direction = CssDirections.forward;
            end
        end
        
        function self = ChangeInstrument(self, newInstrument, calibration)
            self.instrument = newInstrument;
            self.compute_reduced_g(calibration);
            self.reduced_g = self.reduced_g + self.benchmark.offset.*0.3086;
        end
        
        function self = compute_reduced_g(self, calibration)
            % from Kevin's function reduceLCR.m
            %LCR   Reduction of Lacoste Romberg gravimeter readings
            % This function applies the calibration curve established by the manufacturer 
            % for each gravimeter, and also applies a gravity tide correction

            g_reading = mean(self.raw_data);
            self.reading = g_reading;
            
            np=length(calibration);
            gx=100*(0:1:np-1);       % the grvmtr readings corresponding to the values in the calib curve
            if or(g_reading < gx(1), g_reading > gx(end))
                errStruct.identifier = 'CssObservation:compute_reduced_g:calibration_error';
                errStruct.message = ['The g_reading falls outside the calibrated range of this instrument! (' num2str(gx(1)) ' to ' num2str(gx(end)) '; entered value: ' num2str(g_reading) ')'];
                error(errStruct)
            end

            % apply the calibration curve
            gg=spline(gx, calibration, g_reading);
            
            % compute gravity tide
            tide = self.gtide(self);

            % apply tidal correction
            self.reduced_g = gg + tide; % note that the sign convention for the gravity tide
                                        % means it must be added rather than subtracted
        end
    end
    
    methods(Static)
%%%%%% START OF WHAT FRANCO ADDED ON FEB 21, 2024
        function jd = AJD(y,m,d,ut)
            %AJD  Julian Day (standard astronomical definition)
            % using a fast algorithm after Montenbruck (1984) that is 
            % valid for dates between March 1900 and Feb 2100.
            %
            % USAGE
            %
            %    jd = AJD(y,m,d)
            %
            %  gives the JD for calender year y, month m, and integer or real
            %  valued day d  (in UT), or
            %
            %   jd = AJD(y,m,d,t)   
            %  
            %  as before but d is integer day, and t is time in UT hours
            %
            %  EXAMPLE   jd = AJD(1988,6,19.5)  results in jd = 2447332
            
            % Version 1.0        29 Nov 98      Mike Bevis
            
            if nargin==3
               ut=rem(d,1)*24;
               d=floor(d);
            end
            if nargin==4
               if rem(d,1)~=0
                  error('give four input arguments, d must be integer')
               end
            end
            if m<=2
               m=m+12;
               y=y-1;
            end
            jd=floor(365.25*y) + floor(30.6001*(m+1)) + d + ut/24 + 1720981.5;
        end

        function epoch = JD1900(y,m,d,hr,mn,sc)
            %JD1900  Time in days since Greenwich mean noon on 31 Dec 1899 (Longman,1959)
            % corresponding to the conventional astronomical Julian Day - 2415020
            % 
            % USAGE                                 
            %   epoch = JD1900(y,m,d,hr,mn,sc)       
            %   epoch = JD1900(y,m,d,hr,mn)
            %   epoch = JD1900(y,m,d,hr)
            %   epoch = JD1900(y,m,d)
            %  
            % where y is year, m is month, d is day, and the UTC time of day is given 
            % in some combination of hours (hr), minutes (mn) and seconds (sc). In each 
            % of the four patterns of usage, listed above, all the input arguments 
            % must be integers except the last argument which may be (and often 
            % will be) real. To compute the epoch corresponding to 19:30 UTC on 4 May
            % 1988, you could use
            %
            % epoch = JD1900(1988,6,4,19.5)        results in jd = 32297.3125
            %
            % This way of characterizing dates and times is used  with function gtide.m
            % which computes gravity tides.
            %
            % This routine uses function AJD.m to compute conventional Julian Days. 
            % AJD uses a fast algorithm, after Montenbruck (1984) that is valid
            % only for dates between March 1900 and Feb 2100. So this limits function
            % JD1900 to the same range of dates.
            
            % Version 1.2        17 April 2009     Mike Bevis
            
            if y<=1900 | y >=2100
                if y<1900 | (y==1900 & m<3)
                    error('Date must be after 1 March 1900')
                end
                if y>2100 | (y==2100 & m>1)
                    error('Date must be before 1 February 2100')
                end
            end
            if (m<0 | m>12) | (d<0 | d>31)
                error('Invalid value for month or day of month')
            end
            if nargin==6
                ut=hr + mn/60 + sc/3600;
            elseif nargin==5
                ut=hr + mn/60;
            elseif nargin==4
                ut=hr;
            elseif nargin==3
                ut=rem(d,1)*24;
                d=floor(d);
            else
                error('There must be at least 3 input arguments');
            end     
               
            epoch = CssObservation.AJD(y,m,d,ut) - 2415020;
        end

%%%%%% END OF WHAT FRANCO ADDED ON FEB 21, 2024

        function tide = gtide(observation)
            %GTIDE Computes the gravity tide with or w/o an earth tide correction
            % This matlab function follows the general structure of the
            % fortran subroutine TIDE.FOR written by R. Forsberg, which is based
            % on the equations of Longman (1959). The (optional) earth tide correction
            % is that used in the R. Forsberg's fortran code GRREDU.FOR
            %
            % USAGE:   tide = gtide(time,lat,lon,height,etidec)
            %
            % INPUT ARGUMENTS:
            %    time   date and time stated in days since UTC noon 31 Dec 1899 
            %             (see matlab function JD1900.m)
            %     lat   station latitude in degrees (N is +ve)
            %     lon   station longitude in degrees (E is +ve)
            %  height   station height in meters. Technically this should be an
            %             orthometric height (i.e. height above sea level) but
            %             if ellipsoidal height is used instead, the resulting
            %             error is negligible
            % etidec    determines if an earth tide correction (ETC) is applied
            %             if etidec=1     ETC is applied
            %             if etidec=0     ETC is not applied
            %
            % OUTPUT ARGUMENT:
            %    tide   gravity tide in milligals
            %
            % Input Argument Sizes
            %    The physical arguments time,lat,lon and height must follow these
            % rules. They can (1) all be scalar, or (2) one of the can be a vector
            % and the others scalar, or (3) all four arguments can be vectors, but 
            % they must have the same size and shape - in which can tide(i) will
            % be computed at [ lat(i) lon(i) ht(i)] at time(i).
            %
            % See also function JD1900.m

            % Version 1.0           Michael Bevis             17 April 2009

            time = observation.epoch;
            lat  = observation.benchmark.lat;
            lon  = observation.benchmark.lon;
            height = observation.benchmark.height;
            etidec = 1;
            
            if and(etidec ~=0, etidec~=1)    
                error('etidec must be 0 or 1')
            end

            dtr    = pi/180;
            e      = 0.054899720;
            c      = 3.84402e10;
            c1     = 1.495e13;
            aprim  = 1/(c*(1-e*e));
            i      = 0.08979719;
            omega  = 0.4093146162;
            ss     = 1.993e33;
            mm     = 7.3537e25;
            my     = 6.670e-8;
            m      = 0.074804;

            %  computation point
            coslambda  = cos(lat*dtr);
            sinlambda  = sin(lat*dtr);
            r      = 6.378270e8./sqrt(1+0.006738*sinlambda.^2)+ height*100;
            ll     = lon*dtr;

            % some fundamental time-predictable quantities
            tt  = time/36525;
            tt2 = tt.^2;
            tt3 = tt.^3;
            s   = 4.720023438 + 8399.7093*tt + 4.40695e-5*tt2 + 3.29e-8*tt3;
            p   = 5.835124721 + 71.018009*tt - 1.80546e-4*tt2 - 2.181e-7*tt3;
            h   = 4.881627934 + 628.33195*tt + 5.27960e-6*tt2;
            N   = 4.523588570 - 33.757153*tt + 3.67488e-5*tt2 + 3.870e-8*tt3;
            p1  = 4.908229467 + 3.0005264e-2*tt + 7.9024e-6*tt2 + 5.81e-8*tt3;
            e1  = 0.01675104  - 4.18e-5 *tt - 1.26e-7 * tt2;

            %  reciproc distances
            a1prim= 1./(c1*(1-e1.^2));
            resd  = 1/c + aprim*e*cos(s-p) + aprim*e*e*cos(2*(s-p)) ...
                    + 15e0/8*aprim*m*e*cos(s-2*h+p) + aprim*m*m*cos(2*(s-h));
            resdd  = 1/c1 + a1prim.*e1.*cos(h-p1);

            %  longitude of moons ascending node
            cosii = cos(omega)*cos(i)   - sin(omega)*sin(i)*cos(N);
            sinii = sqrt(1-cosii.^2);
            ii    = atan(sinii./cosii);
            ny    = asin(sin(i)*sin(N)./sinii);

            %  longitude and rigth ascension
            t     = 2 * pi * rem(time,1) + ll;
            ksi1  = t + h;
            ksi   = ksi1 - ny;
            l1    = h + 2*e1.*sin(h-p1);
            alfa  = 2 * atan( (sin(omega)*sin(N)./sinii) ./ (1 + cos(N).*cos(ny) ...
                     + sin(N).*sin(ny)*cos(omega)) );
            sigma = s - N + alfa;
            l     = sigma + 2*e*sin(s-p) + 5e0/4*e*e*sin(2*(s-p)) ...
                     + 15e0/4*m*e*sin(s - 2*h + p) + 11e0/8*m*m*sin(2*(s -h));

            %  zenith angles
            costheta = sinlambda.*sinii.*sin(l) + coslambda.*(cos(ii/2).^2.*cos(l-ksi) ...
                        +  sin(ii/2).^2.*cos(l+ksi));
            cosphi   = sinlambda*sin(omega).*sin(l1) ...
                       + coslambda.*(cos(omega/2)^2*cos(l1-ksi1) ...
                       + sin(omega/2)^2.*cos(l1+ksi1));

            %  gravities
            gs    = my.*ss.*r.*resdd.^3.*(3*cosphi.^2 - 1);
            gm    = my.*mm.*r.*resd.^3.*(3*costheta.^2 - 1) ...
                    + (3/2)*my.*mm.*r.^2.*resd.^4.*(5*costheta.^3 - 3*costheta);
            g0    = gm + gs;

            %  transformation from the cgs unit gal to mgal
            tide   = g0 * 1000;

            % perform earth tide correction if desired
            if etidec==1
                delta = 1.14;
                tide = tide*delta + 0.00483 - 0.01573*cos(lat.*dtr).^2;
            end

        end
    end
end

