classdef CssVisualization
    %CSSVISUALIZATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pAxes
        lonlim
        latlim
        maplon
        maplat
        maph
        benchmarks
        lines
        lines_status
        network
        plot_line_handles
        plot_benchmark_handles
        plot_node_handles
        selected_line
        selected_benchmark
        LineNameList
        stdmax
        for_color
        sel_color
        err_color
        nod_color
        residualAxes
    end
    
    methods
        function self = CssVisualization(lines, benchmarks, paxes, for_color, sel_color, nod_color, err_color)
            % visualization class to plot the relevant line and network
            % information
            
            self.pAxes = paxes;
            
            if ~isempty(lines)
                self.lines        = lines;
                self.lines_status = true(length(lines),1);
                self.network      = CssNetwork(lines, benchmarks);
            end
            
            self.selected_line        = [];
            self.selected_benchmark   = [];
            self.LineNameList         = {};
            self.stdmax               = Inf;
            self.benchmarks           = benchmarks;
            
            % selection colors
            self.for_color = for_color;
            self.sel_color = sel_color;
            self.err_color = err_color;
            self.nod_color = nod_color;
           
            self.latlim = [min([benchmarks.lat])-0.1 max([benchmarks.lat]+0.1)];
            self.lonlim = [min([benchmarks.lon])-0.1 max([benchmarks.lon]+0.1)];
            
            % get the topography information
            [self.maph,self.maplon,self.maplat] = m_etopo2([self.lonlim self.latlim]);
            
            self.PlotBaseMap()
        end
        
        function PlotBaseMap(self)
            axes(self.pAxes)
            
            % clear axes
            cla
            % invoke the marcator projection
            m_proj('mercator','lon',self.lonlim,'lat',self.latlim);
            
            % plot basemap
            m_pcolor(self.maplon,self.maplat,self.maph);
            shading flat
            shading interp
            colormap terrain
            % make the topography range fixed
            caxis([-8000 8000])
            % boundaries
            m_gshhs('lc1','color','w')
            m_gshhs('lb1','color','w')
            
            hold on
            
            axis tight
            axis equal
            %m_grid()
            colorbar();
        end
        
        function self = PlotBenchmarks(self)
            % activate the default axes before ploting
            axes(self.pAxes)
            
            reorder_items = [];
            
            for i = 1:length(self.benchmarks)
                lat = self.benchmarks(i).lat;
                lon = self.benchmarks(i).lon;
                
                [x,y] = m_ll2xy(lon,lat);
                if isempty(self.benchmarks(i).absolute_g)
                    self.plot_benchmark_handles{i} = plot(x, y, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', self.for_color);
                else
                    reorder_items = [reorder_items; i];
                    self.plot_benchmark_handles{i} = plot(x, y, '^', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'Color', self.for_color);
                end
            end
            uistack([self.plot_benchmark_handles{reorder_items}],'top')
        end
        
        function self = SelectBenchmark(self, Benchmarks)
            % find the index of benchmark
            BenchmarkIndex = find(ismember(self.benchmarks,Benchmarks));
            
            self.selected_benchmark = BenchmarkIndex;
            
            if ~isempty(self.plot_benchmark_handles)
                for i = 1:length(self.benchmarks)
                    if ~ismember(i, BenchmarkIndex)
                        self.plot_benchmark_handles{i}.Color = self.for_color;
                        self.plot_benchmark_handles{i}.LineWidth = 0.5;
                    else
                        self.plot_benchmark_handles{i}.Color = self.sel_color;
                        self.plot_benchmark_handles{i}.LineWidth = 2;
                    end
                end
            else
                error('No handle to benchmark plot. The benchmark were not drawn yet!')
            end
        end
        
        function self = PickBenchmark(self, sx, sy)
            % pick a benchmark from the plot
            axes(self.pAxes)
            
            x = xlim;
            y = ylim;
            
            % adjust a radius of selection based on the max distance span
            % the map gives
            kmdist = m_xydist(x,y);
            
            % first, try with kmdist/200
            
            [lon, lat] = m_xy2ll(sx,sy);
    
            % find a point that is within 10 km of the clicked location
            d = m_idist(lon,lat,[self.benchmarks.lon], [self.benchmarks.lat]);
            
            selection = self.benchmarks(d < (kmdist*1000)/200);
            
            self = self.SelectBenchmark(selection);
        end
        
        function self = PickLine(self, sx, sy)
            % loop through the lines to find the intersecting line with the
            % crosshairs at sx,sy
            axes(self.pAxes)
            
            xc = xlim;
            yc = ylim;
            
            % crosshairs as a function of the zoom scale
            kmdist = m_xydist(xc,yc);
            xh = (kmdist*1000)/200;
            
            % crosshair lat and lon
            [lon_xh, lat_xh] = m_xy2ll(sx,sy);
            
            % make a romboidal figure fro the crosshair
            [lon(1),lat(1)] = m_fdist(lon_xh,lat_xh,  0,xh);
            [lon(2),lat(2)] = m_fdist(lon_xh,lat_xh, 90,xh);
            [lon(3),lat(3)] = m_fdist(lon_xh,lat_xh,180,xh);
            [lon(4),lat(4)] = m_fdist(lon_xh,lat_xh,270,xh);
            [lon(5),lat(5)] = m_fdist(lon_xh,lat_xh,  0,xh);
            
            if any(lon > 180)
                lon(lon > 180) = lon - 360;
            end
            % convert crosshair to x,y
            [xh, yh] = m_ll2xy(lon',lat');
            
            for i = 1:length(self.lines)
                line_bench_lat = [self.lines(i).benchmarks.lat]';
                line_bench_lon = [self.lines(i).benchmarks.lon]';
                
                [x,y] = m_ll2xy(line_bench_lon,line_bench_lat);
                
                [x0, ~] = self.intersections(x, y, xh, yh);
                
                if ~isempty(x0)
                    % intersection found in this line
                    self = self.SelectLine(i);
                    return
                end
            end
        end
        
        function self = PlotLinesAdjustment(self)
            % this function plots the lines but ignores information related
            % to outliers in residuals. It also assumes that a line without
            % closure (all NaNs in residuals) is actually a deactivated
            % line, so it removes it and marks it with a special color
            
            % activate the default axes before ploting
            axes(self.pAxes)
            
            for i = 1:length(self.lines)
                lat = [self.lines(i).benchmarks.lat];
                lon = [self.lines(i).benchmarks.lon];
                
                if all(isnan(cell2mat(self.lines(i).residuals)))
                    % all gravimeters in the line have no closure!! (drift
                    % = NaN)
                    [x,y] = m_ll2xy(lon,lat);
                    self.plot_line_handles{i} = plot(x, y, '--o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', self.err_color);
                    self.LineNameList(i) = {sprintf('<HTML><BODY bgcolor="%s">%s', 'yellow', [self.lines(i).line_name ' (' strjoin(str(self.lines(i).directions)) ')'])};
                    
                elseif ~self.lines_status(i)
                    % show deactivated line
                    [x,y] = m_ll2xy(lon,lat);
                    self.plot_line_handles{i} = plot(x, y, '-o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', [0 0 1]);
                    self.LineNameList(i) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', [self.lines(i).line_name ' (' strjoin(str(self.lines(i).directions)) ')'])};
                    
                else
                    % normal, activated line
                    [x,y] = m_ll2xy(lon,lat);
                    self.plot_line_handles{i} = plot(x, y, '-o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', self.for_color);
                    self.LineNameList(i) = {[self.lines(i).line_name ' (' strjoin(str(self.lines(i).directions)) ')']};
                end
            end
        end
        
        function self = PlotLines(self, stdmax)
            % plot the lines and output a list of the line_names colored in
            % red for lines with residuals > 3*stdmax
            
            % activate the default axes before ploting
            axes(self.pAxes)
            
            % update the stdmax value
            if ~isempty(stdmax)
                self.stdmax = stdmax;
            end
            
            for i = 1:length(self.lines)
                lat = [self.lines(i).benchmarks.lat];
                lon = [self.lines(i).benchmarks.lon];
                
                if any(abs(cell2mat(self.lines(i).residuals)) > 3*self.stdmax)
                    % line with a residual > than 3*stdmax
                    [x,y] = m_ll2xy(lon,lat);
                    self.plot_line_handles{i} = plot(x, y, '-o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', self.err_color);
                    self.LineNameList(i) = {sprintf('<HTML><BODY bgcolor="%s">%s', 'red', self.lines(i).line_name)};
                elseif all(isnan(cell2mat(self.lines(i).residuals)))
                    % all gravimeters in the line have no closure!! (drift
                    % = NaN)
                    [x,y] = m_ll2xy(lon,lat);
                    self.plot_line_handles{i} = plot(x, y, '--o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', self.err_color);
                    self.LineNameList(i) = {sprintf('<HTML><BODY bgcolor="%s">%s', 'yellow', self.lines(i).line_name)};
                else
                    % line with OK residuals
                    [x,y] = m_ll2xy(lon,lat);
                    self.plot_line_handles{i} = plot(x, y, '-o', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'Color', self.for_color);
                    self.LineNameList(i) = {self.lines(i).line_name};
                end
            end
        end
        
        function self = UpdateStdMax(self, stdmax)
            % update the standard deviation and update the plot
            self.stdmax = stdmax;
            self.SelectLine(self.selected_line);
            
            % update the display list
            for i = 1:length(self.lines)
                if any(abs(cell2mat(self.lines(i).residuals)) > 3*self.stdmax)
                    % line with a residual > than 3*stdmax
                    self.LineNameList(i) = {sprintf('<HTML><BODY bgcolor="%s">%s', 'red', self.lines(i).line_name)};
                elseif all(isnan(cell2mat(self.lines(i).residuals)))
                    % all gravimeters in the line have no closure!! (drift
                    % = NaN)
                    self.LineNameList(i) = {sprintf('<HTML><BODY bgcolor="%s">%s', 'yellow', self.lines(i).line_name)};    
                else
                    % line with OK residuals
                    self.LineNameList(i) = {self.lines(i).line_name};
                end
            end
        end
        
        function self = PlotNodes(self)
            % plot the nodes and absolute benchmarks            
            for i = 1:length(self.benchmarks)
                lat = self.benchmarks(i).lat;
                lon = self.benchmarks(i).lon;
                
                if and(ismember(self.benchmarks(i).name, {self.network.nodes.name}), isempty(self.benchmarks(i).absolute_g))
                    % is a node, but not an ansolute benchmark
                	self.plot_node_handles{i} = m_plot(lon, lat, 's', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'Color', self.nod_color);
                    
                elseif ~isempty(self.benchmarks(i).absolute_g)
                    % is a node and/or an absolute benchmark
                    self.plot_node_handles{i} = m_plot(lon, lat, '^', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'Color', self.nod_color);
                end
            end
        end
        
        function self = PlotResiduals(self, residuals, outliers)
            axes(self.pAxes)
            
            c = ajet;
            m = size(c,1);
            cmin = nanmin(residuals.residual);  % Minimum color value
            cmax = nanmax(residuals.residual);  % Maximum color value

            stdev = 3*nanstd(residuals.residual(outliers));
            
            for i = 1:length(residuals.start_benchmark)
                if and(~isnan(residuals.residual(i)), abs(residuals.residual(i)) > stdev)
                    v = min(m,round((m-1)*(residuals.residual(i)-cmin)/(cmax-cmin))+1);
                    
                    [x,y] = m_ll2xy([residuals.start_benchmark(i).lon residuals.end_benchmark(i).lon],[residuals.start_benchmark(i).lat residuals.end_benchmark(i).lat]);
                    plot(x,y,'color',c(v,:), 'LineWidth',2, 'LineStyle', '--');
                    %self.plot_benchmark_handles{i}.Color = c(v,:);
                    %self.plot_benchmark_handles{i}.MarkerFaceColor = c(v,:);
                %else
                    %self.plot_benchmark_handles{i}.Color = self.for_color;
                    %self.plot_benchmark_handles{i}.MarkerFaceColor = [1 1 1];
                end
            end
            
            % create a colorbar for this data
            if isempty(self.residualAxes)
                ax = axes;
                self.residualAxes = ax;
                colormap(ax, 'ajet');
                axis(ax,'off')
                colorbar(ax,'Position',[.95 .11 .02 .815]);
                uistack(ax,'bottom')
            else
                ax = self.residualAxes;
            end
            caxis(ax,[cmin cmax]);
        end
        
        function self = ToggleLineActivation(self, LineIndex)
            % function to activate and deactivate a line
            
            if ~all(isnan(cell2mat(self.lines(LineIndex).residuals)))
                % do not activate of deactivate if the line is already
                % deactivated by having no deltas
                self.lines_status(LineIndex) = ~self.lines_status(LineIndex);

                % repaint lines and regenerate list for listbox
                if self.lines_status(LineIndex)
                    self.plot_line_handles{LineIndex}.Color = self.for_color;
                    self.plot_line_handles{LineIndex}.LineWidth = 0.5;
                    self.LineNameList(LineIndex) = {[self.lines(LineIndex).line_name ' (' strjoin(str(self.lines(LineIndex).directions)) ')']};
                else
                    self.plot_line_handles{LineIndex}.Color = [0 0 1];
                    self.plot_line_handles{LineIndex}.LineWidth = 2;
                    self.LineNameList(LineIndex) = {sprintf('<HTML><strong><FONT color="%s">%s', 'blue', [self.lines(LineIndex).line_name ' (' strjoin(str(self.lines(LineIndex).directions)) ')'])};
                end
            end
        end
        
        function self = SelectLineAdjustment(self, LineIndex)
            % paint the line corresponding to LineIndex in the select color
            if ~isempty(self.plot_line_handles)
                for i = 1:length(self.lines)
                    if i ~= LineIndex
                        if self.lines_status(i)
                            % if the line is activated, paint normal
                            % otherwise, paint in blue thicker
                            self.plot_line_handles{i}.Color = self.for_color;
                            self.plot_line_handles{i}.LineWidth = 0.5;
                        else
                            self.plot_line_handles{i}.Color = [0 0 1];
                            self.plot_line_handles{i}.LineWidth = 2;
                        end
                    else
                        self.plot_line_handles{i}.Color = self.sel_color;
                        self.plot_line_handles{i}.LineWidth = 2;
                        self.selected_line = LineIndex;
                    end
                end
                
                % check if the line has any nodes and/or absolute gravity
                % benchmarks
                for i = 1:length(self.plot_node_handles)
                    self.plot_node_handles{i}.Color = self.nod_color;
                    self.plot_node_handles{i}.LineWidth = 0.5;
                end
                
                % find the indecies of the line benchmarks in the
                % benchmarks array
                bi    = ismember({self.benchmarks.name}, {self.lines(LineIndex).benchmarks.name});
                
                % find which ones have a node_handle object
                bi   = find(and(~isempty(self.plot_node_handles), bi));
                
                for i = 1:length(bi)
                    self.plot_node_handles{bi(i)}.Color = self.sel_color;
                    self.plot_node_handles{bi(i)}.LineWidth = 2;
                end
            else
                error('No handle to lines plot. The lines were not drawn yet!')
            end
        end
        
        function self = SelectLine(self, LineIndex)
            % paint the line corresponding to LineIndex in the select color
            if ~isempty(self.plot_line_handles)
                for i = 1:length(self.lines)
                    if i ~= LineIndex
                        if any(abs(cell2mat(self.lines(i).residuals)) > 3*self.stdmax)
                            % the line has some outlier deltas, paint in
                            % err_color
                            self.plot_line_handles{i}.Color = self.err_color;
                            self.plot_line_handles{i}.LineWidth = 0.5;
                        else
                            self.plot_line_handles{i}.Color = self.for_color;
                            self.plot_line_handles{i}.LineWidth = 0.5;
                        end
                    else
                        self.plot_line_handles{i}.Color = self.sel_color;
                        self.plot_line_handles{i}.LineWidth = 2;
                        self.selected_line = LineIndex;
                    end
                end
                
                % check if the line has any nodes and/or absolute gravity
                % benchmarks
                for i = 1:length(self.plot_node_handles)
                    self.plot_node_handles{i}.Color = self.nod_color;
                    self.plot_node_handles{i}.LineWidth = 0.5;
                end
                
                % find the indecies of the line benchmarks in the
                % benchmarks array
                bi    = ismember({self.benchmarks.name}, {self.lines(LineIndex).benchmarks.name});
                
                % find which ones have a node_handle object
                bi   = find(and(~isempty(self.plot_node_handles), bi));
                
                for i = 1:length(bi)
                    self.plot_node_handles{bi(i)}.Color = self.sel_color;
                    self.plot_node_handles{bi(i)}.LineWidth = 2;
                end
            else
                error('No handle to lines plot. The lines were not drawn yet!')
            end
        end
        
        function PlotLineInfo(self, LineIndex, force_plot)
            % force_plot forces the plot even if figures are not open

            h1 = CssVisualization.find_figure_handles(self.lines(LineIndex).line_name);

            if isempty(h1)
                if ~force_plot
                    return
                end
                sc = get(groot,'ScreenSize');
                figure('Position',sc);
            else
                set(0, 'CurrentFigure', h1)
                clf
            end

            set(gcf, 'Name', ['Information for line ' self.lines(LineIndex).line_name]);

            for i = 1:length(self.lines(LineIndex).instruments)
                subplot(length(self.lines(LineIndex).instruments),2,2*i-1)
                % to maximize the space usage
                pos = get(gca, 'Position');
                pos(1) = 0.055;
                pos(3) = 0.42;
                set(gca, 'Position', pos)
                
                % plot deltas
                plot(self.lines(LineIndex), self.lines(LineIndex).instruments{i}, self.stdmax)
            end

            % make another figure for the raw measurements

            for i = 1:length(self.lines(LineIndex).instruments)
                subplot(length(self.lines(LineIndex).instruments),2,2*i)
                % to maximize the space usage
                pos = get(gca, 'Position');
                pos(1) = 0.555;
                pos(3) = 0.42;
                set(gca, 'Position', pos)
                
                % plot raw measurements
                plotRaw(self.lines(LineIndex), self.lines(LineIndex).instruments{i})
            end
        end
        
        function PlotLineTime(self, LineIndex, force_plot)
            h1 = CssVisualization.find_figure_handles(self.lines(LineIndex).line_name);

            if isempty(h1)
                if ~force_plot
                    return
                end
                sc = get(groot,'ScreenSize');
                figure('Position',sc);
            else
                set(0, 'CurrentFigure', h1)
                clf
            end

            set(gcf, 'Name', ['Information for line ' self.lines(LineIndex).line_name]);
            
            % plot the spatial information first
            for i = 1:length(self.lines(LineIndex).instruments)
                subplot(length(self.lines(LineIndex).instruments),1,i)
                % to maximize the space usage
                pos = get(gca, 'Position');
                pos(1) = 0.055;
                pos(3) = 0.90;
                set(gca, 'Position', pos)
                
                % make the time-line plot
                plotLineTime(self.lines(LineIndex), self.lines(LineIndex).instruments{i})
            end
        end
        
        function PlotComparison(self, LineIndex, force_plot)
            h1 = CssVisualization.find_figure_handles(self.lines(LineIndex).line_name);

            if isempty(h1)
                if ~force_plot
                    return
                end
                sc = get(groot,'ScreenSize');
                figure('Position',sc);
            else
                set(0, 'CurrentFigure', h1)
                clf
            end

            set(gcf, 'Name', ['Information for line ' self.lines(LineIndex).line_name]);
            
            % make the time-line plot
            plotComparison(self.lines(LineIndex))
        end
        
        function ZoomInLine(self)
            axes(self.pAxes)
            
            % find selected line's handle
            if ~isempty(self.selected_line)
                x = [min(self.plot_line_handles{self.selected_line}.XData) max(self.plot_line_handles{self.selected_line}.XData)];
                y = [min(self.plot_line_handles{self.selected_line}.YData) max(self.plot_line_handles{self.selected_line}.YData)];
                zoom on;
                
                % set the correct reset point
                set(gca,'ylimMode', 'auto')
                set(gca,'xlimMode', 'auto')
                zoom(gcf,'reset')
                
                % get the aspect ratio of the figure
                a = diff(xlim)/diff(ylim);
                % get the aspect ratio of the line

                
                if diff(x) < diff(y)
                    d = a*diff(y) - diff(x);
                    
                    axis([x(1)-d/2 a*diff(y)+x(1)-d/2 y]);
                else
                    d = 1/a*diff(x) - diff(y);
                    
                    axis([x y(1)-d/2 1/a*diff(x)+y(1)-d/2]);
                end
                

            end
        end
        
    end
    
    methods(Static)
        function h1 = find_figure_handles(line_name)
            fg = get(0,'children');
            h1 = [];

            for i = 1:length(fg)
                if strcmp(fg(i).Name, ['Information for line ' line_name])
                    h1 = fg(i);
                end
            end
        end
        
        function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
            % taken from: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
            %INTERSECTIONS Intersections of curves.
            %   Computes the (x,y) locations where two curves intersect.  The curves
            %   can be broken with NaNs or have vertical segments.
            %
            % Example:
            %   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
            %
            % where X1 and Y1 are equal-length vectors of at least two points and
            % represent curve 1.  Similarly, X2 and Y2 represent curve 2.
            % X0 and Y0 are column vectors containing the points at which the two
            % curves intersect.
            %
            % ROBUST (optional) set to 1 or true means to use a slight variation of the
            % algorithm that might return duplicates of some intersection points, and
            % then remove those duplicates.  The default is true, but since the
            % algorithm is slightly slower you can set it to false if you know that
            % your curves don't intersect at any segment boundaries.  Also, the robust
            % version properly handles parallel and overlapping segments.
            %
            % The algorithm can return two additional vectors that indicate which
            % segment pairs contain intersections and where they are:
            %
            %   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
            %
            % For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
            % (how far along this segment the intersection is).  For example, if I(k) =
            % 45.25 then the intersection lies a quarter of the way between the line
            % segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
            % the vector J and the segments in (X2,Y2).
            %
            % You can also get intersections of a curve with itself.  Simply pass in
            % only one curve, i.e.,
            %
            %   [X0,Y0] = intersections(X1,Y1,ROBUST);
            %
            % where, as before, ROBUST is optional.

            % Version: 2.0, 25 May 2017
            % Author:  Douglas M. Schwarz
            % Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
            % Real_email = regexprep(Email,{'=','*'},{'@','.'})


            % Theory of operation:
            %
            % Given two line segments, L1 and L2,
            %
            %   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
            %   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
            %
            % we can write four equations with four unknowns and then solve them.  The
            % four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
            % L1 and L2, t1 is the distance from the starting point of L1 to the
            % intersection relative to the length of L1 and t2 is the distance from the
            % starting point of L2 to the intersection relative to the length of L2.
            %
            % So, the four equations are
            %
            %    (x1(2) - x1(1))*t1 = x0 - x1(1)
            %    (x2(2) - x2(1))*t2 = x0 - x2(1)
            %    (y1(2) - y1(1))*t1 = y0 - y1(1)
            %    (y2(2) - y2(1))*t2 = y0 - y2(1)
            %
            % Rearranging and writing in matrix form,
            %
            %  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
            %        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
            %   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
            %        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
            %
            % Let's call that A*T = B.  We can solve for T with T = A\B.
            %
            % Once we have our solution we just have to look at t1 and t2 to determine
            % whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
            % line segments cross and we can include (x0,y0) in the output.
            %
            % In principle, we have to perform this computation on every pair of line
            % segments in the input data.  This can be quite a large number of pairs so
            % we will reduce it by doing a simple preliminary check to eliminate line
            % segment pairs that could not possibly cross.  The check is to look at the
            % smallest enclosing rectangles (with sides parallel to the axes) for each
            % line segment pair and see if they overlap.  If they do then we have to
            % compute t1 and t2 (via the A\B computation) to see if the line segments
            % cross, but if they don't then the line segments cannot cross.  In a
            % typical application, this technique will eliminate most of the potential
            % line segment pairs.


            % Input checks.
            if verLessThan('matlab','7.13')
                error(nargchk(2,5,nargin)) %#ok<NCHKN>
            else
                narginchk(2,5)
            end

            % Adjustments based on number of arguments.
            switch nargin
                case 2
                    robust = true;
                    x2 = x1;
                    y2 = y1;
                    self_intersect = true;
                case 3
                    robust = x2;
                    x2 = x1;
                    y2 = y1;
                    self_intersect = true;
                case 4
                    robust = true;
                    self_intersect = false;
                case 5
                    self_intersect = false;
            end

            % x1 and y1 must be vectors with same number of points (at least 2).
            if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
                    length(x1) ~= length(y1)
                error('X1 and Y1 must be equal-length vectors of at least 2 points.')
            end
            % x2 and y2 must be vectors with same number of points (at least 2).
            if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
                    length(x2) ~= length(y2)
                error('X2 and Y2 must be equal-length vectors of at least 2 points.')
            end


            % Force all inputs to be column vectors.
            x1 = x1(:);
            y1 = y1(:);
            x2 = x2(:);
            y2 = y2(:);

            % Compute number of line segments in each curve and some differences we'll
            % need later.
            n1 = length(x1) - 1;
            n2 = length(x2) - 1;
            xy1 = [x1 y1];
            xy2 = [x2 y2];
            dxy1 = diff(xy1);
            dxy2 = diff(xy2);


            % Determine the combinations of i and j where the rectangle enclosing the
            % i'th line segment of curve 1 overlaps with the rectangle enclosing the
            % j'th line segment of curve 2.

            % Original method that works in old MATLAB versions, but is slower than
            % using binary singleton expansion (explicit or implicit).
            % [i,j] = find( ...
            % 	repmat(mvmin(x1),1,n2) <= repmat(mvmax(x2).',n1,1) & ...
            % 	repmat(mvmax(x1),1,n2) >= repmat(mvmin(x2).',n1,1) & ...
            % 	repmat(mvmin(y1),1,n2) <= repmat(mvmax(y2).',n1,1) & ...
            % 	repmat(mvmax(y1),1,n2) >= repmat(mvmin(y2).',n1,1));

            % Select an algorithm based on MATLAB version and number of line
            % segments in each curve.  We want to avoid forming large matrices for
            % large numbers of line segments.  If the matrices are not too large,
            % choose the best method available for the MATLAB version.
            if n1 > 1000 || n2 > 1000 || verLessThan('matlab','7.4')
                % Determine which curve has the most line segments.
                if n1 >= n2
                    % Curve 1 has more segments, loop over segments of curve 2.
                    ijc = cell(1,n2);
                    min_x1 = CssVisualization.mvmin(x1);
                    max_x1 = CssVisualization.mvmax(x1);
                    min_y1 = CssVisualization.mvmin(y1);
                    max_y1 = CssVisualization.mvmax(y1);
                    for k = 1:n2
                        k1 = k + 1;
                        ijc{k} = find( ...
                            min_x1 <= max(x2(k),x2(k1)) & max_x1 >= min(x2(k),x2(k1)) & ...
                            min_y1 <= max(y2(k),y2(k1)) & max_y1 >= min(y2(k),y2(k1)));
                        ijc{k}(:,2) = k;
                    end
                    ij = vertcat(ijc{:});
                    i = ij(:,1);
                    j = ij(:,2);
                else
                    % Curve 2 has more segments, loop over segments of curve 1.
                    ijc = cell(1,n1);
                    min_x2 = CssVisualization.mvmin(x2);
                    max_x2 = CssVisualization.mvmax(x2);
                    min_y2 = CssVisualization.mvmin(y2);
                    max_y2 = CssVisualization.mvmax(y2);
                    for k = 1:n1
                        k1 = k + 1;
                        ijc{k}(:,2) = find( ...
                            min_x2 <= max(x1(k),x1(k1)) & max_x2 >= min(x1(k),x1(k1)) & ...
                            min_y2 <= max(y1(k),y1(k1)) & max_y2 >= min(y1(k),y1(k1)));
                        ijc{k}(:,1) = k;
                    end
                    ij = vertcat(ijc{:});
                    i = ij(:,1);
                    j = ij(:,2);
                end

            elseif verLessThan('matlab','9.1')
                % Use bsxfun.
                [i,j] = find( ...
                    bsxfun(@le,CssVisualization.mvmin(x1),CssVisualization.mvmax(x2).') & ...
                    bsxfun(@ge,CssVisualization.mvmax(x1),CssVisualization.mvmin(x2).') & ...
                    bsxfun(@le,CssVisualization.mvmin(y1),CssVisualization.mvmax(y2).') & ...
                    bsxfun(@ge,CssVisualization.mvmax(y1),CssVisualization.mvmin(y2).'));

            else
                % Use implicit expansion.
                [i,j] = find( ...
                    CssVisualization.mvmin(x1) <= CssVisualization.mvmax(x2).' & CssVisualization.mvmax(x1) >= CssVisualization.mvmin(x2).' & ...
                    CssVisualization.mvmin(y1) <= CssVisualization.mvmax(y2).' & CssVisualization.mvmax(y1) >= CssVisualization.mvmin(y2).');

            end


            % Find segments pairs which have at least one vertex = NaN and remove them.
            % This line is a fast way of finding such segment pairs.  We take
            % advantage of the fact that NaNs propagate through calculations, in
            % particular subtraction (in the calculation of dxy1 and dxy2, which we
            % need anyway) and addition.
            % At the same time we can remove redundant combinations of i and j in the
            % case of finding intersections of a line with itself.
            if self_intersect
                remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
            else
                remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
            end
            i(remove) = [];
            j(remove) = [];

            % Initialize matrices.  We'll put the T's and B's in matrices and use them
            % one column at a time.  AA is a 3-D extension of A where we'll use one
            % plane at a time.
            n = length(i);
            T = zeros(4,n);
            AA = zeros(4,4,n);
            AA([1 2],3,:) = -1;
            AA([3 4],4,:) = -1;
            AA([1 3],1,:) = dxy1(i,:).';
            AA([2 4],2,:) = dxy2(j,:).';
            B = -[x1(i) x2(j) y1(i) y2(j)].';

            % Loop through possibilities.  Trap singularity warning and then use
            % lastwarn to see if that plane of AA is near singular.  Process any such
            % segment pairs to determine if they are colinear (overlap) or merely
            % parallel.  That test consists of checking to see if one of the endpoints
            % of the curve 2 segment lies on the curve 1 segment.  This is done by
            % checking the cross product
            %
            %   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
            %
            % If this is close to zero then the segments overlap.

            % If the robust option is false then we assume no two segment pairs are
            % parallel and just go ahead and do the computation.  If A is ever singular
            % a warning will appear.  This is faster and obviously you should use it
            % only when you know you will never have overlapping or parallel segment
            % pairs.

            if robust
                overlap = false(n,1);
                warning_state = warning('off','MATLAB:singularMatrix');
                % Use try-catch to guarantee original warning state is restored.
                try
                    lastwarn('')
                    for k = 1:n
                        T(:,k) = AA(:,:,k)\B(:,k);
                        [unused,last_warn] = lastwarn; %#ok<ASGLU>
                        lastwarn('')
                        if strcmp(last_warn,'MATLAB:singularMatrix')
                            % Force in_range(k) to be false.
                            T(1,k) = NaN;
                            % Determine if these segments overlap or are just parallel.
                            overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
                        end
                    end
                    warning(warning_state)
                catch err
                    warning(warning_state)
                    rethrow(err)
                end
                % Find where t1 and t2 are between 0 and 1 and return the corresponding
                % x0 and y0 values.
                in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
                % For overlapping segment pairs the algorithm will return an
                % intersection point that is at the center of the overlapping region.
                if any(overlap)
                    ia = i(overlap);
                    ja = j(overlap);
                    % set x0 and y0 to middle of overlapping region.
                    T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
                        min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
                    T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
                        min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
                    selected = in_range | overlap;
                else
                    selected = in_range;
                end
                xy0 = T(3:4,selected).';

                % Remove duplicate intersection points.
                [xy0,index] = unique(xy0,'rows');
                x0 = xy0(:,1);
                y0 = xy0(:,2);

                % Compute how far along each line segment the intersections are.
                if nargout > 2
                    sel_index = find(selected);
                    sel = sel_index(index);
                    iout = i(sel) + T(1,sel).';
                    jout = j(sel) + T(2,sel).';
                end
            else % non-robust option
                for k = 1:n
                    [L,U] = lu(AA(:,:,k));
                    T(:,k) = U\(L\B(:,k));
                end

                % Find where t1 and t2 are between 0 and 1 and return the corresponding
                % x0 and y0 values.
                in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
                x0 = T(3,in_range).';
                y0 = T(4,in_range).';

                % Compute how far along each line segment the intersections are.
                if nargout > 2
                    iout = i(in_range) + T(1,in_range).';
                    jout = j(in_range) + T(2,in_range).';
                end
            end

            % Plot the results (useful for debugging).
            % plot(x1,y1,x2,y2,x0,y0,'ok');
        end

        function y = mvmin(x)
            % Faster implementation of movmin(x,k) when k = 1.
            y = min(x(1:end-1),x(2:end));
        end

        function y = mvmax(x)
            % Faster implementation of movmax(x,k) when k = 1.
            y = max(x(1:end-1),x(2:end));
        end
    end
end

