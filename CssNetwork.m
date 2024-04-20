classdef CssNetwork
    % CSSNODES obtains the nodes (intersection of lines) from an array of 
    % CssGravityLines
    
    properties
        nodes
        names
        design
        adjacency
        network
    end
    
    methods
        function self = CssNetwork(GravityLines, benchmarks)
            % get a list of the benchmark names
            self.names = sort({benchmarks.name});
            
            % create the adjancency matrix 
            self.adjacency = zeros(length(self.names), length(self.names));
            
            for i = 1:length(GravityLines)
                for j = 1:length(GravityLines(i).instruments)
                    
                    % design matrix for this instrument
                    tdesign = GravityLines(i).design{j};
                    
                    if and(~isempty(tdesign),~any(isnan(GravityLines(i).deltas{j})))
                        % if design is not empty, then there are deltas to
                        % put into the adjustment
                        
                        ben = GravityLines(i).benchmarks;
                        
                        % transpose the design matrix to match the find
                        % sequence (along rows).
                        for k = 1:size(tdesign,1)
                            [~, start_benchmarks] = find(tdesign(k,:) ==  1);
                            [~,   end_benchmarks] = find(tdesign(k,:) == -1);

                            index_st = ismember(self.names, {ben(start_benchmarks).name});
                            index_en = ismember(self.names, {ben(end_benchmarks).name});

                            % assign the values. Make matrix symmetric
                            self.adjacency(index_st,index_en) = self.adjacency(index_st,index_en) + 1;
                            self.adjacency(index_en,index_st) = self.adjacency(index_st,index_en);
                        end
                        
                    end
                end
            end
            
            % calculate the relevant information regarding the network
            self.network = graph(self.adjacency, self.names);
            
            % obtain the nodes (all benchmarks with more than 2 degrees)
            nodes_index = degree(self.network) >= 3;
            self.nodes = CssBenchmark.ReturnBenchmark(benchmarks, self.names(nodes_index));
            
        end
        
        function orphans = GetOrphans(self, abs_benchmarks)
            % Obtain a list of the benchmarks that are not connected to an
            % absolute benchmark. 
            % input:
            %    abs_benchmarks -> list of CssBenchmark
            % the orphan list is a function of the connection to absolute
            % benchmarks. Use the absolute benchmarks
            % obtain orphans (lines/points not connected to the network)
            
            bins = conncomp(self.network, 'OutputForm', 'cell');
            
            % FSS: 10/26/2021 the empty array was defined as [] but it
            % should be an empty cell array {}. Otherwise can't compare
            % array using ismember
            orphans = {};
            
            if length(bins) > 1
                abs_names = {abs_benchmarks.name};

                for i = 1:length(bins)
                    bnames = bins{i};
                    if ~any(ismember(abs_names, bnames))
                        orphans = [orphans; bins{i}'];
                    end
                end
            end
        end
        
        function PlotNetwork(self, abs_benchmarks)
            
            % get the names of the absolute benchmarks
            abs_names = {abs_benchmarks.name};
            
            p = plot(self.network, 'NodeLabel', self.network.Nodes.Name, 'Layout','layered', 'sources', abs_names, 'direction', 'right');
            
            % paint abs benchmarks green and make them larger
            highlight(p, abs_names)
            highlight(p, abs_names, 'NodeColor', 'g');
            
            % make the edges proportional to the number of observations
            self.network.Edges.LWidths = 7*self.network.Edges.Weight/max(self.network.Edges.Weight);
            
            p.LineWidth = self.network.Edges.LWidths;
        end
    end
    
end

