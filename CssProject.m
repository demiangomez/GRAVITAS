classdef CssProject
    %CSSPROJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        working_folder
    end
    
    methods
        function self = CssProject(path)
            % try to open the project in 'path'
            self.working_folder = path;
            
            [result, message] = self.CheckFolderStructure();
            
            if ~result
                error(message)
            end
        end
        
        function [result, msg] = CheckFolderStructure(self)
            % this function check the existance of all the necessary folders and
            % structures.

            result = false;

            if ~exist(fullfile(self.working_folder, 'lines'),'dir')
                msg = 'Failed checking lines folder.';
                return
            end

            if ~exist(fullfile(self.working_folder, 'instruments'),'dir')
                msg = 'Failed checking instruments.';
                return
            end

%             if ~exist(fullfile(self.working_folder, 'absolute'),'dir')
%                 msg = 'Failed checking absolute gravity data.';
%                 return
%             end

            if ~exist(fullfile(self.working_folder, 'database.mat'),'file')
                msg = 'Failed checking database mat file.';
                return
            end

            result = true;
            msg   = '';
        end
        
        function out = LoadDatabase(self)
            t = load(fullfile(self.working_folder,'database.mat'));
            out = t.database;
        end
        
        function lines = LoadLines(self)
    
            % load the lines
            % loop through the folder files searching for lines
            files = dir(fullfile(self.working_folder,'lines/*.mat'));
            % remove the directories
            files([files.isdir] == 1) = [];
            files(~cellfun(@isempty,strfind({files.name},'DS_Store'))) = [];
            
            for i = 1:length(files)
                lines(i,1) = CssGravityLine(fullfile(self.working_folder,['lines/' files(i).name]));
            end
            
            if ~exist('lines', 'var')
                lines = [];
            end

            % After a discussion with Kevin (09-05-2017) there seems to be
            % no point in keeping these lines separated from the rest.
%             % load the absolute gravity lines
%             files = dir(fullfile(self.working_folder,'absolute'));
%             % remove the directories
%             files([files.isdir] == 1) = [];
%             n = length(lines);
%             for i = 1:length(files)
%                 load(fullfile(self.working_folder,['absolute/' files(i).name]))
%                 lines(n+i,1) = line_struct;
%             end
        end
        
        function instruments = LoadInstruments(self)
            files = dir(fullfile(self.working_folder,'instruments/*.mat'));
            % remove the directories
            files([files.isdir] == 1) = [];
            files(~cellfun(@isempty,strfind({files.name},'DS_Store'))) = [];
            
            for i = 1:length(files)
                load(fullfile(self.working_folder,['instruments/' files(i).name]))
                instruments(i,1) = instrument_struct;
            end
            
            if ~exist('instruments', 'var')
                instruments = [];
            end
        end
        
        function UpdateBenchmarks(self, nm_benchmarks)
            % function to update the list of benchmarks in the database.mat
            % file and in all lines using the benchmarks passed as
            % nm_benchmarks
            
            database = self.LoadDatabase();
            
            for i = 1:length(nm_benchmarks)
                if ~CssBenchmark.exists(database.benchmarks, nm_benchmarks(i).name)
                    % this is a new benchmarks. It just needs to be added
                    % to the list
                    database.benchmarks = [database.benchmarks; copy(nm_benchmarks(i))];
                else
                    % now explore the lines and replace this benchmark with
                    % the new version
                    lines = self.LoadLines();
                    instruments = self.LoadInstruments();
                    for j = 1:length(lines)
                        if CssBenchmark.exists(lines(j).benchmarks, nm_benchmarks(i).name)
                            % benchmark is member of this line
                            % replace
                            lines(j) = lines(j).UpdateBenchmark(nm_benchmarks(i).name, nm_benchmarks(i), instruments);
                            % save line
                            lines(j).SaveGravityLine(fullfile(self.working_folder, 'lines'), true)
                        end
                    end
                    
                    % update the reference in the database
                    benchmark = CssBenchmark.ReturnBenchmark(database.benchmarks, nm_benchmarks(i).name);
                    benchmark = benchmark.UpdateCoords(nm_benchmarks(i));
                end
            end
            
            % write the database to the disk
            save(fullfile(self.working_folder, 'database.mat'),'database');
        end
    end
    
end

