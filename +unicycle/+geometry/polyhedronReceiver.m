classdef polyhedronReceiver < unicycle.geometry.polyhedron
    properties
        % degrees of freedom (number of parameters solved in numerical integration)
        dgf;
        
        % strain
        e11,e12,e13,e22,e23,e33;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = polyhedronReceiver(varargin)
            % POLYHEDRONRECEIVER is a class representing the geometry and 
            % physical properties of receiver polyhedral volume elements.
            %
            %   src = geometry.polyhedronReceiver('basename',earthModel)
            %
            % or
            %
            %   src = geometry.polyhedronReceiver({'basename1','basename2'},earthModel)
            %
            % where 'basename' is short for the triplet 'basename.ned', 
            % 'basename.tri', and 'basename.phn', with
            %
            %   'basename.ned',     the list of vertex coordinates
            %   'basename.tri',     the list of triangle faces
            %   'basename.phn',     the list of faces
            %
            % and earthModel is an object that provides the facilities to
            % computes the stress and displacement field caused by strain
            % in the volume element.
            %
            % SEE ALSO: unicycle, unicycle.geometry.receiver
            
            import unicycle.geometry.polyhedron;
            
            if isempty(varargin)
                return
            else
                basename=varargin{1};
                obj.earthModel=varargin{2};
            end
            
            if ~iscell(basename)
                basename={basename};
            end
            
            obj.x=[];
            obj.vertices=[];
            obj.faces=[];
            for k=1:length(basename)
                
                % list of vertex coordinates
                fname=[basename{k} '.ned'];
                assert(2==exist(fname,'file'),sprintf('error: can''t find %s',fname));
                [~,x1,x2,x3]=...
                    textread(fname,'%u %f %f %f','commentstyle','shell');
                obj.x=[obj.x;[x2,x1,-x3]];
                assert(0>=max(obj.x(:,3)),'error: all vertices should have positive depth.');
                
                % list of triangle faces
                fname=[basename{k} '.tri'];
                assert(2==exist(fname,'file'),sprintf('error: can''t find %s',fname));
                tri=obj.loadTri(fname);
                obj.vertices=[obj.vertices;tri];
                
                fname=[basename{k} '.phn'];
                assert(2==exist(fname,'file'),sprintf('error: can''t find %s',fname));
                faces=obj.loadVolume(fname);
                obj.faces=[obj.faces;faces];
                
            end
            
            % patch properties
            obj.N=size(obj.faces,1);
            obj.id=1:obj.N;
            
            % center coordinates
            obj.xc=zeros(obj.N,3);
            for k=1:length(obj.faces)
                obj.xc(k,:)=sum((obj.x(obj.vertices(obj.faces{k},1),:)+obj.x(obj.vertices(obj.faces{k},2),:)+obj.x(obj.vertices(obj.faces{k},3),:))/3,1)/numel(obj.faces{k});
            end
            
            obj.convex=true;
            
            % unit vectors (stress in north-east-depth coordinate system)
            obj.sv=[zeros(obj.N,1),ones(obj.N,1),zeros(obj.N,1)];
            obj.nv=[ones(obj.N,1),zeros(obj.N,1),zeros(obj.N,1)];
            obj.dv=[zeros(obj.N,1),zeros(obj.N,1),ones(obj.N,1)];

        end % constructor
        
    end % methods
    
end
