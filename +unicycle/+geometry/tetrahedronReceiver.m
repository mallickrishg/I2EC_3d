classdef tetrahedronReceiver < unicycle.geometry.tetrahedron
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
        function obj = tetrahedronReceiver(varargin)
            % TETRAHEDRONRECEIVER is a class representing the geometry and 
            % physicalproperties of receiver tetrahedral volume elements.
            %
            %   src = geometry.tetrahedronReceiver('basename')
            %
            % or
            %
            %   src = geometry.tetrahedronReceiver({'basename1','basename2'})
            %
            % where 'basename' is short for the pair 'basename.tet' and 'basename.ned'.
            %
            % SEE ALSO: unicycle, unicycle.geometry.receiver
            
            import unicycle.geometry.tetrahedron;
            
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
            for k=1:length(basename)
                fname=[basename{k} '.ned'];
                assert(2==exist(fname,'file'),sprintf('error: can''t find %s',fname));
                [~,x1,x2,x3]=...
                    textread(fname,'%u %f %f %f','commentstyle','shell');
                obj.x=[obj.x;[x2,x1,-x3]];
                
                assert(0>=max(obj.x(:,3)),'error: all vertices should have positive depth.');
                
                fname=[basename{k} '.tet'];
                assert(2==exist(fname,'file'),sprintf('error: can''t find %s',fname));
                    
                [e11,e12,e13,e22,e23,e33,mesh]=obj.loadVolume(fname);
                obj.vertices=[obj.vertices;mesh];
                obj.e11=[obj.e11;e11(:)];
                obj.e12=[obj.e12;e12(:)];
                obj.e13=[obj.e13;e13(:)];
                obj.e22=[obj.e22;e22(:)];
                obj.e23=[obj.e23;e23(:)];
                obj.e33=[obj.e33;e33(:)];
            end
            
            % patch properties
            obj.N=size(obj.vertices,1);
            obj.id=1:obj.N;
            
            % center of fault patch
            obj.xc=[(obj.x(obj.vertices(:,1),1)+ ...
                     obj.x(obj.vertices(:,2),1)+ ...
                     obj.x(obj.vertices(:,3),1)+ ...
                     obj.x(obj.vertices(:,4),1))/4, ...
                    (obj.x(obj.vertices(:,1),2)+ ...
                     obj.x(obj.vertices(:,2),2)+ ...
                     obj.x(obj.vertices(:,3),2)+ ...
                     obj.x(obj.vertices(:,4),2))/4, ...
                    (obj.x(obj.vertices(:,1),3)+ ...
                     obj.x(obj.vertices(:,2),3)+ ...
                     obj.x(obj.vertices(:,3),3)+ ...
                     obj.x(obj.vertices(:,4),3))/4];
                
            % circumsphere radius of tetrahedron
            obj.circumRadius=zeros(obj.N,1);
            
            % minimum radius of triangle faces incircle
            obj.minInscribedCircleRadius=zeros(obj.N,1);
            
            for k=1:obj.N
                A=obj.x(obj.vertices(k,1),:)';
                B=obj.x(obj.vertices(k,2),:)';
                C=obj.x(obj.vertices(k,3),:)';
                D=obj.x(obj.vertices(k,4),:)';
                
                a=det([[A',1];[B',1];[C',1];[D',1]]);
                O(1)=det([[sum(A.^2),A([2,3])',1];[sum(B.^2),B([2,3])',1];[sum(C.^2),C([2,3])',1];[sum(D.^2),D([2,3])',1]])/2/a;
                O(2)=det([[sum(A.^2),A([1,3])',1];[sum(B.^2),B([1,3])',1];[sum(C.^2),C([1,3])',1];[sum(D.^2),D([1,3])',1]])/2/a;
                O(3)=det([[sum(A.^2),A([1,2])',1];[sum(B.^2),B([1,2])',1];[sum(C.^2),C([1,2])',1];[sum(D.^2),D([1,2])',1]])/2/a;
                
                % radius of circumsphere
                obj.circumRadius(k)=norm(O(:)-A(:));
                
                % minimum radius of triangle faces incircle
                obj.minInscribedCircleRadius(k)=min( ...
                    [norm(cross(B-A,C-A))./(norm(B-A)+norm(C-A,2)+norm(B-C,2)), ...
                     norm(cross(B-A,D-A))./(norm(B-A)+norm(D-A,2)+norm(B-D,2)), ...
                     norm(cross(C-A,D-A))./(norm(C-A)+norm(D-A,2)+norm(C-D,2)), ...
                     norm(cross(C-B,D-B))./(norm(C-B)+norm(D-B,2)+norm(C-D,2))]);
                
            end
            
            % volume of tetrahedron
            obj.volume=tetrahedron.computeVolume(obj.x,obj.vertices);
            
            % unit vectors (stress in north-east-depth coordinate system)
            obj.sv=[zeros(obj.N,1),ones(obj.N,1),zeros(obj.N,1)];
            obj.nv=[ones(obj.N,1),zeros(obj.N,1),zeros(obj.N,1)];
            obj.dv=[zeros(obj.N,1),zeros(obj.N,1),ones(obj.N,1)];

        end % constructor
        
    end % methods
    
end
