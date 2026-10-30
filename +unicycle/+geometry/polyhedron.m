classdef polyhedron < handle
    properties
        % number of polyhedra
        N;
        % identification number
        id;
        % edge coordinates
        x;
        % triangle faces
        vertices;
        % polyhedron
        faces;
        % is convex?
        convex;
        % center position
        xc;
        % unit vectors in strike, dip and normal directions
        sv;
        dv;
        nv;
        % volume of polyhedron
        volume;
        % earth model
        earthModel;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = polyhedron()
            % POLYHEDRON is a meta class representing the geometry of 
            % polyhedron volume elements, defined in terms of the vertices
            % and faces, as illustrated below:
            %
            %                                + E
            %                 N (x1)    .   ..
            %                /     .       . .
            %               / .           .C .
            %   x1,x2,x3 ->A ----------- +   .
            %              |\           /  \ .
            %              : \         /    \.
            %              |  \       /      + D
            %              :   \     /     /
            %              |    \   /   /
            %              :     \ / /
            %              |      +
            %              :       B
            %              Z (x3)
            %
            % SEE ALSO: unicycle, unicycle.geometry.polyhedronReceiver
            
            if (0==nargin)
                return
            end
            
        end % constructor
        
        function [varargout]=tractionKernels(obj,rcv,varargin)
            % TRACTIONKERNELS computes the traction on receivers faults
            % due to strain in volume elements. 
            %
            % rcv - receiver object
            
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.tractionKernels(obj,rcv,varargin{:});
        end
        
        function [varargout]=stressKernels(obj,rcv,varargin)
            % STRESSKERNELS computes the stress on receivers due to 
            % strain in volume elements. 
            %
            % rcv - receiver object
            
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.stressKernels(obj,rcv,varargin{:});
        end
        
        function [varargout]=displacementKernels(obj,x,vecsize,varargin)
            % DISPLACEMENTKERNELS computes the displacement field at x
            % due to strain in volume elements.
            %
            % x       - observations points coordinates
            % vecsize - length of displacement vector
       
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.displacementKernels(obj,x,vecsize,varargin{:});
        end
        
        function plotVolume(obj,varargin)
            % PLOTVOLUME plot contour of volume elements in 3d
            %
            % SEE ALSO: unicycle
            
            if 1>obj.N
                fprintf('unicycle.geometry.polyhedron: nothing to plot\n');
                return
            end
            if nargin > 1
                scatter3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),150,varargin{1},'filled');
                for k=1:obj.N
                    trisurf(obj.vertices(obj.faces{k},:),obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                        varargin{1}(k),'FaceAlpha',0.3), shading flat;
                end
            else
                hold on
                trimesh(obj.vertices(:,[1,2,3]),obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                    0.*obj.x(:,3),'FaceColor','None','LineWidth',1);
            end
        end
        
        function plotVolumeById(obj,Id,varargin)
            % PLOTVOLUMEBYID plot contour of a single volume element in 3d
            %
            % SEE ALSO: unicycle
            
            if 1>obj.N
                fprintf('unicycle.geometry.polyhedron: nothing to plot\n');
                return
            end
            
            if nargin > 2
                trisurf(obj.vertices(obj.faces{Id},:),obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                    varargin{1}), shading flat;
            else
                plot3(obj.xc(Id,1),obj.xc(Id,2),obj.xc(Id,3),'k+');
                
                trisurf(obj.vertices(obj.faces{Id},:),obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                    0.*obj.x(:,3),'FaceColor','None','LineWidth',1);
            end
        end
        
        function [xp,yp,zp,dim]=computeVertexPosition(obj)
            % COMPUTEVERTEXPOSITION computes the position of vertices
            %
            % SEE ALSO: unicycle
            
            xp=[obj.x(obj.vertices(:,1),1),obj.x(obj.vertices(:,2),1),obj.x(obj.vertices(:,3),1),obj.x(obj.vertices(:,4),1)]';
            yp=[obj.x(obj.vertices(:,1),2),obj.x(obj.vertices(:,2),2),obj.x(obj.vertices(:,3),2),obj.x(obj.vertices(:,4),2)]';
            zp=[obj.x(obj.vertices(:,1),3),obj.x(obj.vertices(:,2),3),obj.x(obj.vertices(:,3),3),obj.x(obj.vertices(:,4),3)]';
            dim=3;
        end
        
        function plotUnitVectors(obj,sc)
            % PLOTUNITVECTORS plot normal, strike-slip and dip-slip unit
            % vectors.
            %
            % SEE ALSO: unicycle
            
            hold on
            
            for k=1:size(obj.faces,1)
                
                triangles=obj.vertices(obj.faces{k},:);
                
                % center coordinates
                O=sum((obj.x(triangles(:,1),:) ...
                      +obj.x(triangles(:,2),:) ...
                      +obj.x(triangles(:,3),:))/3,1)/size(triangles,1);
                
                % loop of triangle elements
                for i=1:size(triangles,1)
                
                    A=obj.x(triangles(i,1),:);
                    B=obj.x(triangles(i,2),:);
                    C=obj.x(triangles(i,3),:);
                    
                    % center of face
                    H=(A+B+C)/3;
                    
                    % unit normal vectors
                    n=cross(B-A,C-A); % ABC
                    n=n/norm(n);
                    
                    if obj.convex
                        % check that unit vectors are pointing outward
                        if (n'*(O(:)-(A(:)+B(:)+C(:))/3))>0
                            n=-n;
                        end
                    else
                        % viewed from outside, vertices are given in clockwise order
                        n=-n;
                    end
            
                    plot3(A(1),A(2),A(3),'k+');
                    plot3(B(1),B(2),B(3),'k+');
                    plot3(C(1),C(2),C(3),'k+');
                    
                    quiver3(H(1),H(2),H(3),sc*n(1),sc*n(2),sc*n(3),0,'r');
                end
            end
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                            %
        %        export geometry to GMT .xyz format                  %
        %                                                            %
        % INPUT:                                                     %
        %                                                            % 
        % scale        - scaling factor of geometrical features      %
        % fname        - output file name                            %
        % value        - attribute to plot                           % 
        %                                                            %
        % EXAMPLE:                                                   %
        %                                                            %
        %   shz.exportXYZ(1e-3,'output/tetrahedra.xyz',value)   %
        %                                                            %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportXYZ(o,scale,fname,field)
            if 1>o.N
                return
            end
            
            % vertices
            A=[o.x(o.vertices(:,1),1),o.x(o.vertices(:,1),2),o.x(o.vertices(:,1),3)]*scale;
            B=[o.x(o.vertices(:,2),1),o.x(o.vertices(:,2),2),o.x(o.vertices(:,2),3)]*scale;
            C=[o.x(o.vertices(:,3),1),o.x(o.vertices(:,3),2),o.x(o.vertices(:,3),3)]*scale;
            D=[o.x(o.vertices(:,4),1),o.x(o.vertices(:,4),2),o.x(o.vertices(:,4),3)]*scale;

            fid=fopen(fname,'wt');
            for k=1:length(field)
                if (~isnan(field(k)))
                    fprintf(fid,'> -Z%f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n', ...
                        [field(k);A(k,:)';B(k,:)';C(k,:)';A(k,:)';D(k,:)';C(k,:)';D(k,:)';B(k,:)']);
                end
            end
            fclose(fid);
            
        end
        
        
        
    end % methods
    
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % %
        %                                            %
        %  compute volume from position of vertices  %
        %                                            %
        % % % % % % % % % % % % % % % % % % % % % % %%
        function volume = computeVolume(x,vertices)
            % COMPUTEVOLUME compute the polyhedra volumes.
            %
            %   volume=polyhedron.computeVolume(x,vertices)
            %
            % INPUT:
            %   x        - position of mesh points
            %   vertices - list of vertices forming polyhedra
            %
            % SEE ALSO: unicycle, unicycle.geometry.polyhedron
            
            % vertices
            A=[x(vertices(:,1),1),x(vertices(:,1),2),x(vertices(:,1),3)];
            B=[x(vertices(:,2),1),x(vertices(:,2),2),x(vertices(:,2),3)];
            C=[x(vertices(:,3),1),x(vertices(:,3),2),x(vertices(:,3),3)];
            D=[x(vertices(:,4),1),x(vertices(:,4),2),x(vertices(:,4),3)];
            
            % volume |(D-A) . (B-A) x (C-A)| / 6
            volume=abs(sum((D(:,1)-A(:,1)).*(B(:,2)-A(:,2)).*(C(:,3)-A(:,3))-(B(:,3)-A(:,3)).*(C(:,2)-A(:,2)) ...
                          +(D(:,2)-A(:,2)).*(B(:,3)-A(:,3)).*(C(:,1)-A(:,1))-(B(:,1)-A(:,1)).*(C(:,3)-A(:,3)) ...
                          +(D(:,3)-A(:,3)).*(B(:,1)-A(:,1)).*(C(:,2)-A(:,2))-(B(:,2)-A(:,2)).*(C(:,1)-A(:,1)),2)/6);
            
        end
        
        function tri = loadTri(filename)
            % LOADTRI loads the .tri mesh file
            %
            %   tri = loadTri(filename)
            %
            % where tri is a list of mesh vertices.
            
            [~,i1,i2,i3]=textread(filename,'%u %d %d %d','commentstyle','shell');
            tri = [i1(:),i2(:),i3(:)];
            
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %        load the triangle mesh from file          %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function faces=loadVolume(filename)
            % LOADVOLUME loads the .phn mesh file.
            %
            %   phn.loadVolume(filename)
            %
            % where mesh is a list of eij and mesh vertices, for example:
            %
            % # index n i1 ... in
            %
            
            % count the number of polyhedra
            fid=fopen(filename);
            N=0;
            while 1
                line=fgetl(fid);
                
                % test end of file
                if ~ischar(line), break, end
                
                % trim white space
                line=strtrim(line);
                
                % remove comments
                if (strcmp('#',line(1))), continue, end
                
                % count valid input
                N=N+1;
            end
            fclose(fid);
            
            % initialize number of polyhedra
            faces=cell(N,1);
            
            % read input polyhedron data
            fid=fopen(filename);
            k=1;
            while 1
                line=fgetl(fid);
                % test end of file
                if ~ischar(line), break, end
                % trim white space
                line=strtrim(line);
                % remove comments
                if (strcmp('#',line(1))), continue, end
                
                % read second column
                numberOfFaces=cell2mat(textscan(line,'%*d %d %*[^\n]'));
                % read n faces
                faces{k}=cell2mat(textscan(line,['%*d %*d ' repmat('%d ',[1,numberOfFaces]) ' %*[^\n]']))';
                k=k+1;
            end
            fclose(fid);
        end
        
    end % methods (Static)

end
