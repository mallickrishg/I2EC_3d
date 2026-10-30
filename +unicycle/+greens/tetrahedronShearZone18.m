classdef tetrahedronShearZone18 < unicycle.greens.earthModel
    properties
        % rigidity
        G;
        % Poisson's ratio
        nu;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=tetrahedronShearZone18(G,nu)
            % TETRAHEDRONSHEARZONE18 is a class providing the necessary functions to
            % compute the stress interactions between source and receiver
            % tetrahedral shear zones.
            %
            %   earthModel = greens.tetrahedronShearZone18(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space.
            %
            % SEE ALSO: unicycle
            
            if (0==nargin)
                return
            end
            
            assert(0<=G,'unicycle.greens.tetrahedronShearZone18::rigidity must be positive.')
            assert(nu<=0.5,'unicycle.greens.tetrahedronShearZone18::Poisson''s ratio should be lower than 0.5.')
            assert(-1<=nu,'unicycle.greens.tetrahedronShearZone18::Poisson''s ratio should be greater than -1.')
            
            o.G=G;
            o.nu=nu;
        end
        
        function [varargout]=tractionKernels(obj,src,rcv,varargin)
            % TRACTIONKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % rcv - receiver fault
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            varargout = cell(1,nargout);
            [varargout{:}]=unicycle.greens.computeTractionKernelsTetrahedronShearZone(src,rcv,obj.G,obj.nu,varargin{:});
        end
        
        function [varargout]=stressKernels(obj,src,rcv,varargin)
            % STRESSKERNELS computes the stress on receiver faults due to
            % motion of rectangular dislocations in a half space.
            %
            % rcv - receiver shear zone
            %
            % SEE ALSO: unicycle
            
            varargout = cell(1,nargout);
            [varargout{:}]=unicycle.greens.computeStressKernelsTetrahedronShearZone(src,rcv,obj.G,obj.nu,varargin{:});
        end
        
        function [varargout]=displacementKernels(obj,src,x,vecsize,varargin)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of rectangular dislocations in a half space.
            %
            % src - source fault
            %
            % SEE ALSO: unicycle

            varargout = cell(1,nargout);
            [varargout{:}]=unicycle.greens.computeDisplacementKernelsTetrahedronShearZone(src,obj.nu,x,vecsize,varargin{:});
        end
    end % methods
end % class definition
