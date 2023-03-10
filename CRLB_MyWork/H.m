classdef H
    %H JSTSP中的宽带多径信道
    
    properties        
        alpha
        fc
        lambdac
        r
        N
        theta
        d
        M
        B

        Ch
        dChdt
        dChdt2
        dChdr
        dChdr2
        dChda
        dChda2

        dChdtdr
        dChdtda
        dChdrda
    end
    
    methods
        function obj = H(alpha, fc, r, N, theta, M, B)
            %H 构造此类的实例, 有子载波那个维度
            obj.alpha       = alpha;
            obj.fc          = fc;     
            obj.lambdac     = 3e8/fc;
            obj.r           = r;
            obj.N           = N;
            obj.theta       = theta;
            obj.d           = obj.lambdac/2;
            obj.M           = M;
            obj.B           = B;
            
            delta           = ((2*(0:N-1)-N+1)/2).';
            lambda          = 3e8./(fc+((0:M-1)-M/2+1)/M*B);
            %% 求信道本身
            for i = 1:M
                obj.Ch(:,i) = alpha*exp(-1j*2*pi./lambda(i).*...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta));
            end
            %% 求信道对theta的一阶偏导
            for i = 1:M
                obj.dChdt(:,i)   = obj.Ch(:,i).*...
                                (1j*2*pi/lambda(i)).*...
                                r.*delta.*obj.d./...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta);
            end
            %% 求信道对theta的二阶偏导
            for i = 1:M
                obj.dChdt2(:,i)   = obj.Ch(:,i).*(...
                                (-1*4*pi^2/lambda(i)^2).*...
                                (r^2.*delta.*delta.*obj.d^2)./...
                                (r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta)...
                                +...
                                (1j*2*pi/lambda(i)).*...
                                (r^2.*delta.*delta.*obj.d^2)./...
                                (r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta).^(3/2));
            end
            %% 求信道对r的一阶偏导
            for i = 1:M
                obj.dChdr(:,i)   = obj.Ch(:,i).*...
                                (-1j*2*pi/lambda(i)).*...
                                (r - theta.*delta.*obj.d)./...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta);
            end
            %% 求信道对r的二阶偏导
            for i = 1:M
                obj.dChdr2(:,i)   = obj.Ch(:,i).*...
                                (-4*pi^2/lambda(i)^2).*...
                                (r - theta.*delta.*obj.d).^2./...
                                (r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta)...
                                +...
                                obj.Ch(:,i).*...
                                (-1j*2*pi./lambda(i)).*...
                                (sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta)...
                                - (r - theta.*delta.*obj.d)./...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta))./...
                                (r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta);
            end
            %% 求信道对alpha的一阶偏导
            for i = 1:M
                obj.dChda(:,i)   = exp(-1j*2*pi./lambda(i).*...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta));
            end
            %% 求信道对alpha的二阶偏导
            obj.dChda2      = zeros(N, M);
            %% 求信道对theta的一阶偏导r的一阶偏导                 
            for i = 1:M
                obj.dChdtdr(:,i)   = obj.Ch(:,i).*...
                                (4*pi^2/lambda(i)^2).*...
                                (r - theta.*delta.*obj.d).*(r.*delta.*obj.d)./...
                                (r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta)...
                                +...
                                obj.Ch(:,i).*...
                                (1j*2*pi/lambda(i)).*...
                                (delta.*obj.d.*sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta)-...
                                r.*delta.*obj.d.*(r - theta.*delta.*obj.d)./...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta))./...
                                (r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta);
            end
            %% 求信道对theta的一阶偏导a的一阶偏导
            for i = 1:M
                obj.dChdtda(:,i)   = obj.Ch(:,i)./alpha.*...
                                (1j*2*pi/lambda(i)).*...
                                r.*delta.*obj.d./...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta);
            end
            %% 求信道对r的一阶偏导a的一阶偏导
            for i = 1:M
                obj.dChdrda(:,i)   = obj.Ch(:,i)./alpha.*...
                                (-1j*2*pi/lambda(i)).*...
                                (r - theta.*delta.*obj.d)./...
                                sqrt(r^2 + obj.d^2.*delta.^2 - 2*r*theta*obj.d.*delta);
            end
        end
        
    end
end

