function FR = FisherR(Y, A, HR, sigma2)
%FISHER 此处显示有关此函数的摘要
%   此处显示详细说明
FR  = zeros(3, 3);
M   = size(Y, 2);
for i = 1:M
    FR(1, 1)    = FR(1, 1) + (-1)/sigma2*(-Y(:,i)'*A(:,:,i)*HR.dChdt2(:,i) - ...
                (HR.dChdt2(:,i))'*A(:,:,i)'*Y(:,i) + (HR.dChdt2(:,i))'*(A(:,:,i)'*A(:,:,i))*HR.Ch(:,i) + ...
                2*(HR.dChdt(:,i))'*(A(:,:,i)'*A(:,:,i))*(HR.dChdt(:,i)) + HR.Ch(:,i)'*(A(:,:,i)'*A(:,:,i))*(HR.dChdt2(:,i)));

    FR(1, 2)    = FR(1, 2) + (-1)/sigma2*(-Y(:,i)'*A(:,:,i)*HR.dChdtdr(:,i) - ...
                (HR.dChdtdr(:,i))'*A(:,:,i)'*Y(:,i) + (HR.dChdtdr(:,i))'*(A(:,:,i)'*A(:,:,i))*HR.Ch(:,i) + ...
                2*real((HR.dChdt(:,i))'*(A(:,:,i)'*A(:,:,i))*(HR.dChdr(:,i))) + HR.Ch(:,i)'*(A(:,:,i)'*A(:,:,i))*(HR.dChdtdr(:,i)));

    FR(1, 3)    = FR(1, 3) + (-1)/sigma2*(-Y(:,i)'*A(:,:,i)*HR.dChdtda(:,i) - ...
                (HR.dChdtda(:,i))'*A(:,:,i)'*Y(:,i) + (HR.dChdtda(:,i))'*(A(:,:,i)'*A(:,:,i))*HR.Ch(:,i) + ...
                2*real((HR.dChdt(:,i))'*(A(:,:,i)'*A(:,:,i))*(HR.dChda(:,i))) + HR.Ch(:,i)'*(A(:,:,i)'*A(:,:,i))*(HR.dChdtda(:,i)));    

    FR(2, 2)    = FR(2, 2) + (-1)/sigma2*(-Y(:,i)'*A(:,:,i)*HR.dChdr2(:,i) - ...
                (HR.dChdr2(:,i))'*A(:,:,i)'*Y(:,i) + (HR.dChdr2(:,i))'*(A(:,:,i)'*A(:,:,i))*HR.Ch(:,i) + ...
                2*(HR.dChdr(:,i))'*(A(:,:,i)'*A(:,:,i))*(HR.dChdr(:,i)) + HR.Ch(:,i)'*(A(:,:,i)'*A(:,:,i))*(HR.dChdr2(:,i)));

    FR(2, 3)    = FR(2, 3) + (-1)/sigma2*(-Y(:,i)'*A(:,:,i)*HR.dChdrda(:,i) - ...
                (HR.dChdrda(:,i))'*A(:,:,i)'*Y(:,i) + (HR.dChdrda(:,i))'*(A(:,:,i)'*A(:,:,i))*HR.Ch(:,i) + ...
                2*real((HR.dChdr(:,i))'*(A(:,:,i)'*A(:,:,i))*(HR.dChda(:,i))) + HR.Ch(:,i)'*(A(:,:,i)'*A(:,:,i))*(HR.dChdrda(:,i)));

    FR(3, 3)    = FR(3, 3) + (-1)/sigma2*(-Y(:,i)'*A(:,:,i)*HR.dChda2(:,i) - ...
                (HR.dChda2(:,i))'*A(:,:,i)'*Y(:,i) + (HR.dChda2(:,i))'*(A(:,:,i)'*A(:,:,i))*HR.Ch(:,i) + ...
                2*(HR.dChda(:,i))'*(A(:,:,i)'*A(:,:,i))*(HR.dChda(:,i)) + HR.Ch(:,i)'*(A(:,:,i)'*A(:,:,i))*(HR.dChda2(:,i)));
end
FR(2, 1)    = FR(1, 2);
FR(3, 1)    = FR(1, 3);
FR(3, 2)    = FR(2, 3);
FR  = -FR;
end

