function [ R ] = rMatrix( alpha , theta )

%Function to calculate the rotational matrices

R(:,:) = [cos(theta) , -sin(theta)*cos(alpha) , sin(theta)*sin(alpha);
    sin(theta) , cos(theta)*cos(alpha) , - (cos(theta)*sin(alpha)) ;
    0 , sin(alpha) , cos(alpha)] ;


end

