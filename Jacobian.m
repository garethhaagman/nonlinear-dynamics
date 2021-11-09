%Function to Compute the Jacobian of an arbitrary function f at a given point x
function df = MyJacobian(f,x,h)
n = numel(x);
m = zeros(n,1);
%Recursion of Approximation
for i=1:n
    m(i) = h; 
    df(:,i) = (f(x+m)-f(x-m))/(2*h);
    m(i) = 0; %Zeros
end
end
