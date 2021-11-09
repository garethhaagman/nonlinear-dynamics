%Function to perform arc-length continuation for nonlinear systems
function ylist=MyTrackCurve(userf,userdf,y0,ytan)      
nmax = 2500;     %Maximum number of steps
s = -0.1;        %Step size
tol = 1e-6;             %Tolerance
maxit = 200;            %Maximum iterations
h = 0.01;               %Finite difference parameter 
n = numel(y0);
userdf = MyJacobian(userf,y0,h); 
[Q,~] = qr(userdf');                    %QR decomposition
ytan  = Q(:,end);                    %Xdot for first point
ylist = zeros(n,nmax);
for j = 1:nmax
    ylist(:,j) = y0;
    ypred = y0 + s*ytan;                
    y0 = ypred;                          
    for k = 1:maxit                 
        Fx = userf(y0);
        userdf = [MyJacobian(userf,y0,h); ytan.'];
        if norm(Fx) <= tol, break; end
        y0 = y0 -userdf\[Fx; dot((y0-ypred),ytan)];
    end
    ytan = userdf\[zeros(n-1,1);1];
end
