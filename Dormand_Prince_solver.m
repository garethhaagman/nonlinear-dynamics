%% Dormand-Prince (dopri) Method
% Doubling N divides error by a factor of 32
function [xend,t,xt] = MyIVP(f,x0,tspan,N)
    %Calculate integration timestep size from N
    h=(tspan(2)-tspan(1))/N;
    xt(:,1) = x0;
    k = 2;
    t = tspan(1);
    while t < tspan(2)
        yn = xt(:,k-1);
        tn = t(k-1);
        k1 = h * f(tn, yn);
        k2 = h * f(tn+h/5, yn+k1/5);
        k3 = h * f(tn+h*3/10, yn+k1*3/40 + k2*9/40);
        k4 = h * f(tn+h*4/5, yn + k1*44/45 - k2*56/15 + k3*32/9);
        k5 = h * f(tn+h*8/9, yn + k1*19372/6561 - k2*25360/2187 + k3*64448/6561 - k4*212/729);
        k6 = h * f(tn+h, yn + k1*9017/3168 - k2*355/33 + k3*46732/5247 + k4*49/176 - k5*5103/18656);
        k7 = h * f(tn+h, yn + k1*35/384 + k3*500/1113 + k4*125/192 - k5*2187/6784 + k6*11/84);
        k8 = yn + 5179/57600*k1 + 7571/16695*k3 + 393/640*k4 - 92097/339200*k5 + 187/2100*k6 + 1/40*k7;
        
        xt(:,k) = k8;
        t(k) = t(k-1) + h;
        k = k + 1;
    end
    xend = xt(end);
end
