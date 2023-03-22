function u0 = oscic(x)  
    %mu = 0; sig = 0.25; u0 = 1/sqrt(2*pi*sig^2)*exp(-((x-mu)^2/sig^2)/2); % IC gaussiana
    u0 = sin(pi*x);
    
end

