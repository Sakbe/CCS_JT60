function [PID] = buildPID(params)

% Filters
s = tf('s');

Kp = params.Kp;
Kd = params.Kd;
Ki = params.Ki;
Td = params.Td; 

P = tf(zeros([length(Kp) length(Kp)]));
for i = 1 : length(Kp)
    I(i, i) = Ki(i)/s;
    D(i, i) = Kd(i) * s * Td(i)/(1 + s * Td(i));
    
    PID(i, i) = Kp(i) + I(i, i) + D(i, i);
end

end