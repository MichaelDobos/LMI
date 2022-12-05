% Interval Quadratic Stability
% Example: Pitch Channel Model from Chpt. 11 Duan, Yu

tz=1; % given variable

a1=[1.593 .398]; % interval uncertanties given Table 11.2
ap1=.2;
%ap1=[.285 .043];
a2=[260.559 51.003];
a3=[185.488 53.84];
a4=[1.506 .421];
a5=[.298 .078];

n=3; % matrix dimension
tol=1e-8; % inequality tolerance
P=sdpvar(n); % define variables

F=[]; % define system for constraints using interval uncertainties
F=[F P>=tol*eye(n)]; % symmetric positive matrix constraint
for i=1:2 % define A matrix and 2^k interval constraints
    for j=1:2
        for m=1:2
            for u=1:2
                for v=1:2
                    for w=1:2
                        A=[-a4(v) 1 -a5(w) ; -ap1*a4(v)-a2(m) ap1-a1(i) -ap1*a5(w)-a3(u) ; 0 0 -1/tz];
                        F=[F A'*P+P*A<=-tol*eye(n)];
                    end
                end
            end
        end
    end
end

options=sdpsettings('solver','sedumi','verbose',0); % solver settings
sol=optimize(F,[],options); % optimize for quadratic stability
