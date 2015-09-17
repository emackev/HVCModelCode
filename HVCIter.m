function [w xdyn] = HVCIter(p)
% Runs one iteration of the simulation.  p is a structure of parameters.
%
% Emily Mackevicius 1/14/2015, heavily copied from Hannah Payne's code
% which builds off Ila Fiete's model, with help from Michale Fee and Tatsuo
% Okubo. 

% redefining params that are used often outside the for loop
nsteps = p.nsteps;
n = p.n;
m = p.m;
w = p.w; 
wmax = p.wmax;
Wmax = wmax*m; 
eta = p.eta;
bdyn=p.input;

% initializing variables
xdyn=zeros(n,nsteps);
oldx = zeros(n,1);
oldy = zeros(n,1);

for t = 1:nsteps
    % Adaptation
    y = oldy + 1/p.tau*(-oldy+oldx);
    Aadapt = p.alpha * y; % adaptation
    
    % Net feedforward input.  
    B = bdyn(:,t); % external input
    AE = w*oldx; % excitatory input 
    AIff = p.beta * sum(oldx); % feed forward inhibition
    Anet = AE - AIff - Aadapt + B; % net feed forward input
    Anet(Anet < 0) = 0; % rectify
    
    % recurrent inhibition
    AIrec = p.gamma * sum(Anet); 
    
    % binary output
    x = (Anet - AIrec) > 0; 
    
    % STDP rule (Fiete et al 2010)
    dw_STDP = eta.*(x*(oldx)'-(oldx)*x');
    
    % Hetersynaptic penalty (Fiete et al 2010)
    dw_hLTDpre = eta*ones(n,1)*max(0, sum(w+dw_STDP,1)-Wmax);  % Weights leaving cells (pre)
    dw_hLTDpost = eta*max(0, sum(w+dw_STDP,2)-Wmax)*ones(1,n);  % Weights onto cells (post)
    
    % Update weights
    if eta>0
        dwtotal = dw_STDP-p.epsilon*(dw_hLTDpre+dw_hLTDpost);
        w = w + dwtotal;
        w(w > wmax) = wmax; % hard limit on strength of a single synapse
        w(w < 0) = 0; % weights cannot be negative
        w = w.*(~eye(p.n)); % clamp diagonal
    end
    
    oldx = double(x);
    oldy = y;
    xdyn(:,t)=x;
end




