function Phi=Pilot_Assignment_Orthogonal(K,tau_p)
Phi=zeros(tau_p,K); %Orthogonal 
W=hadamard(tau_p);
for kk=1:K
    Phi(:,kk)=W(:,kk);
    Phi(:,kk)=1/sqrt(tau_p)*Phi(:,kk);
end
end