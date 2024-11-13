function Phi=Pilot_Assignment_Orthogonal_random(K,tau_p)
Phi=zeros(tau_p,K); %Greedy PA
W=hadamard(tau_p);
for kk=1:tau_p
    Phi(:,kk)=W(:,kk);
    Phi(:,kk)=1/sqrt(tau_p)*Phi(:,kk);
end

for kk=tau_p+1:K
    index_ii=kk-tau_p;
    Phi(:,kk)=W(:,index_ii);
    Phi(:,kk)=1/sqrt(tau_p)*Phi(:,kk);
end


end