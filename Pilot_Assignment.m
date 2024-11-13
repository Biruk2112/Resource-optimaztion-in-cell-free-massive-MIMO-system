function Phi=Pilot_Assignment(K,tau_p)
Phi=zeros(tau_p,K);
% C=PN_Sequences(K,tau_p); Non-orthogonal
C=sign(randn(K,tau_p));
for kk=1:K
    Phi(:,kk)=C(kk,:);
    Phi(:,kk)=1/sqrt(tau_p)*Phi(:,kk);
end
end