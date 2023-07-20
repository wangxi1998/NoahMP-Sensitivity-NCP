function [FOSI,TSI] = Sobol(X,Y)
%% *** Preliminary Statistics *********************************************
[Nsamples,Nparms] = size(Y);
Nparms = (Nparms-2)/2;
Y_single = [Y(:,1);Y(:,2)];
fo = mean(Y_single);
Do = var(Y_single);
% % wangxi modified this function
D = zeros(Nparms,2);
for p = 1:Nparms
    for n = 1:Nsamples
        D(p,1) = D(p,1) + Y(n,2)*(Y(n,Nparms+p+2)-Y(n,1));
        D(p,2) = D(p,2) + (Y(n,1)-Y(n,Nparms+p+2))^2;
    end
end
D = D./Nsamples;
FOSI(:)=D(:,1)./Do ;
TSI(:)=0.5*D(:,2)./Do ;