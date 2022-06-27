function [Maxroot]=stability(Awave)

    [m,mp,k]=size(Awave);
    
%     for tt=1:k
%         A=squeeze(Awave(:,:,tt));
% %         eig_A(tt)= max(abs(eig([A; eye(mp-m,mp) ]))); 
%         eig_A(tt)= max(abs(eig([A; eye(mp-m,mp) ])));    
% 
%     end
%     Maxroot=eig_A;
    
    I = eye(mp-m,mp-m);
    O = zeros(mp-m,m);
    for tt=1:k
        A=squeeze(Awave(:,:,tt));
        eig_A(tt)=max(abs(eig([A; [I O]])));
    end
    Maxroot=eig_A;
    
    
%     figure, plot(eig_A)
end
