%--------------------------------------------------------------------------
%
% Accel_Grav_Harmonic_gottlieb: Computes the acceleration due to the spherical gravity
%                   harmonics field of the primary body
%                
%
% Inputs:
%   mu                  grav parameter of the central body
%   radius              radius of the central body
%   r_MI                Satellite position vector in the INERTIAL Frame
%                           centered on the body of interest 
%   C(n,m), S(n,m)      Spherical gravity harmonics coefficients (normalized)
%   degree              Maximum degree
%   order               Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a_GHAR           Acceleration (a=d^2r/dt^2) in INERTIAL frame
%
% Ref: "Normalization and Implementation of Three Gravitational Acceleration Models" -
% NASA - Reviewed by Robert G.Gottlieb
%
% Last modified:   8/Apr/2024
% 
%--------------------------------------------------------------------------


function a_GHAR = Accel_Grav_Harmonic_gottlieb(mu, radius, r_MI, C, S, degree, order, rot_mat)

%Extract the relevant Normalization Parameters computed before hand to save
%computation time
% norm1 = norm_param_gottlieb{1};
% norm2 = norm_param_gottlieb{2};
% norm11 = norm_param_gottlieb{3};
% normn10 = norm_param_gottlieb{4};
% norm1m = norm_param_gottlieb{5};
% norm2m = norm_param_gottlieb{6};
% normn1 = norm_param_gottlieb{7};

norm1 = zeros(degree+1);
norm2 = zeros(degree+1);
norm11 = zeros(degree+1);
normn10 = zeros(degree+1);
norm1m = zeros(degree+1,degree+1);
norm2m = zeros(degree+1,degree+1);
normn1 = zeros(degree+1,degree+1);
for n = 2:degree+1 
    norm1(n) = sqrt((2*n+1)/(2*n-1)); 
    norm2(n) = sqrt((2*n+1)/(2*n-3)); 
    norm11(n) = sqrt((2*n+1)/(2*n))/(2*n-1); 
    normn10(n) = sqrt((n+1)*n/2); 
    for m = 1:n 
        norm1m(n,m) = sqrt((n-m)*(2*n+1)/((n+m)*(2*n-1))); 
        norm2m(n,m) = sqrt((n-m)*(n-m-1)*(2*n+1)/((n+m)*(n+m-1)*(2*n-3))); 
        normn1(n,m) = sqrt((n+m+1)*(n-m)); 
    end 
end

x=rot_mat*r_MI; %Compute the position vector into Body-Fixed frame
r = sqrt(x(1)^2+x(2)^2+x(3)^2); %Norm of the position
ri=1/r;
xor=x(1)*ri;
yor=x(2)*ri;
zor=x(3)*ri;
ep=zor;
reor=radius*ri;
reorn=reor;
muor2=mu*ri*ri;

%Computation of the Associated Legendre Functions
p = zeros(degree+3,degree+3);
p(1,1) = 1; %RAE
p(1,2) = 0; %RAE
p(1,3) = 0; %RAE
p(2,2) = sqrt(3); %RAE %norm
p(2,3) = 0; %RAE
p(2,4) = 0; %RAE

for n = 2:degree %RAE
    ni = n+1; %RAE
    p(ni,ni) = norm11(n)*p(n,n)*(2*n-1); %RAE %norm
    p(ni,ni+1) = 0; %RAE
    p(ni,ni+2) = 0; %RAE
end

ctil = zeros(degree+1);
stil = zeros(degree+1);
ctil(1)=1; %RAE
stil(1)=0; %RAE
ctil(2)=xor; %RAE
stil(2)=yor; %RAE
sumh=0;
sumgm=1;
sumj=0;
sumk=0;
p(2,1) = sqrt(3)*ep; %RAE %norm

for n=2:degree
    ni=n+1; %RAE
    reorn=reorn*reor;
    n2m1=n+n-1;
    nm1=n-1;
    np1=n+1;
    p(ni,n) = normn1(n,n-1)*ep*p(ni,ni); %RAE %norm
    p(ni,1) = (n2m1*ep*norm1(n)*p(n,1)-nm1*norm2(n)*p(nm1,1))/n; %RAE %norm
    p(ni,2) = (n2m1*ep*norm1m(n,1)*p(n,2)-n*norm2m(n,1)*p(nm1,2))/(nm1); %RAE %norm
    sumhn=normn10(n)*p(ni,2)*C(ni,1); %norm %RAE
    sumgmn=p(ni,1)*C(ni,1)*np1; %RAE
    if (order>0)
        for m = 2:n-2
            mi = m+1; %RAE
            p(ni,mi) = (n2m1*ep*norm1m(n,m)*p(n,mi)-...
                (nm1+m)*norm2m(n,m)*p(nm1,mi))/(n-m); %RAE %norm
        end
        sumjn=0;
        sumkn=0;
        ctil(ni)=ctil(2)*ctil(ni-1)-stil(2)*stil(ni-1); %RAE
        stil(ni)=stil(2)*ctil(ni-1)+ctil(2)*stil(ni-1); %RAE
        if(n<order)
            lim=n;
        else
            lim=order;
        end
        for m=1:lim
            mi=m+1; %RAE
            mm1=mi-1; %RAE
            mp1=mi+1; %RAE
            mxpnm=m*p(ni,mi); %RAE
            bnmtil=C(ni,mi)*ctil(mi)+S(ni,mi)*stil(mi); %RAE
            sumhn=sumhn+normn1(n,m)*p(ni,mp1)*bnmtil; %RAE %norm
            sumgmn=sumgmn+(n+m+1)*p(ni,mi)*bnmtil; %RAE
            bnmtm1=C(ni,mi)*ctil(mm1)+S(ni,mi)*stil(mm1); %RAE
            anmtm1=C(ni,mi)*stil(mm1)-S(ni,mi)*ctil(mm1); %RAE
            sumjn=sumjn+mxpnm*bnmtm1;
            sumkn=sumkn-mxpnm*anmtm1;
        end
        sumj=sumj+reorn*sumjn;
        sumk=sumk+reorn*sumkn;
    end
    sumh = sumh+reorn*sumhn;
    sumgm = sumgm+reorn*sumgmn;
end

lambda=sumgm+ep*sumh;
g = zeros(3,1);
g(1,1)=-muor2*(lambda*xor-sumj);
g(2,1)=-muor2*(lambda*yor-sumk);
g(3,1)=-muor2*(lambda*zor-sumh);

a_GHAR=rot_mat'*g; %RAE

end




