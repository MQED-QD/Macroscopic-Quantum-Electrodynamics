function Ideal_V_Coupling(Angle_Set)

c = 2.99792458e8;

% hbar = 1.0545718e-34;
% eV = 1.602176634*1e-19;

epsilon0 = 8.8541878128 * 1e-12;
eV = 1;

N = 1e3;

eps = zeros(N,1);          % site energies (可改成你的)
J   = zeros(N-1,1);        % bond couplings H(m,m+1), m=1..N-1

% 例：每条 bond 不同（你可替换成来自 beta 的表达式）
% g1 = (5e-4*eV)^2;
% g1 = 0;

g0 = (5e-4*eV)^2;
g0 = 0;

g1 = (1e-3*eV)^2;

eps = eps + sqrt(g0)*randn(N,1);

beta = sqrt(g1)*randn(N-1,1);   % 只存上对角的无序

J = (-1e-3*eV) + beta;             % 举例：J_m = J0 + beta_m

% J = (0*eV) + beta;     


% 如果你还需要非厄米对角项，例如 eps = eps - 1i*gamma/2; 也可直接写进 eps

saveFileName = sprintf('Parameter_Set%.0f.mat', Angle_Set);
save(saveFileName, 'eps','J','-v7.3');  % v7.3 用 HDF5，Python 读取很方便

 
% Check_g1_Single = diag(Vab,1);

% Var_2 = var(Check_g1_Single)/eV^2
% mean(Check_g1_Single)/eV
% histogram(Check_g1,50)


end