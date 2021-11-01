function A = SISOsubid(y,u,n)
  y = y(:); y = y'; ny = length(y); % turn y into row vector
  u = u(:); u = u'; nu = length(u); % turn u into row vector
  i = 2*n; % #rows in Hankel matrices. Typically: i = 2 * (max order)
  twoi = 4*n;           

  if ny ~= nu, error('y and u must be same size'); end
  if ((ny-twoi+1) < twoi); error('Not enough data points'); end

  % Determine the number of columns in the Hankel matrices
  j = ny-twoi+1;

  % Make Hankel matrices Y and U
  Y=zeros(twoi,j); U=zeros(twoi,j);
  for k=1:2*i
    Y(k,:)=y(k:k+j-1); U(k,:)=u(k:k+j-1);
  end
  % Compute the R factor
  R = triu(qr([U;Y]'))'; % R factor
  R = R(1:4*i,1:4*i); 	 % Truncate

  % ------------------------------------------------------------------
  % STEP 1: Calculate oblique and orthogonal projections
  % ------------------------------------------------------------------
  Rf = R(3*i+1:4*i,:);              % Future outputs
  Rp = [R(1:1*i,:);R(2*i+1:3*i,:)]; % Past inputs and outputs
  Ru  = R(1*i+1:2*i,1:twoi); 	      % Future inputs
  % Perpendicular future outputs 
  Rfp = [Rf(:,1:twoi) - (Rf(:,1:twoi)/Ru)*Ru,Rf(:,twoi+1:4*i)]; 
  % Perpendicular past inputs and outputs
  Rpp = [Rp(:,1:twoi) - (Rp(:,1:twoi)/Ru)*Ru,Rp(:,twoi+1:4*i)]; 

  % The oblique projection is computed as (6.1) in VODM, page 166.
  % obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
  % The extra projection on Ufp (Uf perpendicular) tends to give 
  % better numerical conditioning (see algo on VODM page 131)

  % Funny rank check (SVD takes too long)
  % This check is needed to avoid rank deficiency warnings
  if (norm(Rpp(:,3*i-2:3*i),'fro')) < 1e-10
    Ob = (Rfp*pinv(Rpp')')*Rp; 	% Oblique projection
  else
    Ob = (Rfp/Rpp)*Rp;
  end

  % ------------------------------------------------------------------
  % STEP 2: Compute weighted oblique projection and its SVD
  %         Extra projection of Ob on Uf perpendicular
  % ------------------------------------------------------------------
  WOW = [Ob(:,1:twoi) - (Ob(:,1:twoi)/Ru)*Ru,Ob(:,twoi+1:4*i)];
  [U,S,~] = svd(WOW);
  ss = diag(S);

  % ------------------------------------------------------------------
  % STEP 3: Partitioning U into U1 and U2 (the latter is not used)
  % ------------------------------------------------------------------
  U1 = U(:,1:n); % Determine U1

  % ------------------------------------------------------------------
  % STEP 4: Determine gam = Gamma(i) and gamm = Gamma(i-1) 
  % ------------------------------------------------------------------
  gam  = U1*diag(sqrt(ss(1:n)));
  gamm = gam(1:(i-1),:);
  gam_inv  = pinv(gam); 			% Pseudo inverse of gam
  gamm_inv = pinv(gamm); 			% Pseudo inverse of gamm

  % ------------------------------------------------------------------
  % STEP 5: Determine A matrix (also C, which is not used) 
  % ------------------------------------------------------------------
  Rhs = [gam_inv*R(3*i+1:4*i,1:3*i),zeros(n,1); R(i+1:twoi,1:3*i+1)];
  Lhs = [gamm_inv*R(3*i+1+1:4*i,1:3*i+1); R(3*i+1:3*i+1,1:3*i+1)];
  sol = Lhs/Rhs;    % Solve least squares for [A;C]
  A = sol(1:n,1:n); % Extract A
return