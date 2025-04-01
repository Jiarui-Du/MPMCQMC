function u = seqfunMH(sequence,N,d,n)
switch sequence
    case "iid"
        u = rand(N,d);
    case "lattice"
        latticeseq_b2 init0
        qu = latticeseq_b2(d, N);
        % random shift
        w = rand(1,d);
        W = repmat(w,N,1);
        u = mod(qu'+W,1);
    case "halton"
        halton = haltonset(d);
        qu = net(halton,N);
        % random shift
        w = rand(1,d);
        W = repmat(w,N,1);
        u = mod(qu+W,1);
    case "sobol-owen"
        sobol = sobolset(d);
        % qu = net(sobol,N);
        % random shift
        % w = rand(1,d);
        % W = repmat(w,N,1);
        % u = mod(qu+W,1);
        u = net(scramble(sobol,'MatousekAffineOwen'),N);
    case "sobol-Liao"
        sobol = sobolset(d);
        uu = net(sobol,N);
        p = randperm(N);
        qu = uu(p,:);
        % random shift
        w = rand(1,d);
        W = repmat(w,N,1);
        u = mod(qu+W,1);
    case "FELFSR"
        [~,bi] = FELFSR(n,N);
        % digital shift
        w = rand(1,d);
        k = n;
        bb = dec2base(floor(w*2^k),2,k);
        b = logical(double(bb)-48);
        mm = N-1;
        s = gcd(d,mm);
        u = zeros(mm,d);
        powertimes = 2.^(-(1:k));
        ds = d-1+s-1;
        if ds <= mm
            bi = [bi,bi(:,1:ds)];
        else
            dsm = mod(ds,mm);
            k = (ds-dsm)/mm;
            bi = [repmat(bi,1,k+1),bi(:,1:dsm)];
        end
        bb = (N-1)/s;
        for i = 1:s
            for j = 1:bb
                temp1 = mod((j-1)*d+1,mm);
                temp1 = temp1+(temp1==0)*mm;
                temp1 = temp1+i-1;
                temp2 = temp1+d-1;
                tempb = bi(:,temp1:temp2);
                qub = xor(tempb,b');
                u((i-1)*bb+j,:) = powertimes*qub;
            end
        end
        u = [w;u];
        if ~all(all(u)) 
            u(u==0) = rand*2^(-(n+1)); 
        end
    case "HaraseF2"
        [~,bi] = HaraseF2(n,N);
        % digital shift
        w = rand(1,d);
        k = 32;
        bb = dec2base(floor(w*2^k),2,k);
        b = logical(double(bb)-48);
        mm = N-1;
        s = gcd(d,mm);
        u = zeros(mm,d);
        powertimes = 2.^(-(1:k));
        ds = d-1+s-1;
        if ds <= mm
            bi = [bi,bi(:,1:ds)];
        else
            dsm = mod(ds,mm);
            k = (ds-dsm)/mm;
            bi = [repmat(bi,1,k+1),bi(:,1:dsm)];
        end
        bb = (N-1)/s;
        for i = 1:s
            for j = 1:bb
                temp1 = mod((j-1)*d+1,mm);
                temp1 = temp1+(temp1==0)*mm;
                temp1 = temp1+i-1;
                temp2 = temp1+d-1;
                tempb = bi(:,temp1:temp2);
                qub = xor(tempb,b');
                u((i-1)*bb+j,:) = powertimes*qub;
            end
        end
        u = [w;u];
        if ~all(all(u)) 
            u(u==0) = 1e-10; 
        end
end