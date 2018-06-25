% Can't do it by else d=1
% Use function found in reference of paper
% try only one point one range: settled at one position and can't be moved
% sloved by set different sub d
% principle: one IF THEN -> two constraints in terms of one flag variable 
%            different IF THEN use different flag variables

clear
debug = 1;
% Params
cyl_Num = 1;
cyl_r = 1;
Nl = 1;
Nr = 1;
cubeNumx = 1;
cubex = 2*cubeNumx+1;
cubeNumy = 1;
cubey = 2*cubeNumy+1;
cubeNumz = 1;
cubez = 2*cubeNumz+1;
cubeN = cubex*cubey*cubez; 
cube_side = 1;
ss_Num = (Nr+1)^Nl;
v_one_Num = 3;
v_Num = v_one_Num*Nl;  
constr_Num = 4; %linearize by 4 planars
theta = [20 0 -20]; theta = theta*pi/180;
dir = -ones(Nr+1,Nr);

for uptodown_r = 1:Nr+1
    if uptodown_r<=Nr
        dir(uptodown_r,uptodown_r:end) = -dir(uptodown_r,uptodown_r:end);
    end
end
DIR = zeros(ss_Num,Nl*Nr);
stack_N = ss_Num/(Nr+1); Rep_N = 1;
for uptodown_l = 1:Nl
    block = [];
    
    for j = 1:Nr+1
        block = [block; repmat(dir(j,:),Rep_N,1)];
    end
    
    DIR_ = repmat(block,stack_N,1);
    DIR(:,(uptodown_l-1)*Nr+1:uptodown_l*Nr) = DIR_;
    
    stack_N = stack_N/(Nr+1); Rep_N = Rep_N*(Nr+1);
end

CUBE = [];
%save cube pos in a matrix
for ca = 1:cubex
    for cb = 1:cubey
        for cc = 1:cubez
            x = ca-cubeNumx-1;
            y = cb-cubeNumy-1;
            z = cc-cubeNumz-1;
            cyl_f = 0;
                for j = 1:cyl_Num
                    if sqrt(x^2+y^2)<cyl_r*j+cube_side/2 && sqrt(x^2+y^2)>cyl_r*j-cube_side/2
                        cyl_f = 1;
                        break;
                    end
                end
            if cyl_f
            CUBE = [CUBE;[ca-cubeNumx-1 cb-cubeNumy-1 cc-cubeNumz-1 j]];
            end
        end
    end
end
% CUBE = CUBE;


M = 200;
pe = 0.01;

    clear model;

    
A5 = zeros(0,v_Num);
R5 = [];

S_A5 = sparse(A5);
IND_Z = [];
for cyl = 1:cyl_Num
    cyl_ind = find(CUBE(:,4)==cyl);
    sCUBE = CUBE(cyl_ind,:);
    
    ind_Z = [];
    A3 = [];
    R3 = [];
    S_A3 = sparse(A3);
    S_R3 = sparse(R3);
    pre_off_Z = 0;
    for s=1:ss_Num
    % s=1
        ud = DIR(s,:);

        A2 = [];
        R2 = [];
        S_A2 = sparse(A2);
        S_R2 = sparse(R2);
        
        cubeN = size(sCUBE,1);
        for c = 1:cubeN
    % c=1
            cx = sCUBE(c,1);
            cy = sCUBE(c,2);
            cz = sCUBE(c,3);

            % one cube in one subspace
            A4 = [];
            R4 = [];
            ind_l = [];
            for l = 1:Nl
                ud_one = ud((l-1)*Nr+1:l*Nr);
                A1 = [];
                R1 = [];
                for k = 1:constr_Num
                    A1_ = zeros(2*Nr,v_Num+constr_Num*(Nr+1)+1); %2:one IFTHEN->2constr % 1:one judge fk %1: summary judge F
                    F1_ = zeros(2,size(A1_,2));
                    R1_ = zeros(2*Nr,1);

                    for r = 1:Nr
                        B = ud_one(r);  % B==1 -> UP

        %                 B = 1;
                        angle = tan(theta(r));
                        if k==1 % cz > cx && cz > -cx
                           A1_((r-1)*2+1:r*2,(l-1)*v_one_Num+1:l*v_one_Num) = [-angle*B   0   B;
                                                                                angle*B   0  -B];
                           R1_((r-1)*2+1:r*2,1) = [(M-cx*angle*B+cz*B);
                                                   (cx*angle*B-cz*B-pe)];
                        elseif k==2 % cz > -cx && cz > cx
                           A1_((r-1)*2+1:r*2,(l-1)*v_one_Num+1:l*v_one_Num) = [ angle*B   0   B;
                                                                                -angle*B   0  -B];
                           R1_((r-1)*2+1:r*2,1) = [(M+cx*angle*B+cz*B);
                                                   (-cx*angle*B-cz*B-pe)];
                        elseif k==3 % cz > cy && cz > -cy
                           A1_((r-1)*2+1:r*2,(l-1)*v_one_Num+1:l*v_one_Num) = [0   -angle*B   B;
                                                                               0    angle*B  -B];
                           R1_((r-1)*2+1:r*2,1) = [(M-cy*angle*B+cz*B);
                                                   (cy*angle*B-cz*B-pe)];
                        elseif k==4 % cz > -cy && cz > cy
                           A1_((r-1)*2+1:r*2,(l-1)*v_one_Num+1:l*v_one_Num) = [0    angle*B   B;
                                                                               0   -angle*B  -B];
                           R1_((r-1)*2+1:r*2,1) = [(M+cy*angle*B+cz*B);
                                                   (-cy*angle*B-cz*B-pe)];
                        end
                    end
                    A1_(:,v_Num+(k-1)*(Nr+1)+1:v_Num+k*(Nr+1)-1) = kron(eye(Nr), [M;-(M+pe)]);
                    % Add two lines to do: IF a cube survives all the lasers for
                    %  kth side THEN fk is 1
                    F1_(1,v_Num+(k-1)*(Nr+1)+1:v_Num+k*(Nr+1)-1) = ones(1,Nr); F1_(1,v_Num+k*(Nr+1)) = -1;  % d1+d2+..+dNr-f<=Nr-1+pe
                    F1_(2,v_Num+(k-1)*(Nr+1)+1:v_Num+k*(Nr+1)-1) = -ones(1,Nr); F1_(2,v_Num+k*(Nr+1)) = M;  % -d1-d2-..-dNr+Mf<=M-Nr+pe
                    % x(v_Num+k*(Nr+1)) is fk that judges whether this cube satisfy
                    % all the lasers for kth side
                    A1 = [A1;A1_;F1_];
                    R1 = [R1;R1_;Nr-1+pe;M-Nr+pe];
                end



                % Add two lines to do: F = f1 or f2 or... or fconstr_Num
                F1_summary = zeros(2,v_Num+constr_Num*(Nr+1)+1);
                F1_summary(1,v_Num+Nr+1:(Nr+1):end) = -ones(1,constr_Num); F1_summary(1,end) = 1; %-f1-f2...-fconstrNum+F<=pe
                F1_summary(2,v_Num+Nr+1:(Nr+1):end) =  ones(1,constr_Num); F1_summary(2,end) = -M; %f1+f2+...+fconstrNum-MF<=pe
                A1 = [A1;F1_summary];
                R1 = [R1;pe;pe];

                % end of one lidar

                [H_onelidar,W_onelidar] = size(A1);
                W_onelidar = W_onelidar-v_Num;
                A4_ = zeros(H_onelidar,v_Num+ Nl*W_onelidar+1);%1: for flag of one cube one ss
                A4_(:,1:v_Num) = A1(:,1:v_Num);
                A4_(:,v_Num+(l-1)*W_onelidar+1:v_Num+l*W_onelidar) = A1(:,v_Num+1:end);
                A4 = [A4;A4_];
                R4 = [R4;R1];
                ind_l = [ind_l (size(A1,2)-v_Num)*l+v_Num]; % flag for one lidar
            end
            F4_summary = zeros(2,size(A4,2));
            F4_summary(1,ind_l) = ones(1,length(ind_l))*debug; F4_summary(1,end) =  -1*debug;
            F4_summary(2,ind_l) =-ones(1,length(ind_l))*debug; F4_summary(2,end) =  M*debug;
            A4 = [A4;F4_summary];
            R4 = [R4; pe+Nl-1; pe+M-Nl];
            S_A4 = sparse(A4);
            S_R4 = sparse(R4);
            A4_length = size(A4,2)-v_Num;
            % end for onecube in one subspace

            [H_onecube,W_onecube] = size(A4);
            W_onecube = W_onecube - v_Num;
            A2_ = zeros(H_onecube,v_Num+ cubeN*W_onecube);

            A2_(:,1:v_Num) = A4(:,1:v_Num);
            A2_(:,v_Num+(c-1)*W_onecube+1:v_Num+c*W_onecube) = A4(:,v_Num+1:end);
            S_A2_ = sparse(A2_);
    %         A2 = [A2;A2_];
    %         R2 = [R2 ;R4];
            S_A2 = [S_A2;S_A2_];
            S_R2 = [S_R2 ;S_R4];


        end
        % Add a line for Z>sum(cube in ss flag)
        F2_ = zeros(1,size(A2_,2));
        F2_(1,v_Num+A4_length:A4_length:end) = ones(1,cubeN);
        S_F2_ = sparse(F2_);
    %     A2 = [A2; F2_];
    %     R2 = [R2;0];
        S_A2 = [S_A2; S_F2_];
        S_R2 = [S_R2;sparse(0)];
        ind_Z = [ind_Z size(S_A2,1)+pre_off_Z];
        %end of one subspace

        [H_onesubspace,W_onesubspace] = size(S_A2);
        W_onesubspace = W_onesubspace-v_Num;
    %     A3_ = zeros(H_onesubspace,v_Num+ ss_Num*W_onesubspace);
    %     A3_(:,1:v_Num) = A2(:,1:v_Num);
    %     A3_(:,v_Num+(s-1)*W_onesubspace+1:v_Num+s*W_onesubspace) = A2(:,v_Num+1:end);
        S_A3_ = sparse(H_onesubspace,v_Num+ ss_Num*W_onesubspace);
        S_A3_(:,1:v_Num) = S_A2(:,1:v_Num);
        S_A3_(:,v_Num+(s-1)*W_onesubspace+1:v_Num+s*W_onesubspace) = S_A2(:,v_Num+1:end);

        S_A3 = [S_A3;S_A3_];
        S_R3 = [S_R3;S_R2];
        pre_off_Z = size(S_A3,1);
    end
    
    ind_Z = ind_Z+size(S_A5,1);
    
    Ling1 = sparse(size(S_A5,1),size(S_A3,2)-v_Num);
    Conf = S_A3(:,1:v_Num);
    Ling2 = sparse(size(Conf,1),size(S_A5,2)-v_Num);
    NLing = S_A3(:,v_Num+1:end);

    
    S_A5 = [S_A5 Ling1;
            Conf Ling2 NLing];
    
    IND_Z = [IND_Z ind_Z];
    
    R5 = [R5;full(S_R3)];
end
    
    S_A = S_A5;
    S_Z_col = sparse(size(S_A,1),1);
    S_Z_col(IND_Z) = sparse(-1);

%     
    S_A = [S_A S_Z_col];
    
%     AA = full(S_A);
    for i =1:length(S_Z_col)
        if S_Z_col(i)
            fprintf("%d\n",i);
        end
    end
    
    RHS = R5;
    
%     obj = zeros(1,size(A2_,2));
%     obj(1,v_Num+A1_length:A1_length:end) = ones(1,cubeN);
    obj = zeros(1,size(S_A,2));
    obj(end) = 1;
    vtype = [];
    lb = [];
    ub = [];
    for i = 1:Nl
    vtype = [vtype;'C';'C';'C'];
    lb = [lb;-cubeNumx;-cubeNumy;-cubeNumz];
    ub = [ub;cubeNumx;cubeNumy;cubeNumy];
    end
    state_L = size(S_A,2)-v_Num;
    for i=1:state_L
        vtype = [vtype;'B'];
        lb = [lb;0];
        ub = [ub;1];
    end
    vtype(end) = 'C';
    lb(end) = 0;
    ub(end) = cubeN;
    
    model.A = S_A;
    model.obj = obj;
    model.rhs = RHS;
    model.sense = '<';
    model.vtype = vtype;
    model.modelsense = 'min';
    model.lb = lb;
    model.ub = ub;
    
    gurobi_write(model, '/Library/gurobi752/mac64/examples/build/lidar_config2.lp');
    clear params;
    params.outputflag = 0;
    
    tic
    result = gurobi(model,params);
    toc
%  format short;  
    disp(result)

%     for v=1:length(model.obj)
    for v=1:v_Num
        fprintf(' %f\n', result.x(v));
    end

    fprintf('Obj: %f\n', result.objval);
    