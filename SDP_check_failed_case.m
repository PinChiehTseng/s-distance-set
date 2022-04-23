warning('off','all')


n = 5;
cases = [2, 4];
order = 10^(0);
num_dis = 2;
start = 0;
error_order = 10^(0);
d = -1; %test : compute A(n,d)
if d > -1
    if d == 0
        cases = 1:n;
    else
        cases = d:n;
    end
end

values = [];


for ind = 1
    dis = [cases, 0];

    [num, box] = num_var(n, dis);


    disp(num);
   
    numbers = 0;
    for i = 0: n
        disp(i);
        b = per(i, 0, 0, n, dis);
        B = sear(b, num, box);
        for j = 0: n
            c = per(0, j, 0, n, dis);
            C = sear(c, num, box);
            for t = 0: n
                a = per(i, j, t, n, dis);
                A = sear(a, num, box);
                if A ~= -1
                    if B ~= -1
                        numbers = numbers + 1;
                        if C ~= -1
                            numbers = numbers + 1;
                        end
                    end
                end
            end
        end
    end
        
    clear i j t b B a A;

    cvx_clear
    cvx_begin sdp
        cvx_precision('default')
        %cvx_precision('high')
        %cvx_precision('best')
        %cvx_precision('medium')
        %cvx_precision('low')
        %cvx_solver sdpt3
        %cvx_solver SeDuMi
        cvx_solver MOSEK
        variable x(num);
        dual variable m{2, (floor(n/2)+1)};
        dual variable Lin{1, numbers};
        
        lst = [];
        object = cvx(zeros(1, n+1));
        for i = 0: n
            c = per(i, 0, 0, n, dis);
            C = sear(c, num, box);
            if C ~= -1
                object(1, i+1) = comb(n,i)*x(C);
                lst = [lst, C];
            end
        end


        maximize(sum(object)*error_order)
        subject to
        %%definition condition
        condition = [x(1)];
        x(1) == 1;
        
        matrix = cell(2, floor(n/2)+1);
        F_tot = cell(2, floor(n/2)+1);
        for k = start: floor(n/2)
            disp(k);
            [matrix_1, matrix_2, F_1, F_2] = generate(x, n, k, num, box, dis, order);
            F_tot{1, k+1} = F_1;
            F_tot{2, k+1} = F_2;
            m{1, k+1}: matrix_2 - matrix_1 >= 0;
            m{2, k+1}: matrix_1 >= 0;
            matrix{1, k+1} = matrix_1;
            matrix{2, k+1} = matrix_2 - matrix_1;
        end
        
        clear k;
        
        number = 0;
        lin = cell(numbers, num);
        for i = 1: numbers
            for j = 1: num
                lin{i, j} = 0;
            end
        end
        clear i j;
        
        for i = 0: n
            disp(i);
            b = per(i, 0, 0, n, dis);
            B = sear(b, num, box);
            for j = 0: n
                c = per(0, j, 0, n, dis);
                C = sear(c, num, box);
                for t = 0: n
                    a = per(i, j, t, n, dis);
                    A = sear(a, num, box);
                    if A ~= -1
                        if B ~= -1
                            number = number + 1;
                            x(A) >= 0; 
                            Lin{1, number}: -x(A) + x(B) >= 0;
                            condition = [condition, x(A), -x(A)+x(B)];
                            lin{number, A} = lin{number, A} - 1;
                            lin{number, B} = lin{number, B} + 1;
                            if C ~= -1
                                number = number + 1;
                                Lin{1, number}: 1+x(A)-x(B)-x(C) >= 0;
                                condition = [condition, 1+x(A)-x(B)-x(C)];
                                lin{number, A} = lin{number, A} + 1;
                                lin{number, B} = lin{number, B} - 1;
                                lin{number, C} = lin{number, C} - 1;
                                lin{number, 1} = lin{number, 1} + 1;
                            end
                        else
                            x(A) == 0;
                            condition = [condition, x(A) == 0];
                        end
                    end
                end
            end
        end
        
        clear i j t b B a A;
        
    cvx_end
    
    values = [values, cvx_optval];
end

disp(max(values));
if num_dis == 2
    if max(values) <= (1/2)*(n^(2)-n+2)
        disp("hold");
    else
        disp("fail");
    end
end
if num_dis == 3
    if max(values) <= n + comb(n, 3)
        disp("hold");
    else
        disp("fail");
    end
end
if num_dis == 4
    if max(values) <= 1 + comb(n, 2) + comb(n, 4)
        disp("hold");
    else
        disp("fail");
    end
end

disp("check primal");
for i = condition
    if i < 0
        disp("primal error");
        disp(i);
    end
end
clear i;
for i = 1: 2
    for j = 1: floor(n/2)+1
        vector = eigs(matrix{i, j});
        for t = vector
            if t < 0
                disp(t);
            end
        end
    end
end
clear i j t vector;
disp("check dual");
error = 0;
check_self = 1;
for t = 1: num
    a = 0;
    for i = 1: 2
        for j = start+1: floor(n/2)+1
            a = a + trace(m{i, j}*transpose(F_tot{i, j}{1, t}));
        end
    end
    clear i j
    for i = 1: numbers
        a = a + Lin{1, i}*lin{i, t};
    end
    if ismember(t, lst)
        check = 0;
        for R = 0: n
            c = per(R, 0, 0, n, dis);
            C = sear(c, num, box);
            if C == t
                check = 1;
                a = a - comb(n, R);
            end
        end
        if check == 0
            a = a;
        end
        if a > 0
            if t ~= 1
                error = error + a;
            else
                check_self = check_self + a;
            end
        end
    end
end
check_self = check_self + 1;
                    
clear i j R c C a t;
        

disp("error");
disp(error);
disp("check_self");
disp(check_self);





function A = krav(x, n, k)
    ans = 0;
    for j = 0: k
        ans = ans + (((-1)^j))*comb(x,j)*comb(n-x, k-j);
    end
    A = ans;
end

function O = sear(it, num, box)
    P = -1;
    if it == "N"
        P = -1;
    else
        for i = 1: num
            if string(it) == box(i)
                P = i;
            end
        end
    end
    O = P;
end

function [O, D] = num_var(n, dis)
    tot = ["[0,0,0]"];
    for i = 0: n
        for j = 0: n
            for t = 0: n
                B = per(i, j, t, n, dis);
                if B ~= "N"
                    if ~ismember(string(B), tot)
                        tot = [tot, string(B)];
                    end
                end
            end
        end
    end
    C = size(tot);
    O = C(1, 2);
    D = tot;
end

function Y = per(i, j, t, n, dis)
    O = "1";
    if ~ismember(i, dis)
        O = "N";
    end
    if ~ismember(j, dis)
        O = "N";
    end
    if ~ismember(i+j-2*t, dis)
        O = "N";
    end
    if i+j-t > n
        O = "N";
    end
    if i-t < 0 || j-t < 0 || i+j-2*t < 0 || i+j-2*t > n
        O = "N";
    end
    V = perms([i j i+j-2*t]);
    a = [];
    item = 0;
    for p = 1: 6
        vector = V(p, :);
        m = vector(1, 1);
        l = vector(1, 2);
        k = ((vector(1, 3) - m - l) / (-2)) ;
        check_1 = 0: n;
        if (ismember(k, check_1) && (ismember(m, check_1) && ismember(l, check_1)))
            item = item + 1;
            a = [a; [m, l, k]];
        end
    end
    if item == 1
        [A, B, C] = a;
    end
    first = [];
    for q = 1: item
        first = [first, a(q, 1)];
    end
    one = min(first);
    f = [];
    f_num = 0;
    second = [];
    for r = 1: item
        if a(r, 1) == one
            f = [f; a(r, :)];
            f_num = f_num + 1;
            second = [second, a(r, 2)];
        end
    end
    if f_num == 1
        A = f(1);
        B = f(2);
        C = f(3);
    end
    two = min(second);
    g = [];
    g_num = 0;
    third = [];
    for s = 1: f_num
        if f(s, 2) == two
            g_num = g_num + 1;
            g = [g; f(s, :)];
            third = [third, f(s, 3)];
        end
    end
    if g_num == 1
        A = g(1);
        B = g(2);
        C = g(3);
    end
    three = min(third);
    h = [];
    h_num = 0;
    for l = 1: g_num
        if f(l, 3) == three
            h = f(l, :);
            h_num = 1;
        end
    end
    if h_num == 1
        A = h(1);
        B = h(2);
        C = h(3);
    end
    if O == "1"
        O = '['+string(A)+',' +string(B) +','+string(C) +']';
    end
    Y = O;
end

function B = comb_num(i, j, t, k, n)
    ans = 0;
    for u = 0:n
        ans = ans + ((-1)^(u-t))*comb(u,t)*comb(n-2*k,u-k)*comb(n-k-u,i-u)*comb(n-k-u,j-u);
    end
    B = ans;
end

function A = comb(a, b)
    ans = 0;
    if a >= 0 && b >= 0 && b <= a
        ans = nchoosek(a, b);
    end
    A = ans;
end

function [matrix_1, matrix_2, F_1, F_2] = generate(x, n, k, num, box, dis, order)
    goal_1 = cvx(zeros(n-2*k+1, n-2*k+1));
    goal_2 = cvx(zeros(n-2*k+1, n-2*k+1));
    F_1 = cell(1, num);
    F_2 = cell(1, num);
    for f = 1: num+1
        for g = 1: 2
            F_1{g, f} = zeros(n-2*k+1, n-2*k+1);
            F_2{g, f} = zeros(n-2*k+1, n-2*k+1);
        end
    end
    for i = k: n-k
        for j = k: n-k
            const = (comb(n-2*k, i-k)*comb(n-2*k, j-k))^(1/2);
            a = 0;
            b = 0;
            if ismember(i, dis) && ismember(j, dis)
                for t = 0:n
                    a_1 = per(i, j, t, n, dis);
                    A = sear(a_1, num, box);
                    if A ~= -1
                        a = a + comb_num(i, j, t, k, n) * x(A);  
                        F_2{1, A}(i-k+1, j-k+1) = F_2{1, A}(i-k+1, j-k+1) + comb_num(i, j, t, k, n) * order / const;
                        F_1{1, A}(i-k+1, j-k+1) = F_1{1, A}(i-k+1, j-k+1) - comb_num(i, j, t, k, n) * order / const;
                    end
                end
            end
            for t = 0:n
                if ismember(i+j-2*t, dis)
                    if i+j-2*t >= 0 && i + j - 2*t <= n
                        b_1 = per(i+j-2*t, 0, 0, n, dis);
                        B = sear(b_1, num, box);
                        if B ~= -1
                            b = b + comb_num(i, j, t, k, n) * x(B);
                            F_1{1, B}(i-k+1, j-k+1) = F_1{1, B}(i-k+1, j-k+1) + comb_num(i, j, t, k, n) * order / const;
                        end
                    end
                end
            end
            goal_1(i-k+1, j-k+1) = a*order / const;
            goal_2(i-k+1, j-k+1) = b*order / const;
        end
    end
    matrix_1 = goal_1;
    matrix_2 = goal_2;
 end