warning('on','all')

n = 4;
d = [2, 4];
w = 16;
d = [d, 0];
order = 10^(-3);
avoid = 0;
dis = -1; %test if dis > -1
if dis > -1
    d = [];
    for k = 0: n
        if (k == 0 || k >= dis) && mod(k, 2) == 0
            d = [d, k];
        end
    end
    clear k;
end
        
[num, box] = num_var(n, w, d);
disp(num)
value = [];
for n_1 = w
    n_2 = n-n_1;


    clear a b c;

    disp(num);

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

        object = cvx(zeros(1, min(n_1, n-n_1)+1));
        for i = 0: min(n_1, n_2)
            c = per(i, 0, 0, 0, n, n_1, d);
            C = sear(c, num, box);
            if C ~= -1
                object(1, i+1) = comb(n_1,i)*comb(n_2, i)*x(C);
            end
        end

        clear i c C;


        maximize(sum(object))
        subject to
        %%definition condition
        x(1) == 1;
        for k_1 = 0: floor(n_1/2)
            for k_2 = 0: floor(n_2/2)
                if avoid == 1
                    if k_1 ~= 0 || k_2 ~= 0
                        disp([k_1, k_2]);
                        [matrix_1, matrix_2] = generate(x, n, k_1, k_2, num, box, n_1, d);
                        (matrix_2 - matrix_1) * order >= 0;
                        matrix_1 * order >= 0;
                    end
                else
                    disp([k_1, k_2]);
                    [matrix_1, matrix_2] = generate(x, n, k_1, k_2, num, box, n_1, d);
                    (matrix_2 - matrix_1) * order >= 0;
                    matrix_1 * order >= 0;
                end
            end
        end

        clear k_1 k_2;

        for i = 0: min(n_1, n-n_1)
            disp(i);
            b = per(i, 0, 0, 0, n, n_1, d);
            B = sear(b, num, box);
            for j = 0: min(n_1, n-n_1)
                c = per(j, 0, 0, 0, n, n_1, d);
                C = sear(c, num, box);
                for t = 0: min(n_1, n-n_1)
                    for s = 0: min(n_1, n-n_1)
                        a = per(i, j, t, s, n, n_1, d);
                        A = sear(a, num, box);
                        if A ~= -1
                            if B ~= -1
                                x(A) >= 0; 
                                -x(A) + x(B) >= 0;
                            end
                        end
                    end
                end
            end
        end
        clear i j t s a b A B;

    cvx_end
    value = [value, cvx_optval];
end

disp(min(value));





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

function [O, D] = num_var(n, w, d)
    tot = ['[0,0,0,0]'];
    for i = 0: min(w, n-w)
        for j = 0: min(w, n-w)
            for t = 0: min(w, n-w)
                for s = 0: min(w, n-w)
                    B = per(i, j, t, s, n, w, d);
                    if B ~= "N"
                        if ~ismember(string(B), tot)
                            tot = [tot, string(B)];
                        end
                    end
                end
            end
        end
    end
    C = size(tot);
    O = C(1, 2);
    D = tot;
end

function Y = per(i, j, t, s, n, w, d)
    O = "1";
    if ~ismember(2*i, d)
        O = "N";
    end
    if ~ismember(2*j, d)
        O = "N";
    end
    if ~ismember(2*(i+j-(t+s)), d)
        O = "N";
    end
    if i+j-(t+s) < 0 || i+j-(t+s) > min(w, n-w)
        O = "N";
    end
    if i-t < 0 || j-t < 0 || i-s < 0 || j-s < 0 || i+j-t > w || i+j-s > n-w || i+j-w+s > n
        O = "N";
    end
    if i+j-2*t > w || i+j-2*s > n-w
        O = "N";
    end
    V = perms([i j i+j-(t+s)]);
    a = [];
    K = t-s;
    item = 0;
    for p = 1: 6
        vector = V(p, :);
        m = vector(1, 1);
        l = vector(1, 2);
        k = (-(vector(1, 3) - m - l) + K) / 2;
        z = (-(vector(1, 3) - m - l) - K) / 2;
        check_1 = 0: min(w, n-w);
        if (ismember(k, check_1) && (ismember(m, check_1) && ismember(l, check_1) && ismember(z, check_1)))
            item = item + 1;
            a = [a; [m, l, k, z]];
        end
    end
    if item == 0
        Q = "N";
    end
    if item == 1
        [A, B, C, D] = a;
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
        D = f(4);
    end
    two = min(second);
    g = [];
    g_num = 0;
    third = [];
    for L = 1: f_num
        if f(L, 2) == two
            g_num = g_num + 1;
            g = [g; f(L, :)];
            third = [third, f(L, 3)];
        end
    end
    if g_num == 1
        A = g(1);
        B = g(2);
        C = g(3);
        D = g(4);
    end
    three = min(third);
    h = [];
    h_num = 0;
    fourth = [];
    for l = 1: g_num
        if f(l, 3) == three
            h = [h; g(l, :)];
            h_num = h_num+1;
            fourth = [fourth, g(l, 4)];
        end
    end
    if h_num == 1
        A = h(1);
        B = h(2);
        C = h(3);
        D = h(4);
    end
    four = min(fourth);
    I = [];
    P = 0;
    for X = 1: h_num
        if h(X, 4) == four
            P = 1;
            I = h(X, :);
        end
    end
    if P == 1
        A = I(1);
        B = I(2);
        C = I(3);
        D = I(4);
    end
    if item == 0
        O = "N";
    end
    if O == "1"
        O = '['+string(A)+',' +string(B) +','+string(C) + ','+string(D) +']';
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

function [matrix_1, matrix_2] = generate(x, n, k_1, k_2, num, box, w, d)
    v = n-w;
    K = min(w-k_1, v-k_2);
    N = max(k_1, k_2);
    goal_1 = cvx(zeros((K-N+1), (K-N+1)));
    goal_2 = cvx(zeros((K-N+1), (K-N+1)));
    n = w+v;
    for i = N: K
        for j = N: K
            a = 0;
            b = 0;
            const = comb(w-2*k_1, i-k_1)*comb(w-2*k_1, j-k_1)*comb(v-2*k_2, i-k_2)*comb(v-2*k_2, j-k_2);
            for t_1 = 0: min(w, v)
                for t_2 = 0: min(w, v)
                    if ismember(2*i, d) && ismember(2*j, d)
                        a_1 = per(i, j, t_1, t_2, n, w, d);
                        A = sear(a_1, num, box);
                        if A ~= -1
                            a = a + comb_num(i, j, t_1, k_1, w)*comb_num(i, j, t_2, k_2, v)*x(A);  
                        end
                    end
                    if ismember(2*(i+j-t_1-t_2), d)
                        if i+j-(t_1+t_2) >= 0 && i + j - (t_1+t_2) <= min(w, v)
                            b_1 = per(i+j-(t_1+t_2), 0, 0, 0, n, w, d);
                            B = sear(b_1, num, box);
                            if B ~= -1
                                b = b + comb_num(i, j, t_1, k_1, w)*comb_num(i, j, t_2, k_2, v)*x(B);
                            end
                        end
                    end
                end
            end
            goal_1(i-N+1, j-N+1) = a * (const^(-1/2));
            goal_2(i-N+1, j-N+1) = b * (const^(-1/2));
        end
    end
    matrix_1 = goal_1;
    matrix_2 = goal_2;
 end