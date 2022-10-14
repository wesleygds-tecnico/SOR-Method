function [] = SOR(A,B,x,it,omega)
n =  length(A);
E = A;
betha = zeros(n,1);
betha(1,1) = 1;
for i = 1:n
    soma = 0;
    for j = 1:(i-1)
        soma = soma + abs(A(i,j))*betha(j);
    end
    for j = i+1:n
        soma = soma + abs(A(i,j));
    end
    betha(i,1) = soma/A(i,i);
    if betha(1:n,1) < 1 
        a = 1;
    else 
        fprintf('Não converge \n')
        a = 0;
        break
    end
end
c = 1;
for k = 1:n
    F = eig(A);
    C(1:k,1:k) = A(1:k,1:k);
    if (A(1,n) == 0 && A(n,1) == 0 && F(k,1)>0 && det(C)>0)
    else 
        c = 0;
    end
end
if a == 1
    fprintf('Converge \n')
    if (c == 1)
        L = zeros(n);
        U = zeros(n);
        for j = 1:(n-1)
            for i = 2:n
                L(i,1:(i-1)) = -E(i,1:(i-1));
                U(j,(j+1):n) = -E(j,(j+1):n);
            end
        end
        D =E+(L+U);
        Tj = (L+U)/D;
        c = eig(Tj);
        for k = 1:n
            if((abs(c(k,1)) > 10^-3) && ((c(k,1)+conj(c(k,1)))/2) == c(k,1))
                o = c(k,1);
                break
            end
        end
        omega = 2/(1+((1-(o^2))^(1/2))) %fator de relaxação
    end
    for k = 1:1000000
        for i = 1:n
            y=x;
            soma = 0;
            for j = 1:(i-1)
                soma = soma + A(i,j)*x(j,1);
            end
            for j = i+1:n
                soma = soma + A(i,j)*x(j,1);
            end
            x(i,1) = ((1-omega)*y(i,1))+((omega*(B(i,1) - soma))/A(i,i));
        end
        erro = max(abs(x-y))/max(max(abs(x)),1);
        if (erro<10^-it)
            break
        end
    end
    fprintf('Número de iterações \n')
    disp(k)
    fprintf('Solução \n')
    disp(x)
    fprintf('Verificação \n')
    a = A*x;
    disp(a)
end
end