function[x]=sorL(A,b,xold,erromax,fr)
n = length(A);
% Calculando B e g
for i = 1:n
    B(i,:) = A(i,:)/A(i,i); 
    g(i,1) = b(i)/A(i,i);
end
B = eye(n)-B;
% Verificando se é convergente
% Critério de Sassenfeld
beta = zeros(n,1);
for i = 1:n
    soma = 0.;
    for j = 1:i-1
    soma = soma + beta(j) * abs(B(i,j));
    end
    for j = i+1:n
    soma = soma + abs(B(i,j));
    end
beta(i) = soma;
end
beta
max(beta)
if (max(beta) < 1)
    disp('Converge')
else
    disp('Não converge')
    return
end
disp('Continuando')
x = xold;
for itera = 1:1000
    for i = 1:n
    soma = 0.;
    for j = 1:n
    soma = soma + B(i,j)*x(j);
    end
    %x(i) = (1.-fr) * x(i) + fr*(soma + g(i));
    x(i) = x(i) + fr*(soma + g(i)-x(i));
end
x;
erro = max(abs(x-xold))/max(max(abs(x)),1)
if (erro<erromax)
    itera
    break;
end
xold = x;
end
end