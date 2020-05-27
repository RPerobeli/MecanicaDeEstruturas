clear all
close all
clc

%%Código para o programa de Análise Estrutural Utilizando Método da Rigidez
% Abrir o aqrquivo txt de entrada
entrada = fopen('D:\UFJF\Mecânica das Estruturas\Exercicio 6.txt','r'); 

% Parte 1
% Leitura e Armazenamento das Informações de Estrutura
% Leitura e Armazenamento do Titulo
titulo = fscanf(entrada,'%s',1); 

% Verificar qual tipo de estrutura será analisada
if (titulo == 'Trel')
   
    % Parte 1
    % Leitura e Armazenamento das Informações de Estrutura
    % Leitura e Armazenamento do número de nós
    num_no = fscanf(entrada,'%i',1); 
    
    % Leitura e Armazenamento do numero de elementos
    num_el = fscanf(entrada,'%i',1); 
    
    % Armazenamento das coordenadas dos nós (nome do no, coordenadas)
    coordenada = fscanf(entrada,'%f',[4,num_no]);  
    % Alteração o formato da Matriz para o adequado a ser usado em calculos
    coordenada = coordenada'; 

    % Armazenamento das restrições dos nós(nome do no, restrição)
    restricoes = fscanf(entrada,'%i',[4,num_no]);
    restricoes = restricoes';

    % Armazenamento dos carregamentos
    carregamento = fscanf(entrada,'%f',[4,num_no]);
    carregamento= carregamento';

    % Conectividades (numero do elemento,no_saida,no_chegada, tipo de material)
    conectividade = fscanf(entrada,'%i',[4,num_el]);
    conectividade = conectividade';

    %Armazenamento das características do material 
    %Tipo de material,Seção transversal,Mod Elasticidade, Mod Transversal, Momento de inercia
    num_material = fscanf(entrada,'%i',1); %Quantidade de tipos de materiais
    material = fscanf(entrada,'%f',[5,num_material]);
    material = material';

    %Fechar arquivo de entrada
    fclose(entrada);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parte 2
    % Renumeração da matriz de restrição
    num_eq = 0;
    for i = 1:num_no
        for j = 2:4
           if(restricoes(i,j) == 1)
               restricoes(i,j)=0;
           else
               restricoes(i,j)= num_eq + 1;
               num_eq = num_eq+1;
           end
        end
    end

    % Renumeração dos nós restritos
    dimensao_K_G = num_eq; %contador que define a dimensão na matriz Global K_G
    for i = 1:num_no
        for j = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
           if(restricoes(i,j)==0)
               restricoes(i,j) = dimensao_K_G+1;
               dimensao_K_G = dimensao_K_G+1;
           end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parte 3
    % Montagem do vetor de carregamentos externos
    F = zeros(dimensao_K_G,1);
    for i=1:dimensao_K_G
        for j = 1:num_no
            for k = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
                if(restricoes(j,k) == i)
                    F(i) = carregamento(j,k);
                end    
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Parte 4
    % Matriz de Rigidez Global e Vetor LM
    L = zeros(num_el,1);        % Comprimento do elemento
    R_L = zeros(6,6,num_el);    % Matriz de Rotação
    K_L = zeros(6,6,num_el);    % Matriz de Rigidez no Referencial Local
    K_R = zeros(6,6,num_el);    % Matriz de Rigidez no Referencial Local Rotacionada
    K_g = zeros(dimensao_K_G,dimensao_K_G);         % Matriz de Rigidez Global Auxiliar
    K_G = zeros(dimensao_K_G,dimensao_K_G);         % Matriz de Rigidez Global
    LM = zeros(num_el,6);       % Vetor LM

    for i = 1:num_el 
        %Abre um loop para cada elemento, a cada iteração será calculado o
        %comprimento do elemento, a matriz de rigidez local, a matriz de
        %rotação, matriz de rigidez rotacionada, o vetor LM, e parte da matriz
        %global. Ao final do loop, a matriz global K_G é obtida.

        % Cálculo do comprimento de cada elemento
        x_i = coordenada(conectividade(i,2),2);
        y_i = coordenada(conectividade(i,2),3);
        z_i = coordenada(conectividade(i,2),4);
        x_f = coordenada(conectividade(i,3),2);
        y_f = coordenada(conectividade(i,3),3);
        z_f = coordenada(conectividade(i,3),4);
        L(i) = sqrt((x_f-x_i)^2+(y_f-y_i)^2+(z_f-z_i)^2);

        % Definição da matriz de rigidez local da treliça espacial  
        M = (material(conectividade(i,4),2))*(material(conectividade(i,4),3))/L(i);
        K_L(1,1,i)= M; K_L(1,4,i)=-M; K_L(4,1,i)=-M; K_L(4,4,i)=M;

        % Definição da matriz de rotação local da treliça espacial    
        Cx = (x_f-x_i)/L(i);
        Cy = (y_f-y_i)/L(i);
        Cz = (z_f-z_i)/L(i);
        if(Cx<0)
            Cxz = -sqrt((Cx)^2 + (Cz)^2);
        else
            Cxz = sqrt((Cx)^2 + (Cz)^2);
        end
        if(Cx~=0)
            %Matriz de rotação para treliça espacial (GERE E WEAVER,1965) 
            
            R_L(1,1,i)=Cx;                  R_L(4,4,i)=R_L(1,1,i); 
            R_L(2,2,i)=Cxz;                 R_L(5,5,i)=R_L(2,2,i);
            R_L(3,3,i)=Cx/Cxz;              R_L(6,6,i)=R_L(3,3,i); 
            R_L(1,2,i)=Cy;                  R_L(4,5,i)=R_L(1,2,i);
            R_L(1,3,i)=Cz;                  R_L(4,6,i)=R_L(1,3,i);
            R_L(2,3,i)=-Cy*Cz/Cxz;          R_L(5,6,i)=R_L(2,3,i);
            R_L(2,1,i)=-Cx*Cy/Cxz;          R_L(5,4,i)=R_L(2,1,i);
            R_L(3,1,i)=-Cz/Cxz;             R_L(6,4,i)=R_L(3,1,i);
            R_L(3,2,i)=0;                   R_L(6,5,i)=R_L(3,2,i);
        else
            %Caso o elemento esteja na vertical, faz-se necessário o uso de uma
            %matriz adaptada
            R_L(1,1,i)=0;                   R_L(4,4,i)=R_L(1,1,i); 
            R_L(2,2,i)=0;                   R_L(5,5,i)=R_L(2,2,i);
            R_L(3,3,i)=1;                   R_L(6,6,i)=R_L(3,3,i); 
            R_L(1,2,i)=Cy;                  R_L(4,5,i)=R_L(1,2,i);
            R_L(1,3,i)=0;                   R_L(4,6,i)=R_L(1,3,i);
            R_L(2,3,i)=0;                   R_L(5,6,i)=R_L(2,3,i);
            R_L(2,1,i)=-Cy;                 R_L(5,4,i)=R_L(2,1,i);
            R_L(3,1,i)=0;                   R_L(6,4,i)=R_L(3,1,i);
            R_L(3,2,i)=0;                   R_L(6,5,i)=R_L(3,2,i);
        end


        % Matriz de rigidez local rotacionada para o referencial global
        K_R(1,1,i)=M*Cx^2; K_R(4,4,i)=M*Cx^2; 
        K_R(2,2,i)=M*Cy^2; K_R(5,5,i)=M*Cy^2;
        K_R(3,3,i)=M*Cz^2; K_R(6,6,i)=M*Cz^2; 
        K_R(1,2,i)=M*Cx*Cy; K_R(2,1,i)=M*Cx*Cy; K_R(4,5,i)=M*Cx*Cy; K_R(5,4,i)=M*Cx*Cy;
        K_R(1,3,i)=M*Cx*Cz; K_R(3,1,i)=M*Cx*Cz; K_R(6,4,i)=M*Cx*Cz; K_R(4,6,i)=M*Cx*Cz;
        K_R(2,3,i)=M*Cz*Cy; K_R(3,2,i)=M*Cz*Cy; K_R(5,6,i)=M*Cz*Cy; K_R(6,5,i)=M*Cz*Cy;
        K_R(4,1,i)=-M*Cx^2; K_R(1,4,i)=-M*Cx^2; 
        K_R(5,2,i)=-M*Cy^2; K_R(2,5,i)=-M*Cy^2;
        K_R(6,3,i)=-M*Cz^2; K_R(3,6,i)=-M*Cz^2;
        K_R(5,1,i)=-M*Cx*Cy; K_R(1,5,i)=-M*Cx*Cy; K_R(4,2,i)=-M*Cx*Cy; K_R(2,4,i)=-M*Cx*Cy;
        K_R(6,1,i)=-M*Cx*Cz; K_R(1,6,i)=-M*Cx*Cz; K_R(4,3,i)=-M*Cx*Cz; K_R(3,4,i)=-M*Cx*Cz;
        K_R(6,2,i)=-M*Cz*Cy; K_R(2,6,i)=-M*Cz*Cy; K_R(5,3,i)=-M*Cz*Cy; K_R(3,5,i)=-M*Cz*Cy;


        % Construção do Vetor LM
        for j=1:6
            if(j <= 3)
                LM(i,j) = restricoes(conectividade(i,2),j+1);
            end
            if(j>3 && j<=6)
                LM(i,j) = restricoes(conectividade(i,3),j-2);   
            end            
        end

        % Matriz de rigidez global
        K_g = zeros(dimensao_K_G,dimensao_K_G);
        for j=1:6
            for k=1:6
                K_g(LM(i,j),LM(i,k))=K_R(j,k,i); 
            end
        end
        K_G = K_G + K_g;
        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PARTE 5
    % Cálculo dos deslocamentos desconhecidos
    K = zeros(num_eq,num_eq);
    F_D = zeros(num_eq,1);

    %Está pegando a parte de K_G e de F relacionadaa aos deslocamentos desconhecidos
    for i=1:num_eq
        for j=1:num_eq
            K(i,j) = K_G(i,j);
        end
        F_D(i,1) = F(i,1);
    end
    d = K\F_D; %solução do sistema de equações

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PARTE 6
    % Cálculo das ações de extremidades das barras
     A_ML = zeros(6,num_el);    % Açoes de engastamento perfeito 
     d_T = zeros(3*num_no,1);   % Deslocamentos totais -> numerados em ordem dos deslocamentos desconhecidos para conhecidos
     d_L = zeros(6,num_el);     % Deslocamento de cada elemento, a partir de d_T -->de acordo com o LM
     C1 = zeros(6,num_el);      % Vetor intermediário que representa R_L*d_L
     A_M = zeros(6,num_el);     % Vetor das ações nas extremidades

    % Adequa o vetor d no vetor de deslocamentos totais d_T
     for i = 1:num_eq
         d_T(i) = d(i);
     end 

    % Atraves do LM adquire-se um vetor d_L de deslocamentos locais 
     for i = 1:num_el
         for j = 1:6
             d_L(j,i) = d_T(LM(i,j));
         end
     end

    % Efetuação do calculo de C1
     for i = 1:num_el
         for j = 1:6
             for k = 1:6
                 C1(j,i) = C1(j,i)+ R_L(j,k,i)*d_L(k,i);
             end
         end     
     end

    % Cálculo do A_M somando as ações de engastamento perfeito
    for i = 1:num_el
         for j = 1:6
             for k = 1:6
                 A_M(j,i) = A_M(j,i)+ K_L(j,k,i)*C1(k,i);
             end
             A_M(j,i)=A_ML(j,i)+A_M(j,i);
         end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Parte 7
    %Calculo das reações de apoio
    aux = (dimensao_K_G-num_eq);
    K_RD = zeros(aux,num_eq);
    F_RL = zeros(aux,1);
    F_R = zeros(aux,1); 
    %Está pegando a parte de K_G e de F relacionadaa aos deslocamentos conhecidos
    for i=num_eq+1:dimensao_K_G 
        for j=1:num_eq
            K_RD(i-num_eq,j) = K_G(i,j);
        end
        F_RL(i-num_eq,1) = F(i,1);
    end

    for i=1:aux 
        for j=1:num_eq
            F_R(i) = F_R(i)+ F_RL(i)+ K_RD(i,j)*d(j);
        end
    end
end

if ((titulo == 'Port')|(titulo == 'Viga'))
    
    % Parte 1
    % Leitura e Armazenamento das Informações de Estrutura
    % Leitura e Armazenamento do número de nós
    num_no = fscanf(entrada,'%i',1); 
    
    % Leitura e Armazenamento do numero de elementos
    num_el = fscanf(entrada,'%i',1); 
    
    % Armazenamento das coordenadas dos nós (nome do no, coordenadas)
    coordenada = fscanf(entrada,'%f',[4,num_no]);  
    coordenada = coordenada'; %Altera o formato da Matriz para o adequado a ser usado em calculos

    % Armazenamento das restrições dos nós(nome do no, restrição)
    restricoes = fscanf(entrada,'%i',[4,num_no]);
    restricoes = restricoes';

    % Armazenamento dos carregamentos dos nós (num do no, Fx, Fy, Mz)
    carregamento_no = fscanf(entrada,'%f',[4,num_no]);
    carregamento_no= carregamento_no';

    % Armazenamento dos carregamentos distribuidos em y (numero do elemento, posiçao
    %de inicio do carregamento(xm), posição final do carregamento(xm),
    %qi_y, qf_y
    carr_disty = fscanf(entrada, '%f', [5,num_el]);
    carr_disty = carr_disty';
    
    % Armazenamento dos carregamentos distribuidos em x (numero do elemento, posiçao
    %de inicio do carregamento(ym), posição final do carregamento(ym),
    %qi_x, qf_x
    carr_distx = fscanf(entrada, '%f', [5,num_el]);
    carr_distx = carr_distx';

    % Conectividades (numero do elemento,no_saida,no_chegada, tipo de material)
    conectividade = fscanf(entrada,'%i',[4,num_el]);
    conectividade = conectividade';

    % Armazenamento das características do material 
    % Tipo de material ,Seção transversal,Mod Elasticidade, Mod Transversal, Momento de inercia
    num_material = fscanf(entrada,'%i',1); %Quantidade de tipos de materiais
    material = fscanf(entrada,'%f',[5,num_material]);
    material = material';
    
    % Consideração apoio elástico
    % num do apoio, nó que ta localizado, valor x, valor y, valor z
    % num_apoios = fscanf(entrada,'%i',1);
    % apoio_elastico = fscanf(entrada,'%f',[5,num_apoios]);  
    % apoio_elastico = apoio_elastico';

    %Fechar arquivo de entrada
    fclose(entrada);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parte 2
    % Renumeração da matriz de restrição
    num_eq = 0;
    for i = 1:num_no
        for j = 2:4
           if(restricoes(i,j) == 1)
               restricoes(i,j)=0;
           else
               restricoes(i,j)= num_eq + 1;
               num_eq = num_eq+1;
           end
        end
    end

    % Renumeração dos nós restritos
    dimensao_K_G = num_eq; %contador que define a dimensão na matriz Global K_G
    for i = 1:num_no
        for j = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
           if(restricoes(i,j)==0)
               restricoes(i,j) = dimensao_K_G+1;
               dimensao_K_G = dimensao_K_G+1;
           end
        end
    end

    for i=1:dimensao_K_G
        for j = 1:num_no
            for k = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
                if(restricoes(j,k) == i)
                    F(i) = carregamento_no(j,k);
                end    
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Parte 3
    % Matriz de Rigidez Global e Vetor LM
    L = zeros(num_el,1);                    % Comprimento do elemento
    R_L = zeros(6,6,num_el);                % Matriz de Rotação
    K_L = zeros(6,6,num_el);                % Matriz de Rigidez no Referencial Local
    K_R = zeros(6,6,num_el);                % Matriz de Rigidez no Referencial Local Rotacionada
    K_g = zeros(dimensao_K_G,dimensao_K_G); % Matriz de Rigidez Global Auxiliar
    K_G = zeros(dimensao_K_G,dimensao_K_G); % Matriz de Rigidez Global
    LM = zeros(num_el,6);                   % Vetor LM
    
    for i = 1:num_el 
        %Abre um loop para cada elemento, a cada iteração será calculado o
        %comprimento do elemento, a matriz de rigidez local, a matriz de
        %rotação, matriz de rigidez rotacionada, o vetor LM, e parte da matriz
        %global. Ao final do loop, a matriz global K_G é obtida.

        % Cálculo do comprimento de cada elemento
        x_i = coordenada(conectividade(i,2),2);
        y_i = coordenada(conectividade(i,2),3);
        z_i = coordenada(conectividade(i,2),4);
        x_f = coordenada(conectividade(i,3),2);
        y_f = coordenada(conectividade(i,3),3);
        z_f = coordenada(conectividade(i,3),4);
        L(i) = sqrt((x_f-x_i)^2+(y_f-y_i)^2+(z_f-z_i)^2);

        % Constantes simplificadoras da matriz rotacionada
        A = (material(conectividade(i,4),2))*(material(conectividade(i,4),3))/L(i);
        B = 12*(material(conectividade(i,4),5))*(material(conectividade(i,4),3))/(L(i))^3;
        C = 6*(material(conectividade(i,4),5))*(material(conectividade(i,4),3))/(L(i))^2;
        D = 4*(material(conectividade(i,4),5))*(material(conectividade(i,4),3))/L(i);

        % Definição da matriz de rigidez local do pórtico plano 
        K_L(1,1,i) = (material(conectividade(i,4),2))*(material(conectividade(i,4),3))/L(i);
        K_L(4,4,i) = K_L(1,1,i);
        K_L(1,4,i) = -K_L(1,1,i);
        K_L(4,1,i) = -K_L(1,1,i);

        K_L(2,2,i) = 12*(material(conectividade(i,4),5))*(material(conectividade(i,4),3))/(L(i))^3;
        K_L(5,5,i) = K_L(2,2,i);
        K_L(2,5,i) = -K_L(2,2,i);
        K_L(5,2,i) = -K_L(2,2,i);

        K_L(3,3,i) = 4*(material(conectividade(i,4),5))*(material(conectividade(i,4),3))/L(i);
        K_L(6,6,i) = K_L(3,3,i);
        K_L(6,3,i)= K_L(3,3,i)/2;
        K_L(3,6,i)= K_L(3,3,i)/2;

        K_L(3,2,i) = 6*(material(conectividade(i,4),5))*(material(conectividade(i,4),3))/(L(i))^2;
        K_L(2,3,i)= K_L(3,2,i);
        K_L(2,6,i)= K_L(3,2,i);
        K_L(6,2,i)= K_L(3,2,i);
        K_L(3,5,i)= -K_L(3,2,i);
        K_L(5,3,i)= -K_L(3,2,i);
        K_L(6,5,i)= -K_L(3,2,i);
        K_L(5,6,i)= -K_L(3,2,i);


        % Definição da matriz de rotação local da treliça espacial    
        Cx = (x_f-x_i)/L(i);
        Cy = (y_f-y_i)/L(i);
        Cz = (z_f-z_i)/L(i);

        %Matriz de rotação para treliça espacial (GERE E WEAVER,1965)  
        R_L(1,1,i) = Cx;                 R_L(4,4,i) = Cx; 
        R_L(2,2,i) = Cx;                 R_L(5,5,i) = R_L(2,2,i);
        R_L(3,3,i) = 1;                  R_L(6,6,i) = R_L(3,3,i); 
        R_L(2,1,i) = -Cy;                R_L(1,2,i) = Cy;
        R_L(5,4,i) = -Cy;                R_L(4,5,i) = Cy;

        % Matriz de rigidez local rotacionada para o referencial global
        K_R(1,1,i) = A*(Cx^2) + B*(Cy^2);       K_R(4,4,i) = K_R(1,1,i);               
        K_R(2,2,i) = A*(Cy^2) + B*(Cx^2);       K_R(5,5,i) = K_R(2,2,i);
        K_R(3,3,i) = D;                         K_R(6,6,i) = K_R(3,3,i);    
        K_R(1,2,i) = (A-B)*Cx*Cy;               K_R(2,1,i) = K_R(1,2,i);
        K_R(1,3,i) = -C*Cy;                     K_R(1,4,i) = -K_R(1,1,i);
        K_R(1,5,i) = -K_R(1,2,i);               K_R(1,6,i) = K_R(1,3,i);
        K_R(2,3,i) = C*Cx;                      K_R(2,4,i) = -K_R(2,1,i);
        K_R(2,5,i) = -K_R(2,2,i);               K_R(2,6,i) = K_R(2,3,i);
        K_R(3,1,i) = K_R(1,3,i);                K_R(3,2,i) = K_R(2,3,i);
        K_R(3,4,i) = -K_R(3,1,i);               K_R(3,5,i) = -K_R(3,2,i);
        K_R(3,6,i) = K_R(3,3,i)/2;              K_R(4,1,i) = -K_R(1,1,i);
        K_R(4,2,i) = -K_R(2,1,i);               K_R(4,3,i) = -K_R(1,3,i);
        K_R(4,5,i) = -K_R(1,5,i);               K_R(4,6,i) = -K_R(1,6,i);
        K_R(5,1,i) = -K_R(1,2,i);               K_R(5,2,i) = -K_R(2,2,i);
        K_R(5,3,i) = -K_R(2,3,i);               K_R(5,4,i) = -K_R(2,4,i);
        K_R(5,6,i) = -K_R(2,6,i);               K_R(6,1,i) = K_R(3,1,i);
        K_R(6,2,i) = K_R(3,2,i);                K_R(6,3,i) = K_R(3,3,i)/2;
        K_R(6,4,i) = K_R(3,4,i);                K_R(6,5,i) = K_R(3,5,i);



        % Construção do Vetor LM
        for j=1:6
            if(j <= 3)
                LM(i,j) = restricoes(conectividade(i,2),j+1);
            end
            if(j>3 && j<=6)
                LM(i,j) = restricoes(conectividade(i,3),j-2);   
            end            
        end

        % Matriz de rigidez global
        K_g = zeros(dimensao_K_G,dimensao_K_G);
        for j=1:6
            for k=1:6
                K_g(LM(i,j),LM(i,k))=K_R(j,k,i); 
            end
        end
        K_G = K_G + K_g; 
    end
    
    % Consideração de apoios elásticos
%     for i=1:num_apoios
%         w = apoio_elastico(i,2);
%         j = apoio_elastico(i,3);
%         k = apoio_elastico(i,4);
%         l = apoio_elastico(i,5);
%         if(j~=0)
%             z = restricoes(w,2);
%             K_G(z,z) = K_G(z,z) +  apoio_elastico(i,3);
%         end
%         if(k~=0)
%             z = restricoes(w,3);
%             K_G(z,z) = K_G(z,z) +  apoio_elastico(i,4);
%         end
%         if(l~=0)
%             z = restricoes(w,4);
%             K_G(z,z) = K_G(z,z) +  apoio_elastico(i,5);
%         end
%     end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parte 4
    % Montagem do vetor de carregamentos externos
    F = zeros(dimensao_K_G,1);
    F_nodal = zeros(dimensao_K_G,1);
    f_disty = zeros(6,num_el); % Vetor de carregamentos distribuidos local 
    F_disty = zeros(dimensao_K_G,1); %Vetor de carregamentos distribuidos global
    f_distx = zeros(6,num_el); % Vetor de carregamentos distribuidos local 
    F_distx = zeros(dimensao_K_G,1); %Vetor de carregamentos distribuidos global
    
    % Montagem do vetor de engastamento perfeito
    for i=1:num_el
        x_i = coordenada(conectividade(i,2),2);
        y_i = coordenada(conectividade(i,2),3);
        z_i = coordenada(conectividade(i,2),4);
        x_f = coordenada(conectividade(i,3),2);
        y_f = coordenada(conectividade(i,3),3);
        z_f = coordenada(conectividade(i,3),4);
        L(i) = sqrt((x_f-x_i)^2+(y_f-y_i)^2+(z_f-z_i)^2);
        
        Cx = (x_f-x_i)/L(i);
        Cy = (y_f-y_i)/L(i);
        Cz = (z_f-z_i)/L(i);
        
        if (Cy == 0)
            f_disty(2,i) = L(i)/2*(carr_disty(i,4) + 3/10*(carr_disty(i,5)-carr_disty(i,4)));
            f_disty(5,i) = L(i)/2*(carr_disty(i,4) + 7/10*(carr_disty(i,5)-carr_disty(i,4)));
            f_disty(3,i) = (L(i)^2)/3*((carr_disty(i,4)/4)+(carr_disty(i,5) - carr_disty(i,4))/10);
            f_disty(6,i) = -(L(i)^2)/4*((carr_disty(i,4)/3)+(carr_disty(i,5) - carr_disty(i,4))/5);

            F_dist_auxy = zeros(dimensao_K_G,1); 
            for j=1:6
                F_dist_auxy(LM(i,j))=f_disty(j,i);
            end
            F_disty = F_disty + F_dist_auxy;
        end
        if (Cx == 0)
            f_distx(1,i) = L(i)/2*(carr_distx(i,4) + 3/10*(carr_distx(i,5)-carr_distx(i,4)));
            f_distx(4,i) = L(i)/2*(carr_distx(i,4) + 7/10*(carr_distx(i,5)-carr_distx(i,4)));
            f_distx(3,i) = (L(i)^2)/3*((carr_distx(i,4)/4)+(carr_distx(i,5) - carr_distx(i,4))/10);
            f_distx(6,i) = -(L(i)^2)/4*((carr_distx(i,4)/3)+(carr_distx(i,5) - carr_distx(i,4))/5);

            F_dist_auxx = zeros(dimensao_K_G,1); 
            for j=1:6
                F_dist_auxx(LM(i,j))=f_distx(j,i);
            end
            F_distx = F_distx + F_dist_auxx;
        end
    end

    % Montagem do vetor de cargas Nodais
    for i=1:dimensao_K_G
        for j = 1:num_no
            for k = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
                if(restricoes(j,k) == i)
                    F_nodal(i) = carregamento_no(j,k);
                end    
            end
        end
    end

    % Construção do vetor de cargas global
    F = F_nodal - F_disty - F_distx; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PARTE 5
    % Cálculo dos deslocamentos desconhecidos
    K = zeros(num_eq,num_eq);
    F_D = zeros(num_eq,1);

    %Está pegando a parte de K_G e de F relacionadaa aos deslocamentos desconhecidos
    for i=1:num_eq
        for j=1:num_eq
            K(i,j) = K_G(i,j);
        end
        F_D(i,1) = F(i,1);
    end
    d = K\F_D; %solução do sistema de equações

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PARTE 6
    % Cálculo das ações de extremidades das barras
     A_ML = zeros(6,num_el); %Açoes de engastamento perfeito 
     d_T = zeros(3*num_no,1); %Deslocamentos totais -> numerados em ordem dos deslocamentos desconhecidos para conhecidos
     d_L = zeros(6,num_el); %Deslocamento de cada elemento, a partir de d_T -->de acordo com o LM
     C1 = zeros(6,num_el); %Vetor intermediário que representa R_L*d_L
     A_M = zeros(6,num_el); %Vetor das ações nas extremidades
     C2 = zeros(6,num_el); %Vetor intermediário que representa R_L*d_L
     f_distx2 = zeros(6,num_el);
     
     for i=1:num_el
            y_i = coordenada(conectividade(i,2),3);
            y_f = coordenada(conectividade(i,3),3);
            Cy = (y_f-y_i)/L(i);

         if(Cy<0)         
            f_distx2(1,i)=f_distx(2,i);
            f_distx2(2,i)=f_distx(1,i);
            f_distx2(4,i)=f_distx(5,i);
            f_distx2(5,i)=f_distx(4,i);
            f_distx2(3,i)=f_distx(3,i);
            f_distx2(6,i)=f_distx(6,i);
         end
        if(Cy>0)         
            f_distx2(1,i)=f_distx(2,i);
            f_distx2(2,i)=-f_distx(1,i);
            f_distx2(4,i)=f_distx(5,i);
            f_distx2(5,i)=-f_distx(4,i);
            f_distx2(3,i)=f_distx(3,i);
            f_distx2(6,i)=f_distx(6,i);
        end
     end
     
     for i=1:num_el
         for j=1:6
            A_ML(j,i) = f_disty(j,i) + f_distx2(j,i);
         end
     end

    % Adequa o vetor d no vetor de deslocamentos totais d_T
     for i = 1:num_eq
         d_T(i) = d(i);
     end 

    % Atraves do LM adquire-se um vetor d_L de deslocamentos locais 
     for i = 1:num_el
         for j = 1:6
             d_L(j,i) = d_T(LM(i,j));
         end
     end

    % Efetuação do calculo de C1
     for i = 1:num_el
         for j = 1:6
             for k = 1:6
                 C1(j,i) = C1(j,i)+ R_L(j,k,i)*d_L(k,i);
             end
         end     
     end

   % Cálculo do A_M somando as ações de engastamento perfeito
    for i = 1:num_el
         for j = 1:6
             for k = 1:6
                 C2(j,i) = C2(j,i)+ K_L(j,k,i)*C1(k,i);
             end
             A_M(j,i)=A_ML(j,i)+A_M(j,i)+C2(j,i);
         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   % Parte 7
   % Calculo das reações de apoio
    aux = (dimensao_K_G-num_eq);
    K_RD = zeros(aux,num_eq);
    F_RL = zeros(aux,1);
    F_R = zeros(aux,1); 

   % Está pegando a parte de K_G e de F relacionadaa aos deslocamentos conhecidos
    for i=num_eq+1:dimensao_K_G 
        for j=1:num_eq
            K_RD(i-num_eq,j) = K_G(i,j);
        end
        F_RL(i-num_eq,1) = F_disty(i,1)+ F_distx(i,1);
    end
    AUX=zeros(aux,1);
    for i=1:aux 
        for j=1:num_eq
            AUX(i) = AUX(i) + K_RD(i,j)*d(j);
        end
        F_R(i) = F_R(i)+ F_RL(i)+ AUX(i);
    end    
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
