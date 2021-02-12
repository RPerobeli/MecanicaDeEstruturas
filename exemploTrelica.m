clear all
close all
clc

%Exemplo de método dos deslocamentos para disciplina Mecânica das estruturas

%Exemplo de trelica

% Abrir o aqrquivo txt de entrada
entrada = fopen('G:\Meu Drive\MESTRADO\Trimestre 3\Mecanica das Estruturas\trelica.txt','r');
% Leitura e Armazenamento do Titulo (no caso, viga)
titulo = fscanf(entrada,'%s',1);

%% Parte 1

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
    
%% Parte 2
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
    
%% Parte 3

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
    
%% Parte 4
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
        
%         K_R(:,:,i) = R_L(:,:,i)*K_L(:,:,i);

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
%% PARTE 5
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

format long
disp(d);

