clear all
close all
clc

%Exemplo de método dos deslocamentos para disciplina Mecânica das estruturas

%Exemplo de Viga

% Abrir o aqrquivo txt de entrada
entrada = fopen('G:\Meu Drive\MESTRADO\Trimestre 3\Mecanica das Estruturas\Viga.txt','r');
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

%% Parte 2

% Renumeração da matriz de restrição
% A renumeração é necessária pois os nós desconhecidos devem ser os
% primeiros nós a serem numerados. para posterior solução
% se a deslocabilidades desconhecidas são d5 e d6, elas serao numeradas
% como d1 e d2.

    num_eq = 0;
    for i = 1:num_no
        for j = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
           if(restricoes(i,j) == 1)
               restricoes(i,j)=0;
           else
               num_eq = num_eq+1;
               restricoes(i,j)= num_eq;
           end
        end
    end

    % Renumeração dos nós restritos
    dimensao_K_G = num_eq; %contador que define a dimensão na matriz Global K_G
    for i = 1:num_no
        for j = 2:4 %Começa em 2 pq a primeira coluna representa o numero do nó
           if(restricoes(i,j)==0)
               dimensao_K_G = dimensao_K_G+1;
               restricoes(i,j) = dimensao_K_G;
           end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
 %% Parte 3
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


        % Definição da matriz de rigidez local do pórtico plano TABELADO
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
        Cx = (x_f-x_i)/L(i); %cossenos diretores
        Cy = (y_f-y_i)/L(i);
        Cz = (z_f-z_i)/L(i);

        %Matriz de rotação para portico plano (GERE E WEAVER,1965)  
        R_L(1,1,i) = Cx;                 R_L(4,4,i) = Cx; 
        R_L(2,2,i) = Cx;                 R_L(5,5,i) = R_L(2,2,i);
        R_L(3,3,i) = 1;                  R_L(6,6,i) = R_L(3,3,i); 
        R_L(2,1,i) = -Cy;                R_L(1,2,i) = Cy;
        R_L(5,4,i) = -Cy;                R_L(4,5,i) = Cy;

        K_R(:,:,i) = R_L(:,:,i)*K_L(:,:,i);


        % Construção do Vetor LM
        % O vetor LM é o responsável por modificar a matriz de rigidez
        % global, de forma que os deslocamentos fiquem nas posições
        % corretas, uma vez que não é necessário saber todas as 
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
 
%% Parte 4
    % Montagem do vetor de carregamentos externos
    F = zeros(dimensao_K_G,1);
    F_nodal = zeros(dimensao_K_G,1);
    f_disty = zeros(6,num_el); % Vetor de carregamentos distribuidos local em y
    F_disty = zeros(dimensao_K_G,1); %Vetor de carregamentos distribuidos global em y
    f_distx = zeros(6,num_el); % Vetor de carregamentos distribuidos local em x
    F_distx = zeros(dimensao_K_G,1); %Vetor de carregamentos distribuidos global em y
    
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