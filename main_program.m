
% Exercício Programa - Introdução a Sistemas Elétricos de Potência

% Versão : 19.06.16.2
% Data da última edição: 16/06/2019 - 24:00

% Autores:  Gustavo Gransotto Ribeiro       9300557
%           Pedro Emanuel Rodrigues Castro  98
%           Pedro Nunes Tadeu Rodrigues     90


% Limpa variáveis anteriores
clear all;
clc;

% Inicia o timer para medir desempenho do programa
tic();

% Definição de variáveis
global Vnom;
global Vth;
global Zth;
global Ith
global Yth;
global Rmax;
global elapsed_time;
global a1;
global a2;
global Vg;
global Ig;
% Operador alpha (a1) e alpha quadrado (a2)
a1 = 1*exp(deg2rad(120)*1i);
a2 = 1*exp(deg2rad(-120)*1i);

% Tensão nominal de linha
Vnom = 13800;
% Fonte de tensao do equivalente de Thevenin do no da subestacao
Vth = Vnom / sqrt (3) ;
Zth = (Vth^2) /(((1/3)*10E8 ) *( cos (80*pi/180) - 1i* sin (80*pi/180) )); % Impedancia equivalente do Thevenin do no da subestacao
Ith = Vth / Zth ; % Fonte de corrente do equivalente de Norton do no da subestacao
Yth = 1 / Zth ; % Admitancia do equivalente de Norton do no da subestacao
Rmax = 10; % Valor maximo da resistencia de falta

Vg(1,1) = Vth;
Vg(2,1) = Vth * a2;
Vg(3,1) = Vth * a1;

Ig(1,1) = Vg(1,1) / Zth;
Ig(2,1) = Vg(2,1) / Zth;
Ig(3,1) = Vg(3,1) / Zth;


% Carrega os arquivos do enunciado

% Arquivo com a matriz 'cargas' = <no onde a carga esta conectada >< potencia [kW] da carga >< fator de potencia ( indutivo )>
CAR039;


% Arquivo com a matriz 'topologia' = <no pai ><no filho >< comprimento [m]>< resistencia do trecho [ ohms /m]>< reatancia do trecho [ ohms /m]>
TOP039;

% Arquivo com a matriz 'Emedido' = <numero do evento de curto - circuito >< parte real da tensao [V]>< parte imaginaria da tensao [V]>
VOL039;


fprintf('Caso , Nó 1 , Nó 2, Dist , Res , Funcao\n');

topologiaBackup = topologia;

% Calcula impedâncias médias próprias e mútuas no trecho a partir do arquivo de topologia (Só serão usadas nos 5 primeiros casos de simulação)

% ZtopMed = < no pai > < no filho > < comprimento [m] > < impedancia propria media no trecho [ohms/m] > < impedancia mutua media no trecho [ohms/m] >
tamanhoTOP = size(topologia,1);

Resultados_Simulacao = zeros(10 * tamanhoTOP,6);
Resultados_Localizacao = zeros(10,6);

ZtopMed = zeros(tamanhoTOP,5);

for i = 1:tamanhoTOP
    ZtopMed(i,1) = topologia(i,1);
    ZtopMed(i,2) = topologia(i,2);
    ZtopMed(i,3) = topologia(i,3);
    % Imp. própria média = Zaa + Zbb + Zcc / 3
    ZtopMed(i,4) = ( (topologia(i,4) + topologia(i,12)+ topologia(i,20)) + 1i *(topologia(i,5) + topologia(i,13) + topologia(i,21)) ) / 3;

    % Imp. mutua média = Zab + Zac + Zba + Zbc + Zca + Zcb / 6
    ZtopMed(i,5) = (topologia(i,6) + topologia(i,8)+ topologia(i,10)) + 1i *(topologia(i,7) + topologia(i,9) + topologia(i,11));
    ZtopMed(i,5) = ZtopMed(i,5) + (topologia(i,14) + topologia(i,16)+ topologia(i,18)) + 1i *(topologia(i,15) + topologia(i,17) + topologia(i,19));
    ZtopMed(i,5) = ZtopMed(i,5) / 6;
endfor

ZtopMedBackup = ZtopMed;
ZtopMed(tamanhoTOP+1,1) = 1000; % Coloca na nova linha a falta como nó de ligação 1

% Cria a matriz de admitancias Primitivas Ypr contendo a impedancia própria e mútua nos trechos (CASO EQUILIBRADO)

Ypr = zeros(tamanhoTOP+1,tamanhoTOP+1);
for i = 1 : tamanhoTOP
    % No caso de circuitos trifásicos equilibrados, o equivalente monofásico da rede considera a subtração da impedância mútua média pois a1 + a2 = -1
    Ypr(i,i) = 1 / ( ZtopMed(i,3) * ( ZtopMed(i,4) - ZtopMed(i,5) ) );
endfor
YprBackup = Ypr;

% Calcula a impedância na carga a partir da tensão nominal da linha
tamanhoCAR = size(cargas,1);
Zcarga = zeros(tamanhoCAR+1,2);
for i = 1 : tamanhoCAR
    Zcarga(i,1) = cargas(i,1);
endfor
Zcarga(tamanhoCAR+1,1) = 1000;

Zcarga = sort(Zcarga);

for i = 1 : tamanhoCAR
    for j = 1 : tamanhoCAR
      if (Zcarga(i,1) == cargas(j,1))
        Zcarga(i,2) = (Vnom^2) /(1000*(cargas (j ,2) - 1i*cargas(j ,2).*tan ( acos ( cargas (j,3) ) ) ) ) ;
      endif
    endfor
endfor

% Último espaço do vetor é reservado para a impedância da carga


% Cria um vetor com todos os nós distintos do problema
nosDistintos(1) = topologiaBackup(1,1);

for w = 1:2

    for i = 2:tamanhoTOP
        numDeNosDistintos = size(nosDistintos,2);
        flag = 0;

        for j = 1 : numDeNosDistintos

            if (nosDistintos(j) == topologiaBackup(i,w))
                flag = 1;
                break;
            endif

        endfor

        if (flag == 0)
            nosDistintos(numDeNosDistintos+1) = topologiaBackup(i,w);

        endif

    endfor

endfor
numDeNosDistintos = size(nosDistintos,2);
nosDistintos(numDeNosDistintos+1) = 1000;
nosDistintos = sort(nosDistintos);


% Cria matriz de incidências com o tamanho da matriz de Nós distintos + 1 (Nó da falta)
% Uma forma rápida de checar se essa matriz está certa é somar todos os seus elementos. Deve ser 2*numDeLinhas da matriz topologia

Matriz_Incidencias = zeros(tamanhoTOP+1, numDeNosDistintos+1);

for i = 1 : tamanhoTOP

    for j = 1 : numDeNosDistintos

        if ( ZtopMed(i,1) == nosDistintos(j) )

            Matriz_Incidencias(i,j) = 1; % Corrente no ramo i sai do nó j

        elseif ( ZtopMed(i,2) == nosDistintos(j) )

            Matriz_Incidencias(i,j) = -1; % Corrente no ramo i entra no nó j

        else

            Matriz_Incidencias(i,j) = 0; % Ramo i não está conectado ao nó j

        endif

    endfor

endfor

Matriz_Incidencias_Backup = Matriz_Incidencias;

Inos = zeros( numDeNosDistintos+1 , 1 );
Inos(1,1) = Ith;

% for casoSimulacao = 1 : 5
%     funcaoLocaliz = Inf;
%     distFaltaLocaliz = 0;
%     ResFaltaLocaliz = 0;
%     no_1_Localiz = 0;
%     no_2_Localiz = 0;
%
%     for trechoDaFalta = 1 : tamanhoTOP
%
%         noDeLigacao_1_Falta = topologiaBackup(trechoDaFalta, 1);   % Nó
%         noDeLigacao_2_Falta = topologiaBackup(trechoDaFalta, 2);   % Nó
%         funcao = Inf;
%         distFaltaCalc = 0;
%         ResFaltaCalc = 0;
%
%
%         for i = 2 : size(ZtopMedBackup,2)
%             ZtopMed(tamanhoTOP+1,i) = ZtopMedBackup(trechoDaFalta,i);
%         endfor
%
%         ZtopMed(trechoDaFalta,2) = 1000; % Coloca no lugar do nó de ligação 2 original do trecho o nó da falta
%
%
%         for distancia_no_1_Falta = 1 : topologiaBackup(trechoDaFalta,3)-1 % [m]
%
%             ZtopMed(tamanhoTOP+1,3) = ZtopMedBackup(trechoDaFalta,3) - distancia_no_1_Falta;
%             ZtopMed(trechoDaFalta,3) = distancia_no_1_Falta;
%
%             Ypr (trechoDaFalta,trechoDaFalta) =  1 / ( ZtopMed(trechoDaFalta,3) * ( ZtopMed(trechoDaFalta,4) - ZtopMed(trechoDaFalta,5) ) );
%             Ypr (tamanhoTOP+1,tamanhoTOP+1) =  1 / ( ZtopMed(tamanhoTOP+1,3) * ( ZtopMed(tamanhoTOP+1,4) - ZtopMed(tamanhoTOP+1,5) ) );
%
%             for j = 1 : numDeNosDistintos+1
%
%                     if ( ZtopMed(trechoDaFalta,1) == nosDistintos(j) )
%
%                         Matriz_Incidencias(trechoDaFalta,j) = 1; % Corrente no ramo da falta sai do nó j
%
%                     elseif ( ZtopMed(trechoDaFalta,2) == nosDistintos(j) )
%
%                         Matriz_Incidencias(trechoDaFalta,j) = -1; % Corrente no ramo da falta entra no nó j
%
%                     else
%
%                         Matriz_Incidencias(trechoDaFalta,j) = 0; % Ramo da falta não está conectado ao nó j
%
%                     endif
%
%                     if ( ZtopMed(tamanhoTOP+1, 1) == nosDistintos(j) )
%
%                         Matriz_Incidencias(tamanhoTOP+1,j) = 1; % Corrente no ramo da falta sai do nó j
%
%
%                     elseif ( ZtopMed(tamanhoTOP+1, 2) == nosDistintos(j) )
%
%                         Matriz_Incidencias(tamanhoTOP+1, j) = -1; % Corrente no ramo da falta entra no nó j
%
%
%                     else
%
%                         Matriz_Incidencias(tamanhoTOP+1,j) = 0; % Ramo da falta não está conectado ao nó j
%
%
%                     endif
%             endfor
%
%
%             % Cria a matriz de admitâncias nodais inserindo as admitâncias da linha
%             Ynos = transpose(Matriz_Incidencias) * Ypr * Matriz_Incidencias;
%
%             Ynos(1,1) = Ynos(1,1) + Yth; % Insere a admitância equivalente de Thevenin
%             YnosBackup = Ynos;
%
%             for i = 1 : numDeNosDistintos-1
%                 Ynos(i+1,i+1) = Ynos(i+1,i+1) + ( 1 / Zcarga(i,2) );
%             endfor
%
%             for resistenciaDaFalta = 0.1 : 0.1 : Rmax  % [ohms]
%
%                 Zcarga(tamanhoCAR+1,2) = resistenciaDaFalta;
%
%                 Ynos(numDeNosDistintos+1, numDeNosDistintos+1) = YnosBackup(numDeNosDistintos+1, numDeNosDistintos+1) + ( 1 / Zcarga(tamanhoCAR+1,2) );
%
%                 % Calcula a tensão nos nós a partir da matriz de admitâncias e da corrente de thevenin calculada
%
%                 Ecalc = inv(Ynos) * Inos;
%                 E10calc = Ecalc(1);
%                 E10med = Emedido(casoSimulacao,2) + 1i * Emedido(casoSimulacao,3); % Pega somente caso de simulação 1 e tensões de fase em A
%
%
%                 funcao_old = abs ( E10med - E10calc ) / abs ( E10med ) ;
%                 if funcao_old < funcao
%                   distFaltaCalc = distancia_no_1_Falta ;
%                   ResFaltaCalc = resistenciaDaFalta ;
%                   funcao = funcao_old ;
%                 endif
%
%             endfor % resistenciaDaFalta = 0.1 : 0.1 : Rmax  % [ohms]
%
%         endfor % distancia_no_1_Falta = 1 : topologiaBackup(trechoDaFalta,3)-1 % [m]
%
%         fprintf('%02.f , %03.f , %03.f , %03.f , %2.1f , %2.3f\n', casoSimulacao, noDeLigacao_1_Falta, noDeLigacao_2_Falta, distFaltaCalc, ResFaltaCalc, funcao);
%         Resultados_Simulacao((casoSimulacao-1)*tamanhoTOP + trechoDaFalta, : ) = [casoSimulacao, noDeLigacao_1_Falta, noDeLigacao_2_Falta, distFaltaCalc, ResFaltaCalc, funcao];
%         OUT_ID = fopen('OUT039.csv','a+');
%         fprintf(OUT_ID,'%02.f, %03.f, %03.f, %03.f, %2.1f, %2.3f\n', casoSimulacao, noDeLigacao_1_Falta, noDeLigacao_2_Falta, distFaltaCalc, ResFaltaCalc, funcao);
%         fclose(OUT_ID);
%
%         if funcao < funcaoLocaliz
%           distFaltaLocaliz = distFaltaCalc;
%           ResFaltaLocaliz = ResFaltaCalc;
%           no_1_Localiz = noDeLigacao_1_Falta;
%           no_2_Localiz = noDeLigacao_2_Falta;
%           funcaoLocaliz = funcao;
%         endif
%
%         for i = 1 : size(ZtopMed,2)
%             ZtopMed(trechoDaFalta,i) = ZtopMedBackup(trechoDaFalta,i); %Reseta a linha alterada da matriz ZtopMed com o backup
%         endfor
%
%         Ypr(trechoDaFalta, : ) = YprBackup(trechoDaFalta, : );
%         Matriz_Incidencias(trechoDaFalta, : ) = Matriz_Incidencias_Backup (trechoDaFalta, : );
%
%
%     endfor % trechoDaFalta = 1 : size(topologiaBackup,1)
%
%     Resultados_Localizacao ( casoSimulacao, : ) = [casoSimulacao, no_1_Localiz, no_2_Localiz, distFaltaLocaliz, ResFaltaLocaliz, funcaoLocaliz];
%     REL_ID = fopen('REL039.csv','a+');
%         fprintf(REL_ID,'%02.f, %03.f, %03.f, %03.f, %2.1f, %2.3f\n',casoSimulacao, no_1_Localiz, no_2_Localiz, distFaltaLocaliz, ResFaltaLocaliz, funcaoLocaliz);
%     fclose(REL_ID);
% endfor % casoSimulacao = 1 : 1 % 5




%%%%%%%%%%%%%%%%%%%%%%%%%% CASOS DE SIMULAÇÃO COM REDE DESEQUILIBRADA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Cria a matriz de admitancias Primitivas Ypr contendo a impedancia própria e mútua nos trechos (CASO TRIFÁSICO DESEQUILIBRADO)

Ypr = zeros(3*(tamanhoTOP+1),3*(tamanhoTOP+1));

for i = 0 : tamanhoTOP-1
  j = 1 + 3*i;

  Ypr(j,j) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,4) + 1i * topologiaBackup(i+1,5) ) ); % Insere admitância própria da fase A
  Ypr(j,j+1) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,6) + 1i * topologiaBackup(i+1,7) )); % Insere admitância mútua da fase A em relação à B
  Ypr(j,j+2) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,8) + 1i * topologiaBackup(i+1,9) )); % Insere admitância mútua da fase A em relação à C

  Ypr(j+1,j) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,10) + 1i * topologiaBackup(i+1,11) ) ); % Insere admitância mútua da fase B em relação à A
  Ypr(j+1,j+1) = 1 / (  topologiaBackup(i+1,3) * ( topologiaBackup(i+1,12) + 1i * topologiaBackup(i+1,13) ) ); % Insere admitância própria da fase B
  Ypr(j+1,j+2) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,14) + 1i * topologiaBackup(i+1,15) ) ); % Insere admitância mútua da fase B em relação à C

  Ypr(j+2,j) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,16) + 1i * topologiaBackup(i+1,17) ) ); % Insere admitância mútua da fase C em relação à A
  Ypr(j+2,j+1) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,18) + 1i * topologiaBackup(i+1,19) ) ); % Insere admitância mútua da fase C em relação à B
  Ypr(j+2,j+2) = 1 / ( topologiaBackup(i+1,3) * ( topologiaBackup(i+1,20) + 1i * topologiaBackup(i+1,21) ) ); % Insere admitância própria da fase C

endfor

YprBackup = Ypr;

% Cria matriz de Incidências (CASO TRIFÁSICO DESEQUILIBRADO)

% Cria matriz de incidências com o tamanho da matriz de Nós distintos + 1 (Nó da falta)
% Uma forma rápida de checar se essa matriz está certa é somar todos os seus elementos. Deve ser 2*numDeLinhas da matriz topologia

Matriz_Incidencias = zeros(3*(tamanhoTOP+1), 3*(numDeNosDistintos+1));

for i = 0 : tamanhoTOP-1

    k =  1 + 3*i;

    for j = 0 : numDeNosDistintos-1

        w = 1 + 3*j;

        if ( topologiaBackup(i+1,1) == nosDistintos(j+1) )

            Matriz_Incidencias(k,w) = 1; % Corrente no ramo i sai do nó j
            Matriz_Incidencias(k+1,w+1) = 1; % Corrente no ramo i sai do nó j
            Matriz_Incidencias(k+2,2+2) = 1; % Corrente no ramo i sai do nó j

        elseif ( topologiaBackup(i+1,2) == nosDistintos(j+1) )

            Matriz_Incidencias(k,w) = -1; % Corrente no ramo i entra no nó j
            Matriz_Incidencias(k+1,w+1) = -1; % Corrente no ramo i entra no nó j
            Matriz_Incidencias(k+2,w+2) = -1; % Corrente no ramo i entra no nó j

        endif

    endfor

endfor

Matriz_Incidencias_Backup = Matriz_Incidencias;

Inos = zeros( 3*(numDeNosDistintos+1), 1 );
Inos(1,1) = Ig(1,1);
Inos(2,1) = Ig(2,1);
Inos(3,1) = Ig(3,1);

for casoSimulacao = 6 : 6
    funcaoLocaliz = Inf;
    distFaltaLocaliz = 0;
    ResFaltaLocaliz = 0;
    no_1_Localiz = 0;
    no_2_Localiz = 0;

    for trechoDaFalta = 0 : tamanhoTOP-1

        trechoDaFaltaTri = 1 + 3*trechoDaFalta;
        noDeLigacao_1_Falta = topologiaBackup(trechoDaFalta+1, 1);   % Nó
        noDeLigacao_2_Falta = topologiaBackup(trechoDaFalta+1, 2);   % Nó
        funcao = Inf;
        distFaltaCalc = 0;
        ResFaltaCalc = 0;


        for i = 2 : size(topologiaBackup,2)
            topologia(tamanhoTOP+1,i) = topologiaBackup(trechoDaFalta+1,i);
        endfor

        topologia(trechoDaFalta+1,2) = 1000; % Coloca no lugar do nó de ligação 2 original do trecho o nó da falta


        for distancia_no_1_Falta = 1 : topologiaBackup(trechoDaFalta+1,3)-1 % [m]

            topologia(tamanhoTOP+1,3) = topologiaBackup(trechoDaFalta+1,3) - distancia_no_1_Falta;
            topologia(trechoDaFalta+1,3) = distancia_no_1_Falta;

            % Reescreve a impedância primitiva no trecho em que a falta foi inserida
            Ypr(trechoDaFaltaTri,trechoDaFaltaTri) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,4) + 1i * topologia(trechoDaFalta+1,5) ) ); % Insere admitância própria da fase A
            Ypr(trechoDaFaltaTri,trechoDaFaltaTri +1) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,6) + 1i * topologia(trechoDaFalta+1,7) )); % Insere admitância mútua da fase A em relação à B
            Ypr(trechoDaFaltaTri,trechoDaFaltaTri +2) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,8) + 1i * topologia(trechoDaFalta+1,9) )); % Insere admitância mútua da fase A em relação à C

            Ypr(trechoDaFaltaTri +1,trechoDaFaltaTri +1) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,10) + 1i * topologia(trechoDaFalta+1,11) )); % Insere admitância mútua da fase B em relação à A
            Ypr(trechoDaFaltaTri +1,trechoDaFaltaTri +1) = 1 / (  topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,12) + 1i * topologia(trechoDaFalta+1,13) )); % Insere admitância própria da fase B
            Ypr(trechoDaFaltaTri +1,trechoDaFaltaTri +2) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,14) + 1i * topologia(trechoDaFalta+1,15) )); % Insere admitância mútua da fase B em relação à C

            Ypr(trechoDaFaltaTri +2,trechoDaFaltaTri +1) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,16) + 1i * topologia(trechoDaFalta+1,17) )); % Insere admitância mútua da fase C em relação à A
            Ypr(trechoDaFaltaTri +2,trechoDaFaltaTri +1) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,18) + 1i * topologia(trechoDaFalta+1,19) )); % Insere admitância mútua da fase C em relação à B
            Ypr(trechoDaFaltaTri +2,trechoDaFaltaTri +2) = 1 / ( topologia(trechoDaFalta+1,3) * ( topologia(trechoDaFalta+1,20) + 1i * topologia(trechoDaFalta+1,21) )); % Insere admitância própria da fase C

            % Reescreve a impedância primitiva no trecho da falta
            Ypr(3*tamanhoTOP +1, 3*tamanhoTOP +1) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,4) + 1i * topologia(tamanhoTOP+1,5) ) ); % Insere admitância própria da fase A
            Ypr(3*tamanhoTOP +1, 3*tamanhoTOP +2) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,6) + 1i * topologia(tamanhoTOP+1,7) )); % Insere admitância mútua da fase A em relação à B
            Ypr(3*tamanhoTOP +1, 3*tamanhoTOP +3) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,8) + 1i * topologia(tamanhoTOP+1,9) )); % Insere admitância mútua da fase A em relação à C

            Ypr(3*tamanhoTOP +2, 3*tamanhoTOP +1) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,10) + 1i * topologia(tamanhoTOP+1,11) )); % Insere admitância mútua da fase B em relação à A
            Ypr(3*tamanhoTOP +2, 3*tamanhoTOP +2) = 1 / (  topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,12) + 1i * topologia(tamanhoTOP+1,13) )); % Insere admitância própria da fase B
            Ypr(3*tamanhoTOP +2, 3*tamanhoTOP +3) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,14) + 1i * topologia(tamanhoTOP+1,15) )); % Insere admitância mútua da fase B em relação à C

            Ypr(3*tamanhoTOP+3, 3*tamanhoTOP+1) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,16) + 1i * topologia(tamanhoTOP+1,17) )); % Insere admitância mútua da fase C em relação à A
            Ypr(3*tamanhoTOP+3, 3*tamanhoTOP+2) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,18) + 1i * topologia(tamanhoTOP+1,19) )); % Insere admitância mútua da fase C em relação à B
            Ypr(3*tamanhoTOP+3, 3*tamanhoTOP+3) = 1 / ( topologia(tamanhoTOP+1,3) * ( topologia(tamanhoTOP+1,20) + 1i * topologia(tamanhoTOP+1,21) )); % Insere admitância própria da fase C


            for j = 0 : numDeNosDistintos

                    w = 1 + 3*j;

                    if ( topologia(trechoDaFalta+1,1) == nosDistintos(j+1) )

                        Matriz_Incidencias(trechoDaFaltaTri , w) = 1; % Corrente no ramo da falta sai do nó j
                        Matriz_Incidencias(trechoDaFaltaTri +1 , w+1) = 1; % Corrente no ramo da falta sai do nó j
                        Matriz_Incidencias(trechoDaFaltaTri +2 , w+2) = 1; % Corrente no ramo da falta sai do nó j

                    elseif ( topologia(trechoDaFalta+1,2) == nosDistintos(j+1) )

                        Matriz_Incidencias(trechoDaFaltaTri , w) = -1; % Corrente no ramo da falta entra no nó j
                        Matriz_Incidencias(trechoDaFaltaTri +1 , w+1) = -1; % Corrente no ramo da falta entra no nó j
                        Matriz_Incidencias(trechoDaFaltaTri +2 , w+2) = -1; % Corrente no ramo da falta entra no nó j

                    else

                        Matriz_Incidencias(trechoDaFaltaTri , w) = 0; % Ramo da falta não está conectado ao nó j
                        Matriz_Incidencias(trechoDaFaltaTri +1 , w+1) = 0; % Ramo da falta não está conectado ao nó j
                        Matriz_Incidencias(trechoDaFaltaTri +2 , w+2) = 0; % Ramo da falta não está conectado ao nó j

                    endif

                    if ( topologia(tamanhoTOP+1, 1) == nosDistintos(j+1) )

                        Matriz_Incidencias(3*tamanhoTOP +1 , w) = 1; % Corrente no ramo da falta sai do nó j
                        Matriz_Incidencias(3*tamanhoTOP +2 , w+1) = 1; % Corrente no ramo da falta sai do nó j
                        Matriz_Incidencias(3*tamanhoTOP +3 , w+2) = 1; % Corrente no ramo da falta sai do nó j

                    elseif ( topologia(tamanhoTOP+1, 2) == nosDistintos(j+1) )

                        Matriz_Incidencias(3*tamanhoTOP +1 , w) = -1; % Corrente no ramo da falta entra no nó j
                        Matriz_Incidencias(3*tamanhoTOP +2 , w+1) = -1; % Corrente no ramo da falta entra no nó j
                        Matriz_Incidencias(3*tamanhoTOP +3 , w+2) = -1; % Corrente no ramo da falta entra no nó j

                    else

                        Matriz_Incidencias(3*tamanhoTOP +1 , w) = 0; % Ramo da falta não está conectado ao nó j
                        Matriz_Incidencias(3*tamanhoTOP +2 , w+1) = 0; % Ramo da falta não está conectado ao nó j
                        Matriz_Incidencias(3*tamanhoTOP +3 , w+2) = 0; % Ramo da falta não está conectado ao nó j

                    endif
            endfor


            % Cria a matriz de admitâncias nodais inserindo as admitâncias da linha
            Ynos = transpose(Matriz_Incidencias) * Ypr * Matriz_Incidencias;

            Ynos(1,1) = Ynos(1,1) + Yth; % Insere a admitância equivalente de Thevenin
            Ynos(2,2) = Ynos(2,2) + Yth; % Insere a admitância equivalente de Thevenin
            Ynos(3,3) = Ynos(3,3) + Yth; % Insere a admitância equivalente de Thevenin

            for i = 0 : numDeNosDistintos-2
                k = 1 + 3*i;
                Ynos(k,k) = Ynos(k,k) + ( 1 / Zcarga(i+1,2) );
                Ynos(k+1,k+1) = Ynos(k+1,k+1) + ( 1 / Zcarga(i+1,2) );
                Ynos(k+2,k+2) = Ynos(k+2,k+2) + ( 1 / Zcarga(i+1,2) );
            endfor

            YnosBackup = Ynos;

            for resistenciaDaFalta = 0.1 : 0.1 : Rmax  % [ohms]

                Zcarga(tamanhoCAR+1,2) = resistenciaDaFalta;

                Ynos(3*numDeNosDistintos +1, 3*numDeNosDistintos +1) = YnosBackup(3*numDeNosDistintos +1, 3*numDeNosDistintos +1) + ( 1 / Zcarga(tamanhoCAR+1,2) );
                Ynos(3*numDeNosDistintos +2, 3*numDeNosDistintos +2) = YnosBackup(3*numDeNosDistintos +2, 3*numDeNosDistintos +2) + ( 1 / Zcarga(tamanhoCAR+1,2) );
                Ynos(3*numDeNosDistintos +3, 3*numDeNosDistintos +3) = YnosBackup(3*numDeNosDistintos +3, 3*numDeNosDistintos +3) + ( 1 / Zcarga(tamanhoCAR+1,2) );
                % Calcula a tensão nos nós a partir da matriz de admitâncias e da corrente de thevenin calculada

                Ecalc = inv(Ynos) * Inos;
                E10calc = [ Ecalc(1), Ecalc(2), Ecalc(3)];
                E10med(1) = Emedido(casoSimulacao,2) + 1i * Emedido(casoSimulacao,3); % Pega somente caso de simulação 1 e tensões de fase em A
                E10med(2) = Emedido(casoSimulacao,4) + 1i * Emedido(casoSimulacao,5); % Pega somente caso de simulação 1 e tensões de fase em B
                E10med(3) = Emedido(casoSimulacao,6) + 1i * Emedido(casoSimulacao,7); % Pega somente caso de simulação 1 e tensões de fase em A

                funcao_old = abs ( E10calc(1) - E10med(1)) / abs ( E10med(1) ) ;
                funcao_old = funcao_old + abs ( E10calc(2) - E10med(2)) / abs ( E10med(2) ) ;
                funcao_old = funcao_old + abs ( E10calc(3) - E10med(3)) / abs ( E10med(3) ) ;

                if funcao_old < funcao
                  distFaltaCalc = distancia_no_1_Falta ;
                  ResFaltaCalc = resistenciaDaFalta ;
                  funcao = funcao_old ;
                endif

            endfor % resistenciaDaFalta = 0.1 : 0.1 : Rmax  % [ohms]

        endfor % distancia_no_1_Falta = 1 : topologiaBackup(trechoDaFalta,3)-1 % [m]

        fprintf('%02.f , %03.f , %03.f , %03.f , %2.1f , %2.3f\n', casoSimulacao, noDeLigacao_1_Falta, noDeLigacao_2_Falta, distFaltaCalc, ResFaltaCalc, funcao);
        %Resultados_Simulacao((casoSimulacao-1)*tamanhoTOP + trechoDaFalta, : ) = [casoSimulacao, noDeLigacao_1_Falta, noDeLigacao_2_Falta, distFaltaCalc, ResFaltaCalc, funcao];
        OUT_ID = fopen('OUT039.csv','a+');
        fprintf(OUT_ID,'%02.f, %03.f, %03.f, %03.f, %2.1f, %2.3f\n', casoSimulacao, noDeLigacao_1_Falta, noDeLigacao_2_Falta, distFaltaCalc, ResFaltaCalc, funcao);
        fclose(OUT_ID);

        if funcao < funcaoLocaliz
          distFaltaLocaliz = distFaltaCalc;
          ResFaltaLocaliz = ResFaltaCalc;
          no_1_Localiz = noDeLigacao_1_Falta;
          no_2_Localiz = noDeLigacao_2_Falta;
          funcaoLocaliz = funcao;
        endif

        for i = 1 : size(topologia,2)
            topologia(trechoDaFalta+1,i) = topologiaBackup(trechoDaFalta+1,i); %Reseta a linha alterada da matriz topologia com o backup
        endfor

        for i = 0 : 2
          Ypr(trechoDaFaltaTri + i, : ) = YprBackup(trechoDaFaltaTri + i, : );
          Matriz_Incidencias(trechoDaFaltaTri + i, : ) = Matriz_Incidencias_Backup (trechoDaFaltaTri + i, : );
        endfor


    endfor % trechoDaFalta = 0 : size(topologiaBackup,1)-1

    %Resultados_Localizacao ( casoSimulacao, : ) = [casoSimulacao, no_1_Localiz, no_2_Localiz, distFaltaLocaliz, ResFaltaLocaliz, funcaoLocaliz];
    REL_ID = fopen('REL039.csv','a+');
    fprintf(REL_ID,'%02.f, %03.f, %03.f, %03.f, %2.1f, %2.3f\n',casoSimulacao, no_1_Localiz, no_2_Localiz, distFaltaLocaliz, ResFaltaLocaliz, funcaoLocaliz);
    fclose(REL_ID);

endfor % casoSimulacao = 6 : 10



elapsed_time = toc();
minutos = fix(elapsed_time/60);
horas = fix(minutos/60);
dias = fix(horas/24);
horas = horas - dias*24;
minutos = minutos - horas*60;
segundos = round(elapsed_time - 60*(minutos +60*(horas + 24*dias)));
fprintf('\n\n Tempo de simulacao(s):\t %2.f d : %2.f h : %2.f m : %2.f s\n\n', dias, horas, minutos, segundos);

clear w;
clear flag;
clear j;
clear k;
clear i;
clear dias;
clear horas;
clear minutos;
clear segundos;
