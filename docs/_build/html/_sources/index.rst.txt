.. Filtering Techniques and Data Analysis of data from CoRoT and Kepler satellites documentation master file, created by
   sphinx-quickstart on Thu Feb 25 13:40:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Study of Filtering Processes of Light Curves's documentation!
==========================================================================

|docs| |last commit| |most used lang|

Resumo
=======

- Dados gerais sobre o projeto:
   1. NSEE, FAPESP
   2. Autoria
   3. Projeto apresenta um pipeline de análise e processamento das curvas de luz

As missões espaciais CoRoT e Kepler observaram juntas mais de 300.000 estrelas e, para cada estrela observada, foi gerada um gráfico que mensura a variabilidade do seu brilho em função do tempo. Esses dados gerados são conhecidos como Curvas de Luz e serão o objeto de estudo desse trabalho. Para cada curva de luz,   há três tipos de classificações possíveis, dependendo dos eclipses (detectados nesses sinais) gerados pela passagem de um corpo celeste entre a estrela observada e o observador. São elas: 

A. Curvas oriundas de **exoplanetas**
B. Curvas oriundas de **binárias eclipsantes**
C. Curvas sem quaisquer eclipses.

Os exoplanetas, detectados nas curvas do tipo A, são definidos como planetas que estão fora do sistema solar e, aplicando uma modelagem teórica sobre essas curvas, pode-se calcular alguns dos parâmetros dos exoplanetas: distância do exoplaneta à estrela que orbita, raio da estrela que orbita, raio do exoplaneta, parâmetro de impacto de passagem do exoplaneta na frente da estrela e perídio orbital do exoplaneta. Já as binárias eclipsantes, detectados nas curvas do tipo B, é definido como um par de estrelas que orbitam o centro de gravidade em comum e que se eclipsam, e possuem outras aplicações no campo da astronomia, mas que não são o foco desse trabalho. Por outro lado, as curvas do tipo C não possuem classificação definidas pois não correspondem a exoplanetas nem a binárias eclipsantes, visto que não foram detectados eclipses. 

Não obstante, as curvas de luz possuem um grande problema que dificulta tanto sua classificação entre as três classes definidas anteriormente quanto a determinação nos parâmetros dos exoplanetas: apresentam ruído de alta frequência. Nesse trabalho, propomos averiguar, usando os dados disponíveis pelos telescópios espaciais CoRoT e Kepler, quais as melhores técnicas de filtragem capazes de remover o ruído de alta frequência das curvas de luz, sem comprometer sua classificação nas três classes definidas e que melhorarão a precisão obtida ao se determinar os parâmetros dos exoplanetas. 

O trabalho se constitui em um projeto de Iniciação Científica, apoiado pelo NSEE - Instituto Mauá de Tecnologia e fomentado pela FAPESP - Fundação de Amparo à Pesquisa do Estado de São Paulo, sendo o aluno Guilherme Samuel de Souza Barbosa, 3º ano de Engenharia de Computação pelo IMT, como o bolsista e o Dr. Roberto Bertoldo Menezes (IAG - USP) como orientador do trabalho. 

Em termos gerais, iremos propor uma pipeline de processamento e análise das curvas de luz que se baseia nas etapas de um projeto de ciência de dados, em que a Etapa 1 - Aquisição dos dados será feita realizando o download do conjunto de curvas do SITE DO COROT E SITE DO KEPLER; a Etapa 2 - Processamento e limpeza dos dados constitui em aplicar as diferentes técnicas de filtragem e eleger a melhor para os propósitos desse trabalho; a Etapa 3 - Extração de features dos dados será feita procurando transformações das curvas de luz em prol de alimentar modelos de inteligência a fim de classificá-las e Etapa 4 - Validação do modelo será feita utilizando um outro conjunto de dados a fim de testar a acurácia do modelo treinado.


Introdução aos telescópios espaciais
=====================================

As curvas de luz, apresentadas na sessão anterior, que serão utilizadas nesse projeto, 
foram obtidas por dois telescópios espaciais: CoRoT e o Kepler. 

CoRoT 
-----
`wikipedia-en <https://en.wikipedia.org/wiki/CoRoT>`_

`wikipedia-br <https://pt.wikipedia.org/wiki/Telesc%C3%B3pio_espacial_CoRoT>`_


Kepler
------


Problema proposto
==================
Explicar o problema


Pipeline
=========
Explicar as etapas para resolver o problema


Resultados
===========
Mostrar os resultados obtidos.



.. toctree::
   :maxdepth: 2
   :caption: Manipulating fits files

   01 - Manipulating fits files

.. toctree::
   :maxdepth: 2
   :caption: Filter

   02 - Filters

.. toctree::
   :maxdepth: 2
   :caption: Filtering processes

   03 - Study of filtering processes

.. toctree::
   :maxdepth: 2
   :caption: Simulation

   04 - Simulation   

.. toctree::
   :maxdepth: 2
   :caption: Feature Engineering

   05 - Feature Engineering

.. toctree::
   :maxdepth: 2
   :caption: Machine Learning

   06 - Machine Learning

.. toctree::
   :maxdepth: 2
   :caption: Code Library

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |docs| image:: https://readthedocs.org/projects/filtering-techniques/badge/?version=latest
   :alt: Documentation Status
   :target: https://readthedocs.org/projects/filtering-techniques/

.. |last commit| image:: https://img.shields.io/github/last-commit/Guilherme-SSB/IC-CoRoT_Kepler   
   :alt: GitHub last commit
   :target: https://github.com/Guilherme/IC-CoRoT_Kepler/

.. |most used lang| image:: https://img.shields.io/github/languages/top/Guilherme-SSB/IC-CoRoT_Kepler  
   :alt: GitHub top language