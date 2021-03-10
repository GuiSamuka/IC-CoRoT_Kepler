.. Filtering Techniques and Data Analysis of data from CoRoT and Kepler satellites documentation master file, created by
   sphinx-quickstart on Thu Feb 25 13:40:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Filtering Techniques of Light Curves Obtained from CoRoT and Kepler satellites's documentation!
===========================================================================================================

Curvas de luz são dados astronômicos frequentemente utilizados no estudo da variabilidade do brilho das estrelas e 
consistem essencialmente em gráficos do brilho estelar em função do tempo. Um exame detalhado desses dados permite detectar, 
entre outros fenômenos, a ocorrência de eclipses, que ocorrem quando um corpo, em órbita ao redor da estrela observada, 
passa entre a estrela e o observador, resultando numa diminuição do brilho estelar, que pode ser detectada na curva de luz. 
Esses eclipses podem ser causados por exoplanetas ou mesmo por outras estrelas orbitando a estrela observada. Nesse último caso, 
tem-se as chamadas binárias eclipsantes. A análise de curvas de luz contendo eclipses causados por exoplanetas permite a determinação 
de parâmetros como o raio e a distância de tais exoplanetas em relação à estrela que orbitam. Observações feitas com telescópios 
espaciais forneceram curvas de luz para um número consideravelmente elevado de estrelas, o que resultou na identificação de milhares 
de exoplanetas e de binárias eclipsantes. O telescópio espacial CoRoT (Convection, Rotation and planetary Transits – Auvergne et al. 2009), 
lançado em 27 de dezembro de 2006, por exemplo, que já não se contra mais ativo, produziu curvas de luz para mais de 160000 estrelas. 
Ao final, foram detectados 34 exoplanetas e 2269 binárias eclipsantes (Deleuil et al. 2018). Outro telescópio espacial relevante para esse 
tipo de pesquisa foi o Kepler (Koch et al. 2010), lançado em 7 de março de 2009, que também já não está mais ativo, mas que produziu curvas de 
luz para mais de 200000 estrelas. Essa grande quantidade de dados resultou na detecção de 4034 candidatos a exoplanetas (Thompson et al. 2018), 
dos quais 2335 foram confirmados como exoplanetas, e de mais de 2400 binárias eclipsantes (Slawson et al. 2011; LaCourse et al. 2015). Um problema 
conhecido das curvas de luz está associado ao ruído de alta frequência, que pode dificultar tanto a detecção de eclipses causados por exoplanetas 
ou binárias eclipsantes quanto a modelagem das curvas de luz visando a determinação, por exemplo, de parâmetros dos exoplanetas. Por conta disso, 
métodos de filtragem, visando a remoção do ruído de alta frequência, costumam ser aplicados a curvas de luz. Nesse trabalho, propomos avaliar, usando 
dados disponíveis dos telescópios CoRoT e Kepler, quais os melhores processos de filtragem a serem aplicados a curvas de luz, capazes de remover ruído 
de alta frequência sem comprometer de maneira significativa os dados observados. Os resultados dessa pesquisa deverão resultar em melhorias consideráveis 
na qualidade das curvas de luz obtidas com diferentes instrumentos, o que, por sua vez, proporcionará a obtenção de parâmetros mais precisos a partir das 
análises de tais curvas de luz.

.. toctree::
   :maxdepth: 2
   :caption: Manipulating fits files:

   01 - Manipulating fits files

.. toctree::
   :maxdepth: 2
   :caption: Butterworth Filter:

   02 - Butterworth Filter

.. toctree::
   :maxdepth: 2
   :caption: Butterworth Math Filtering:

   03 - Butterworth filtering without Scipy



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
