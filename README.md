NenuRaw
=======

NenuRaw est une bibliothèque logicielle qui permet de lire et de traiter des fichiers de données brutes provenant d'observations radioastronomiques. Elle fournit des fonctionnalités pour analyser et
visualiser les données, ainsi que pour effectuer des opérations de traitement telles que le nettoyage des données dans le domaine de fourier, la dedispertion cohérente, l'integration temporel et frequenciel, etc.

Les fichiers de la bibliothèque NenuRaw sont organisés en plusieurs classes, chacune ayant une fonction spécifique. Voici une description des fichiers et des classes les plus importants :

* raw_utils.py : Ce fichier contient la classe **Raw** qui est la classe principale de la bibliothèque. Elle est utilisée pour lire et traiter les fichiers de données brutes. Elle contient des méthodes pour extraire les informations des en-têtes des fichiers, pour effectuer des opérations de traitement sur les données, et pour générer des visualisations.
* main.py : Ce fichier est un exemple d'utilisation de la bibliothèque NenuRaw. Il montre comment créer une instance de la classe **Raw** et comment utiliser ses méthodes pour lire et traiter les fichiers de données brutes.
* wav_utils.py : Ce fichier contient la classe Wav qui est utilisée pour effectuer des opérations de traitement spécifiques sur les données brutes, telles que la démodulation cohérente et le nettoyage des fréquences.
* dynspec_utils.py : Ce fichier contient la classe **Dynspec** qui hérite de la classe **Raw** et ajoute des fonctionnalités spécifiques pour l'analyse et la visualisation des données spectroscopiques.
* rawtf_structure.py : Ce fichier contient la classe Rawtf qui est utilisée pour lire et traiter les fichiers de données brutes au format RAWTF.
* luppi_structure.py : Ce fichier contient la classe LUPPI qui est utilisée pour lire et traiter les fichiers de données brutes au format GUPPI.
* waveolaf_structure.py : Ce fichier contient la classe WaveOlaf qui est utilisée pour lire et traiter les fichiers de données brutes issue du flux des données haut debit du style LOFAR.

Ces fichiers et classes constituent les principaux composants de la bibliothèque NenuRaw et permettent de lire, traiter et analyser les données brutes provenant d'observations radioastronomiques.

Classe Raw
----------

La classe **Raw** dans le fichier raw_utils.py est utilisée pour lire et manipuler des fichiers de données brutes provenant de l'observatoire de Nançay. Voici une description des principales méthodes de la
classe **Raw** :

* __init__(self, obs_files, freq_start=0, freq_end=999, start=None, end=None, duration=None, block_start=0, block_end=1, verbose=False): Le constructeur de la classe **Raw** initialise les attributs de l'objet en fonction des paramètres fournis. Il prend en entrée une liste de fichiers d'observation, ainsi que des paramètres optionnels tels que la fréquence de début et de fin, le temps de début et de fin, la durée, le bloc de début et de fin, et un indicateur de verbosité.
* __try_isot(self, time_string): Cette méthode tente de convertir une chaîne de caractères représentant une date/heure au format ISO en un temps en secondes depuis le début de l'observation. Elle renvoie le temps en secondes si la conversion réussit, sinon elle renvoie False.
* __try_iso(self, time_string): Cette méthode tente de convertir une chaîne de caractères représentant une date/heure au format ISO ou ISO étendu en un temps en secondes depuis le début de l'observation. Elle renvoie le temps en secondes si la conversion réussit, sinon elle renvoie False.
* __nfft(self): Cette méthode calcule le nombre de transformées de Fourier à effectuer en fonction du nombre de canaux, de la longueur de la FFT et du nombre de polarisations.
* __ds(self): Cette méthode calcule la taille du sous-échantillonnage en fonction de la taille de la FFT, de la durée du sous-échantillonnage et de la largeur de bande des canaux.
* __df(self): Cette méthode arrondit la largeur de bande des canaux à une puissance de 2.
* fourier_computation(self, fftlen=16, ds_ms=30, df=1, file_num=0, pol="I"): Cette méthode effectue le calcul de la transformée de Fourier sur les données brutes. Elle prend en entrée les paramètres de la FFT (longueur de la FFT, durée du sous-échantillonnage, largeur de bande des canaux), ainsi que le numéro du fichier et la polarisation à utiliser. Elle renvoie les données transformées.
* fourier_computation_block(self, block_num, file_num=0, pol='I', dm=None): Cette méthode effectue le calcul de la transformée de Fourier sur un bloc spécifique des données brutes. Elle prend en entrée le numéro du bloc, le numéro du fichier, la polarisation et le DM (Dispersion Measure) à utiliser. Elle renvoie les données transformées pour le bloc spécifié.
* get_patidx(self, data, file_num, block_num): Cette méthode renvoie l'indice de paquet (packet index) pour un bloc spécifique des données brutes.
* get_block(self, data, file_num, block_num): Cette méthode renvoie un bloc spécifique des données brutes.

Ces méthodes permettent de lire et de manipuler les données brutes, d'effectuer des calculs de transformée de Fourier et d'obtenir des informations sur les blocs de données.


Classe Dynspec
--------------

La classe **Dynspec** dans le fichier dynspec_utils.py est une classe dérivée de la classe **Raw** de la bibliothèque NenuRaw. Elle ajoute des fonctionnalités spécifiques pour l'analyse et la visualisation des
données spectroscopiques. Voici une description des méthodes associées à la classe **Dynspec** :

* __init__(self, *args, **kwargs): Le constructeur de la classe **Dynspec** qui appelle le constructeur de la classe **Raw** et initialise les paramètres spécifiques à la classe Dynspec.
* dedisperse(self, dm=None): Une méthode qui effectue la dédispersion des données en appliquant un décalage temporel en fonction de la dispersion du signal. Elle prend en paramètre optionnel dm qui représente la dispersion en pc cm^-3.
* clean(self, threshold=3): Une méthode qui effectue le nettoyage des données en supprimant les valeurs aberrantes ou les pics indésirables. Elle prend en paramètre optionnel threshold qui représente le seuil de nettoyage.
* rm_dm0(self): Une méthode qui supprime la dispersion à zéro en ajustant les valeurs temporelles des données.
* rm_baseline(self): Une méthode qui supprime la ligne de base des données en ajustant les valeurs fréquentielles des données.
* rm_pfb(self): Une méthode qui supprime la réponse du filtre passe-bande en ajustant les valeurs fréquentielles des données.
* plot_dynspec(self): Une méthode qui génère une visualisation de la dynamique spectrale des données.
* plot_spectra(self): Une méthode qui génère une visualisation des spectres de puissance des données.
* print_info(self): Une méthode qui affiche les informations sur les données.
* get_spectra(self): Une méthode qui retourne les spectres de puissance des données.
* get_data(self): Une méthode qui retourne les données.
* get_freqarray(self): Une méthode qui retourne un tableau des fréquences.
* get_timearray(self): Une méthode qui retourne un tableau des temps.
* get_start(self): Une méthode qui retourne le temps de début des données.

Les méthodes permettent d'effectuer différentes opérations sur les données spectroscopiques, telles que la dédispersion, le nettoyage, la suppression de la ligne de base, etc. Elles fournissent également
des informations sur les données et permettent d'accéder aux spectres de puissance et aux tableaux de fréquences et de temps.

