#--------JESLI KLASY SA W OSTANIEJ KOLUMNIE------------------
#wywolac wszystkie funkcje z pliku Rocket_funkcje_klasa_w_ostatniej_kolumnie
#--------JESLI KLASY SA W PIERWSZEJ KOLUMNIE-----------------
#wywolac wszystkie funkcje z pliku Rocket_funkcje_klasa_w_pierwszej_kolumnie


#--------przyklad-klasy-w-ostatniej-kolumnie-------------------------------------------

Mw<- read.arff(file="ERingDimension1_TRAIN.arff")
Mx <- read.arff(file="ERingDimension2_TRAIN.arff")
My <- read.arff(file="ERingDimension3_TRAIN.arff")
Mz <- read.arff(file="ERingDimension4_TRAIN.arff")
Tw<- read.arff(file="ERingDimension1_TEST.arff")
Tx <- read.arff(file="ERingDimension2_TEST.arff")
Ty <- read.arff(file="ERingDimension3_TEST.arff")
Tz <- read.arff(file="ERingDimension4_TEST.arff")

X <- list(Mw,Mx, My, Mz) #dane treningowe
Y<-list(Tw, Tx, Ty, Tz) #dane testowe

#UWAGA JESLI TARGETY SA NAPISAMI WYBIERZ typ_klasy=1,
#jesli sa liczbami nic nie podawac lub wpisac typ_klasy=0

#-------TESTY------------------------------------
#Na danych pe³nych
dokladnosc_(X,Y,liczba_ker = 100) 

#Na danych ucietych dodatkowo wybór opcji jak uzupelniamy
dokladnosc_uciete(X,Y,liczba_ker = 100,opcja=1) 




#----------TESTY-dla klas w pierwszej kolumnie
#X lista danych treningowych
#Y lista danych testowych
X<-zmiana_klas(X)
Y<-zmiana_klas(Y)
head(X[[1]])
dokladnosc_(X,Y,liczba_ker = 100) 
dokladnosc_uciete(X,Y,liczba_ker = 100,opcja=0) 