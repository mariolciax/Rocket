################
################ FUNKCJE

# install.packages("tictoc")
# install.packages("foreach")
# install.packages("doParallel")
library(tictoc)
library(foreach)
library(doParallel)
library(foreign)
library(glmnet)

generate_kernels <- function(n_timepoints, num_kernels, n_columns, seed){
  
  if(!missing(seed)){
    set.seed(seed)}
  
  candidate_lengths <- c(7L, 9L, 11L)
  lengths <- sample(x=candidate_lengths, num_kernels, replace = TRUE)
  num_channel_indices <- rep(0L, num_kernels)
  
  for (i in 1:num_kernels) {
    limit <- min(n_columns, lengths[i])
    num_channel_indices[i] <- as.integer(2^(runif(1,0,log2(limit+1))))
  }
  
  channel_indices <- rep(0L, sum(num_channel_indices))
  weights <- rep(0, lengths %*% num_channel_indices)
  biases <- rep(0, num_kernels)
  dilations <- rep(0L, num_kernels)
  paddings <- rep(0, num_kernels)
  a1 <- 1
  a2 <- 1
  
  for (i in 1:num_kernels) {
    length_i <- lengths[i]
    num_channel_indices_i <- num_channel_indices[i]
    weights_i <- rnorm(num_channel_indices_i * length_i, 0, 1)
    b1 <- a1 + (num_channel_indices_i * length_i) - 1
    b2 <- a2 + num_channel_indices_i - 1
    a3 <- 1
   
    for (j in 1:num_channel_indices_i) {
      b3 <- a3 + length_i - 1
      weights_i[a3:b3] <- weights_i[a3:b3] - mean(weights_i[a3:b3])
      a3 <- b3 + 1
    }
    
    weights[a1:b1] <- weights_i
    channel_indices[a2:b2] <- sample(seq(1, n_columns), num_channel_indices_i, replace = FALSE)
    biases[i]<- runif(1,-1,1)
    dilation <- as.integer(2 ^ (runif(1,min = 0, max = log2((n_timepoints)/(length_i)))))
    if(is.na(dilation)){dilation <- 0}
    dilations[i] <- dilation
    padding <- if(sample(c(0,1),1)==1){floor(((length_i) * dilation)/2)} else {0}
    paddings[i] <- padding
    a1 <- b1 + 1
    a2 <- b2 + 1
    
  }
  v_data <- list("weights" = weights, "lengths" = lengths, "biases"=biases, "dilations" = dilations, 
                 "paddings"= paddings, "num_channel_indices" = num_channel_indices, "channel_indices" = channel_indices)
  return(v_data)
}

apply_kernel_univariate <- function(X, weights, length, bias, dilation, padding){
  n_timepoints <- length(X)
  output_length <- n_timepoints + 2 * padding - (length-1) * dilation
  ppv_ <- 0
  max_ <- -Inf
  end <- n_timepoints + padding - (length) * dilation
  
  for (i in -padding:end) {
    sum_ <- bias
    index <- i
    for (j in 1:length) {
      #zmiana -1 na 0 oraz zamiana < na <=
      if(index > 0 && index <= n_timepoints){sum_ <- sum_ + weights[j]*X[index]}
      index <- index + dilation
    }
    if(sum_ > max_){max_ <- sum_}
    if(sum_ > 0){ppv_ <- ppv_ + 1}
  }
  vlist <- c(ppv_/output_length, max_)
  return(vlist)
}

apply_kernel_multivariate <- function(X, weights, length, bias, dilation, padding, num_channel_indices, channel_indices){
  
  n_timepoints <- ncol(X)
  output_length <- n_timepoints + 2 * padding - (length-1) * dilation
  ppv_ <- 0
  max_ <- -Inf
  end <- n_timepoints + padding - (length) * dilation
  
  for (i in -padding:end) {
    sum_ <- bias
    index <- i
    
    for (j in 1:length) {
      
      if(index > 0 && index <= n_timepoints){
        
        for (k in 1:num_channel_indices) 
        {sum_ <- sum_ + weights[k, j] * X[channel_indices[k], index]}
      }
      
      index <- index + dilation
    }
    if(sum_ > max_){max_ <- sum_}
    if(sum_ > 0){ppv_ <- ppv_ + 1}
  }
  vlist <- c(ppv_/output_length, max_)
  return(vlist)
}

X_obs_poj <- function(X, wiersz, n_dim, n_columns) {
  #funckja pomocnicza tworz膮ca macierz zawieraj膮ca wszytskie szeregi czasowe tego 
  #samego eksperymentu o r贸偶nych dimach (zapisane wierszowo, czyli ka偶dy wymiar w innym wierszu)
  #X: macierz szereg贸w czasowych; wiersz: numer potrzebnego wiersza; 
  #n_dim: liczba wymiar贸w; n_timepoints: liczba eksperyment贸w
  
  A <- matrix(NA, nrow = n_dim, ncol = n_columns)
  for (j in 1:n_dim) {
    A[j,] <- as.double(X[[j]][wiersz,])
  }
  return(A)
}

# WOLNA METODA
# apply_kernels <- function(X, kernels){
#   #X to 3d macierz: robimy list(M1, ..., Mn)
#   #kernels: generujemy bior膮c dane z dowolnej macierzy Mi
#   n_instances <- nrow(X[[1]])
#   n_dim <- length(X)
#   n_columns <- ncol(X[[1]])
#   num_kernels <- length(kernels$length)
#   X_ <- matrix(0, ncol = num_kernels * 2 , nrow = n_instances)
#   
#   for (i in 1:n_instances) {
#     a1 <- 1 # for weights
#     a2 <- 1 # for channel_indices
#     a3 <- 1 # for features
#     
#     for (j in 1:num_kernels) {
#       #tic(((i-1)*num_kernels+j)/(num_kernels*n_instances))
#       b1 <- a1 + kernels$num_channel_indices[j] * kernels$lengths[j] - 1
#       b2 <- a2 + kernels$num_channel_indices[j] - 1
#       b3 <- a3 + 2 - 1
#       if(kernels$num_channel_indices[j] == 1){
#         X_[i, a3:b3] = apply_kernel_univariate(as.double(X[[kernels$channel_indices[a2]]][i, ]), kernels$weights[a1:b1], kernels$lengths[j], 
#                                                kernels$biases[j], kernels$dilations[j], kernels$paddings[j])
#       }
#       else {
#         weights_ = matrix(kernels$weights[a1:b1], nrow = kernels$num_channel_indices[j], ncol = kernels$lengths[j])
#         X_[i, a3:b3] = apply_kernel_multivariate(X_obs_poj(X,i,n_dim, n_columns), weights_, kernels$lengths[j], 
#                                                  kernels$biases[j], kernels$dilations[j], kernels$paddings[j], 
#                                                  kernels$num_channel_indices[j], kernels$channel_indices[a2:b2])
#       }
#       a1 <- b1 + 1
#       a2 <- b2 + 1
#       a3 <- b3 + 1
#       #toc()
#     }
#     
#   }
#   
#   
#   return(X_)
#   
# }

#SZYBKA METODA
apply_kernels_fast <- function(X,kernels) {
  library(foreach)
  library(doParallel)
  n_instances <- nrow(X[[1]])
  n_dim <- length(X)
  n_columns <- ncol(X[[1]])
  num_kernels <- length(kernels$length)
  
  myCluster <- makeCluster(detectCores()-1, # number of cores to use
                           type = "PSOCK")
  
  registerDoParallel(myCluster)
  
  output2 <- foreach(i=seq(1,n_instances), .combine = 'rbind', .inorder = FALSE, 
                     .export = c('apply_kernel_univariate', 'apply_kernel_multivariate', 'X_obs_poj')) %dopar% {
                       X_rownolegle <- rep(0, len=2*num_kernels)
                       a1 <- 1 # for weights
                       a2 <- 1 # for channel_indices
                       a3 <- 1 # for features
                       
                       for (j in 1:num_kernels) {
                         #tic(((i-1)*num_kernels+j)/(num_kernels*n_instances))
                         b1 <- a1 + kernels$num_channel_indices[j] * kernels$lengths[j] - 1
                         b2 <- a2 + kernels$num_channel_indices[j] - 1
                         b3 <- a3 + 2 - 1
                         if(kernels$num_channel_indices[j] == 1){
                           X_rownolegle[a3:b3] <- apply_kernel_univariate(as.double(X[[kernels$channel_indices[a2]]][i, ]), kernels$weights[a1:b1], kernels$lengths[j], 
                                                                          kernels$biases[j], kernels$dilations[j], kernels$paddings[j])
                         }
                         else {
                           weights_ = matrix(kernels$weights[a1:b1], nrow = kernels$num_channel_indices[j], ncol = kernels$lengths[j])
                           X_rownolegle[a3:b3] <- apply_kernel_multivariate(X_obs_poj(X,i,n_dim, n_columns), weights_, kernels$lengths[j], 
                                                                            kernels$biases[j], kernels$dilations[j], kernels$paddings[j], 
                                                                            kernels$num_channel_indices[j], kernels$channel_indices[a2:b2])
                         }
                         a1 <- b1 + 1
                         a2 <- b2 + 1
                         a3 <- b3 + 1
                         #toc()
                       }
                       X_rownolegle
                     }
  stopCluster(myCluster)
  return(output2)
}





generate_matrix<-function(X,Y, target_x,target_y, liczba_kerneli, s){
  #X : lista danych treningowych bez targetow, list(Mx, My, Mz)
  #Y : lista danych testowych bez targetow
  #target_x :wektor, target X[1]
  #target_y :wektor, target Y[1]
  #liczba_kerneli: integer
  #s : seed, do generowania
  #output: macierze 1) stworzona z danych treningowych, 2)stworzona z danych testowych
  
  kernels <- generate_kernels(nrow(X[[1]]), liczba_kerneli, length(X), s)
  tic("apply_kernels")
  AK <- apply_kernels_fast(X, kernels)
  toc()
  tic('apply kernel')
  AK_test <- apply_kernels_fast(Y, kernels)
  toc()
  return(list(AK, AK_test))
}


klasyfikacja_czesciowa <-function(AK, AK_test, response_train, response_test){
  #AK-zbior treningowy bez targetow
  #AK_test -zbior testowy bez targetow
  #response_train : targety zbioru treningowego
  #response_test : targety zbioru testowego
  #funkcja tworzy model
  library(glmnet)
  AK<-data.frame(AK)
  AK_test<-data.frame(AK_test)
  
  #dopasowanie modelu ridge regression 
  model <- glmnet(x=AK, y=response_train, alpha = 0)
  #summary(model)
  # k-fold cross-validation zeby znalezc optymalna wartosc lambda 
  cv_model <- cv.glmnet(x=as.matrix(AK), y=response_train, alpha = 0, standarize=T)
  cv_model
  
  #przypisanie optimalnej wartosci lambdy ktora minimalizuje test MSE
  best_lambda <- cv_model$lambda.min
  
  #wykres testu MSE vs lambdy
  #plot(cv_model) 
  
  #tworzymy najlepszy model z dobrana lambda
  best_model <- glmnet(x=as.matrix(AK),y=response_train, alpha = 0,lambda = best_lambda,standarize=T)
  
  
  return(list(best_lambda,best_model))
}

klasyfikacja_<-function(AK, AK_test, response_train, response_test,liczba_klas){
  #budujemy modele dla kazdej klasy osobno
  #AK: macierz wygenerowana przez apply_kernels dla danych treningowych
  #AK_test: macierz wygenerowana przez apply_kernels dla danych testowych
  #response_train: wektor klas zbioru treningowego
  #respone_test: wektor klas zbioru testowego 
  #wynik_klasyfikacji: lista skladajaca sie z dwoch element體, 
  #listy dopasowanych modeli i listy najlepszych lambd dla kazdego modelu
  model_klasa<-c()
  model_lambda<-c()
  
  for (i in 1:liczba_klas){
    a<-response_train
    a[response_train==i]<-1
    a[response_train>i]<- (-1)
    a[response_train<i]<- (-1)
    model_klasa<-c(model_klasa, klasyfikacja_czesciowa(AK, AK_test, a)[2])
    model_lambda<-c(model_lambda,klasyfikacja_czesciowa(AK, AK_test, a)[1] )
  }
  return(list(model_klasa, model_lambda))
  
}

predykcja_<-function( AK_test, liczba_klas, wynik_klasyfikacja){
  #funkcja buduje macierz dla kazdego modelu(liczba modeli =liczbie klas) przewidujac wartosci klas/targetu
  #AK_test: macierz po transformacie przez kernele danych testowych
  #liczba_klas: integer 
  #wynik_klasyfikacji: wektor  z predykcji
  y_predicted<-c()
  for (i in 1:liczba_klas){
    y_predicted<-cbind(y_predicted,predict(wynik_klasyfikacja[[1]][[i]], s =wynik_klasyfikacja[[2]][i], newx = as.matrix(AK_test)))
    #y_predicted <- c(y_predicted, predict(best_model[[i]], s = NULL, newx = as.matrix(AK_test)))
  }
  return(y_predicted)
}



klasyfikacja_ostateczna<-function(wynik_predykcji){
  #wynik predykcji : macierz, gdzie w kazdej kolumnie znajduje sie klasyfikacja 
  #dla kolejnych klas
  #output: wektor klas
  y_predicted<-c()
  wynik_predykcji<-t(wynik_predykcji)
  for (i in 1:ncol(wynik_predykcji)){
    y_predicted<-c(y_predicted,which.max(wynik_predykcji[,i]))
  }
  return(y_predicted)
}


liczba_klas<-function(response_vector){
  #funkcja badajaca liczbe klas/targetow na podstawie listy klas/targetow
  return(length(unique(response_vector)))
}


przygotowanie_klas<-function(X,Y,index_target){
  #X - lista danych treningowych
  #Y -lista danych testowych
  #index_target : numer kolumny w ktorej sa klasy/targety
  id_class<-index_target
  response_train <-X[[1]][,id_class]
  response_train<-as.numeric(as.character(response_train))
  response_test<-Y[[1]][,ncol(Y[[1]])]
  response_test<-as.numeric(as.character(response_test))
  return(list(response_test,response_train))
}

przygotowanie_danych<-function(X,Y, index_target){
  #X - lista danych treningowych
  #Y -lista danych testowych
  #index_target : numer kolumny w ktorej sa klasy/targety
  id_class<-index_target
  #print(ncol(X[[1]]))
  #print(ncol(Y[[1]]))
  for (i in 1:length(X)){
    #X[[i]][,ncol(X[[i]])]<-NULL
    X[[i]]<-X[[i]][,1:ncol(X[[i]])-1]
    #Y[[i]][,ncol(Y[[i]])]<-NULL
    Y[[i]]<-Y[[i]][,1:ncol(Y[[i]])-1]
  }
  return(list(X,Y))
  
}

klasyfikacja<-function(AK,AK_test,response_train,response_test){
  #output: accuracy,macierz_porownania,wektor_najlepszych_lambd
  l_k<-liczba_klas(response_vector = response_test)
  aa<-klasyfikacja_(AK,AK_test,response_train,response_test,l_k)
  best_lambda<-aa[2]
  wynik_predykcja<-predykcja_( as.data.frame(AK_test), l_k, wynik_klasyfikacja=aa)
  przypisanie_klas<-klasyfikacja_ostateczna(wynik_predykcja)
  #porownanie
  macierz_porownanie<-cbind(response_test,przypisanie_klas)
  #accuracy
  accuracy<-sum(macierz_porownanie[,1]==macierz_porownanie[,2])/length(macierz_porownanie[,1])
  return(list(accuracy,macierz_porownanie, best_lambda))
}


wektory_klas_liczby<-function(X,Y, id_target_train=ncol(X[[1]]), id_target_test=ncol(Y[[1]])){
  response_train <-X[[1]][,id_target_train]
  response_train<-as.numeric(as.character(response_train))
  response_test<-Y[[1]][,id_target_test]
  response_test<-as.numeric(as.character(response_test))
  return(list(response_train,response_test))
}

wektory_klas_stringi<-function(X,Y, id_target_train=ncol(X[[1]]), id_target_test=ncol(Y[[1]])){
  response_train <-X[[1]][,id_target_train]
  response_train<-as.numeric(response_train)
  response_test<-Y[[1]][,id_target_test]
  response_test<-as.numeric(response_test)
  return(list(response_train,respone_test))
}

pomocnicza <-function(X,Y,typ_klasy=1, id_target=ncol(X[[1]])){
  #funkcja zamienia typy danych w wektorach klas
  #X - lista danych treningowych
  #Y -lista danych testowych
  id_target_y<-ncol(Y[[1]])
  if (typ_klasy==0){
    for (i in 1:length(X)){
      X[[i]][,id_target]<-as.numeric(X[[i]][,id_target])
      Y[[i]][,ncol(Y[[i]])]<-as.numeric(Y[[i]][,ncol(Y[[i]])])
    }}
  if (typ_klasy==1){
    for (i in 1:length(Y)){
      X[[i]][,id_target]<-as.numeric((X[[i]][,id_target]))
      Y[[i]][,ncol(Y[[i]])]<-as.numeric((Y[[i]][,ncol(Y[[i]])]))
    }}
  return(list(X,Y))
}

badanie_dokladnosci<- function(X,Y, id_target=ncol(X[[1]]), liczba_kerneli=1000,typ=1){
  #X - lista danych treningowych
  #Y -lista danych testowych
  #id_target : numer kolumny w ktorej sa klasy/targety w X
  #response_train: wektor klas z zbioru treningowego
  #response_test:wektor klas z zbioru testowego
  p<-pomocnicza(X,Y, typ_klasy = typ, id_target = id_target)
  X<-p[[1]]
  Y<-p[[2]]
  response_train<-X[[1]][,id_target]
  response_test<-Y[[1]][,ncol(Y[[1]])]
  daane<-przygotowanie_danych(X,Y,index_target = id_target)
  dokladnosc<-c()
  for ( i in 1:10){
    A1<-generate_matrix(X=daane[[1]],Y=daane[[2]],response_train,response_test,liczba_kerneli)
    AK<-A1[1]
    AK_test<-A1[2]
    dokladnosc1<-klasyfikacja(AK,AK_test,response_train,response_test)
    print(dokladnosc1[[1]]) #accuracy
    dokladnosc1[[2]] #macierz porownania
    dokladnosc1[[3]] #lista najlepszych lambd dla kazdego modelu
    
    dokladnosc<-c(dokladnosc,dokladnosc1[[1]])
  }
  return(dokladnosc)
}

dokladnosc_<-function(X,Y, id_target=ncol(X[[1]]),liczba_ker=1000,typ_klasy=1){
  dokladnosc<-badanie_dokladnosci(X,Y,id_target, liczba_kerneli=liczba_ker,typ = typ_klasy)
  #liczymy dane do raportu
  srednia_dokladnosc<-mean(dokladnosc)
  wariancja<-var(dokladnosc)
  dokladnosc<-as.data.frame(dokladnosc)
  dokladnosc_macierz <-rbind(dokladnosc,srednia_dokladnosc)
  dokladnosc_macierz <-rbind(dokladnosc_macierz,wariancja)
  DD <-as.data.frame(dokladnosc_macierz, col.names<- "dokladnosc",row.names = c("doswiadczenie 1","doswiadczenie 2",
                                                                                "doswiadczenie 3","doswiadczenie 4",
                                                                                "doswiadczenie 5","doswiadczenie 6",
                                                                                "doswiadczenie 7","doswiadczenie 8",
                                                                                "doswiadczenie 9","doswiadczenie 10", 
                                                                                "srednia dokladnosc", "wariancja"
  ))
  write.table(DD,file="dokladnosc_zbior.csv",row.names=T)
  return(DD)
}




dokladnosc_uciete<-function(X,Y,wartosc=0,opcja=0,liczba_kerneli=1000,typ_klasy= 1){
  X<-change_data_frames(X)
  Y<-change_data_frames(Y)
  X<-wypelnij_data_frames(X,opcja=opcja,wartosc=wartosc)
  Y<-wypelnij_data_frames(Y,opcja=opcja,wartosc=wartosc)
  dokladnosc_uciete<-dokladnosc_(X,Y,liczba_ker=liczba_kerneli,typ = typ_klasy)
  return(dokladnosc_uciete)
}



change_data_frames<-function(list_data_frames){
  #list_data_frames as list f.ex list(df1,df2)
  df1 <-list_data_frames[[1]]
  class_column <- df1[ , ncol(df1)]
  class_column_name <-names(df1)[ncol(df1)]
  list_classes_df1 <-unique(class_column)
  
  for (i in 1:length(t(list_classes_df1))){
    row_ids <-c()
    class <- t(list_classes_df1)[i]
    for (j in 1: nrow(df1)){
      df1_j<-df1[ j, ncol(df1)]
      
      if (df1_j==class){
        row_ids<-c(row_ids, j)
      }
    }
    list_data_frames<-change_rows(list_data_frames,row_ids)
  }
  return(list_data_frames)
  
}


#wypelnij_data_frames<-function(list_data_frames, wartosc=0){
#list_data_frames as list f.ex list(df1,df2)
#wartosc: wartosc jaka chcemy wypelnic wartos NA
#  l<-length(list_data_frames)
#  for (i in 1:l){
#    list_data_frames[[i]][is.na(list_data_frames[[i]])]=wartosc
#  }
#  return(list_data_frames)
#}

wypelnij_data_frames<-function(list_data_frames, opcja=0, wartosc=0){
  #list_data_frames as list f.ex list(df1,df2)
  #opcja: wybor opcji jak postepujemy z poucinanymi danymi
  #wartosc: (jesli wybrano opcje 0) - wartosc jaka chcemy wypelnic wartos NA 
  #opcja 0: wypelniamy domyslnie zerami, mozna zmienic wartos domyslna
  #opcja 1: wypelnienie srednia z wiersz
  #opcja 2: skrocenie do najkrotszego wiersza
  if (opcja==0){
    list_data_frames<-wypelnij_data_frames0(list_data_frames, wartosc = wartosc)
  }
  if (opcja==1){
    list_data_frames<-wypelnij_data_frames1(list_data_frames)
  }
  if (opcja==2){
    list_data_frames<-wypelnij_data_frames2(list_data_frames)
  }
  
  
  return(list_data_frames)
}


wypelnij_data_frames0<-function(list_data_frames, wartosc=0){
  #list_data_frames as list f.ex list(df1,df2)
  #wartosc: wartosc jaka chcemy wypelnic wartos NA
  l<-length(list_data_frames)
  for (i in 1:l){
    list_data_frames[[i]][is.na(list_data_frames[[i]])]=wartosc
  }
  return(list_data_frames)
}
wypelnij_data_frames1<-function(list_data_frames){
  #wypelnienie srednia z wiersza
  #list_data_frames as list f.ex list(df1,df2)
  
  l<-length(list_data_frames)
  for (i in 1:l){
    liczba_wierszy<-nrow(list_data_frames[[i]])
    liczba_kolumn<-ncol(list_data_frames[[i]])
    df_i<-list_data_frames[[i]]
    for(j in 1:liczba_wierszy){
      wartosc<-mean(as.numeric(df_i[j,1:(liczba_kolumn-1)]),na.rm=TRUE)
      #print(wartosc)
      list_data_frames[[i]][j,][is.na(list_data_frames[[i]][j,])]=wartosc
    }
  }
  return(list_data_frames)
}

  
  wypelnij_data_frames2<-function(list_data_frames){
    #list_data_frames as list f.ex list(df1,df2)
    #uciecie wszystkich do najkrotszego wiersza
    l<-length(list_data_frames)
    koniec<-0
    lista<-c()
    liczba_kolumn<-ncol(list_data_frames[[1]])
    #print(liczba_kolumn)
    liczba_wierszy<-nrow(list_data_frames[[1]])
    for (i in 1:l){
      df_i<-list_data_frames[[i]]
      for (j in 1:liczba_wierszy){
        koniec<-max(koniec, sum(is.na(as.numeric(df_i[j,1:(liczba_kolumn-1)]))))
      }
    }
    for (i in 1:l){
      df_i<-list_data_frames[[i]]
      klasa<-df_i[,ncol(list_data_frames[[i]])]
      endd<-liczba_kolumn-koniec-1
      lista[[i]]<-cbind(list_data_frames[[i]][,1:endd],klasa)
    }
    #print(head(lista[[1]]))
    return(lista)
  }

change_rows<-function(dimensions, row_ids){
  #dimensions as list of data frames
  #row_ids as a vector of row id in the class, f.ex c(1,10,12) (this rows are in first class)
  l<-length(row_ids)
  if (l%%3==0){
    l1<-l/3
    l2<-l1+l1
    l3<-l1+l2
  }
  
  if(l%%3==1){
    l1<-(l-1)/3
    l2<-l1+l1
    l3<-l1+l1+l1+1
  }
  
  if(l%%3==2){
    l1<-(l-2)/3
    l2<-l1+l1+1
    l3<-l1+l2+1
  }
  
  for (i in 1:l1){
    percent1 <-runif(1, min=0.1, max=0.4)
    row_id<-row_ids[i]
    dimensions<-make_Na(dimensions, row_id ,percent1)
  }
  for (i in (l1+1):l2){
    percent2 <-runif(1, min=0.4, max=0.7)
    row_id<-row_ids[i]
    dimensions<-make_Na(dimensions, row_id ,percent2)
  }
  for (i in (l2+1):l3){
    percent3 <- runif(1, min=0.7, max=1)
    row_id<-row_ids[i]
    dimensions<-make_Na(dimensions, row_id ,percent3)
  }
  return(dimensions)
}


make_Na <- function(dimensions, row_id, percent){
  #dimension as list of data frames
  #row_id as vector of row id in the class, f.ex c(1,10,12) (this rows are in first class)
  #percet as double between 0 and 1
  list_dims<-list()
  l<-length(dimensions)
  for (j in 1:l) {
    x<-j
    dim_j<- data.frame(dimensions[x],stringsAsFactors = FALSE) 
    n <- ncol(dim_j)
    lower<-ceiling(percent*(n-1))+1
    upper<-n-1
    if (lower==upper){
      dim_j[row_id,lower]<-NA
      dim_j<-data.frame(dim_j, stringsAsFactors = FALSE)
    } 
    if (lower<upper){
      dim_j[row_id,lower:upper]<-NA
      dim_j<-data.frame(dim_j, stringsAsFactors = FALSE)
    }
    if(lower>upper){
      dim_j[row_id,upper]<-NA
      dim_j<-data.frame(dim_j, stringsAsFactors = FALSE)
    }
    
    list_dims[[x]] <-dim_j
    
  }
  return(list_dims)}


zmiana_kolumny<-function(macierz){
  #macierz: podajemy macierz w kt髍ej chcemy zmienic miejsce kolumny klas
  klasy<-macierz[,1]
  macierz<-cbind(macierz[,2:ncol(macierz)], klasy)
  return(macierz)
  
}

zmiana_klas<-function(lista_data_frames){
  #lista_data_frames : lista macierzy
  len<-length(lista_data_frames)
  
  for (i in 1:len){
    lista_data_frames[[i]]<-zmiana_kolumny(lista_data_frames[[i]])
  }
  return(lista_data_frames)
}


