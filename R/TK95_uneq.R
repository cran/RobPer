TK95_uneq <-
function(tt, alpha=1.5){
     if(!all(tt==sort(tt))) tt<-sort(tt)
     n<- length(tt)     
     aufl<-10*ceiling(1/ min(apply(matrix(tt[1:(2*(n%/%2))], ncol=2, byrow=TRUE),2, diff))) # Auflösung hängt von kleinstem Zeitabstand zwischen drei Zeitpunkten ab
     N   <-aufl*10                                                  # 10-fache Länge 
     rr_komplett<-TK95(N=N, alpha=alpha)                            # Rotes Rauschen generieren
     stueck <- ceiling(runif(1,min=0, max=0.89)*length(rr_komplett))# Wo wird das Zehntel entnommen?
     ausgesucht<- rr_komplett[stueck:(stueck-1+aufl)]               # Zehntel nehmen
     delta_gl <- (max(tt)-min(tt))/(aufl-1)                         # Abstand zwischen zwei aufeinanderfolgenden Beobachtungen
     stellen<-(((tt-min(tt))+ delta_gl/2)%/%delta_gl)+1             # Welche Stellen man im glm. Samplin min(tt), min(tt)+delta,...,max(tt)-delta, max(tt) ablaufen muss 
     rr<-ausgesucht[stellen]                                         # Entsprechenden Wert aus dem Zehntel verwenden
     return(rr)}
