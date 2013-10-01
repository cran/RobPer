TK95 <-
function(N=1000, alpha=1.5){                           # Länge der Zeitreihe, Abstieg des Power law (pink noise: alpha=1)
        f<-seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))] # Fourierfrequenzen
        f_<-1/f^alpha                                            # Power law
        RW<-sqrt(0.5*f_)*rnorm(N/2-1)                            # Realteil der Fouriertrafo 
        IW<-sqrt(0.5*f_)*rnorm(N/2-1)                            # Imaginärteil der Fouriertrafo
        fR<-complex(length.out=N, real=c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), imaginary=c(0, IW, 0, -IW[(N/2-1):1]))
                                                                 # Die komplexen Zahlen, die zurücktransformiert werden sollen,
                                                                 # Reihenfolge der Frequenzen: 0,2pi/N, 2*2pi/N,...,pi,...,2pi-1/N 
                                                                 # So wählen, dass um pi symmetrische Frequenzen konjugiert-komplex
                                                                 # zueinander sind und Frequenzen bei 0 und pi keinen Imaginärteil haben.  
        reihe<-fft(fR, inverse=TRUE)                             # In die Zeit-Domäne überführen
        return(Re(reihe))                                        # Reihe muss nicht länger komplex dargestellt werden, da Imaginärteil 0.
        }
