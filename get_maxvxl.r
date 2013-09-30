get_maxvxl = function(flnm, onm, thr, tp='n'){

    dt = read.table(flnm, header=FALSE, col.names=c('X','Y','Z','val'))

    sr=sort(dt$val,index.return = TRUE)

    if(tp=='r'){
        thr = round((thr/100) * length(dt$val))
    }

    useIX=sr$ix[1:thr]

    ot = data.frame(dt$X[useIX], dt$Y[useIX], dt$Z[useIX])

    write.table(ot, file=onm,row.names=FALSE, col.names=FALSE)
}