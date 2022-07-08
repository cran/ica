ica <- 
  function(X, nc, method = c("fast", "imax", "jade"), ...){
    if(method[1] == "fast"){
      return(icafast(X, nc, ...))
    } else if(method[1] == "imax"){
      return(icaimax(X, nc, ...))
    } else {
      return(icajade(X, nc, ...))
    }
  }