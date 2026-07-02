liveplot <- function(port = 8789){
library(httpgd, lib.loc = "~/R_Library/4.5")
hgd(host = "0.0.0.0", port = port, token = FALSE)  
# host 0.0.0.0 is important — lets connections in from outside the container
  url <- hgd_url()
  node_name <- sub("^https?://([^:/]+).*", "\\1", url)
message(url)
message("ssh -N -f -L ", port, ":", node_name, ":", port, " pawsey")
message("http://localhost:", port, "/live")
}
